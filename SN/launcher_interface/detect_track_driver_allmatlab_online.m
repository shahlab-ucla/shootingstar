%runs nuclear detection 

%input
%parameterfile is full path location of parameter file
% embryodir is full path of directory containing image files e.g.
% L:/disk1/bao/083105/image/tif/
% embryonumber is the prefix of the image file name e.g. 083105_L1
%suffix is an optional string that will be appended to the name to all
%output keeping one processing of a set of image data separate from another
% for example to test multiple parameter sets, or (in conjunction with ROI parameters) 
%if multiple embryos are in one image 
%outputdirectory is where the output goes

%output
%outputs int outputdir/embryonumber_fullmatlabresult.mat a dump of matlab
%detection results. The important part of this is the cell array esequence
%which contains structure fields finalpoints, finaldiams, which
%are the position and size of nuclei.
%output int outputdir is a zip file contaning lineaged results

%test params to run from commandline as script
%parameterfile='C:\Users\shahp2\Documents\MATLAB\SN\launcher_interface\Red_anthonytest_onepfivemicronz.txt';embryodir='L:\kovacevi\Si_obj_trial\20140407_JIM113_SiO60\';embryonumber = '20140407_JIM113_SiO-0.15_1_s1';suffix='_output';outputdirectory='./testoutput/';polygon_points=[];isred=1;

function detect_track_driver_allmatlab_online(parameterfile,embryodir,embryonumber,suffix,outputdirectory,polygon_points,isred,TargetCell)
%Import AceTree classes
import org.rhwlab.acetree.AceTreeLauncher
import org.rhwlab.acetree.AceTree
import org.rhwlab.snight.NucleiMgr

runexpression=true;
newimage=true;

'beginning lineaging'
tic

TargetCell

%cant pass blank string parameter on command line, so comes in as undef
if(~exist('suffix'))
    suffix='';
end


if(~exist('nodata'))
    nodata=true;%whether to use SN data to match online
end
nodatause=true;%whether to use SN data for diameter and stored bottom
savedata=true;
singlevolume=false;


readParameters;

load(distribution_file);

anisotropy=zres/xyres*downsampling;
downsample=downsampling;
%tlist=linspace(1,run_time,(run_time));
firsttimestep=1;
allvalid=[];



eall={};

bottomdata={};

if(singlevolume&~nodata)
    zlevel=embryolevel;
else
    zlevel=slices;
end
%nucleibase=[nucleidir,embryonumber,'\'];

%Generate output xml file
mkdir([outputdirectory,suffix,embryonumber,'/nuclei']);

zipname = [outputdirectory,embryonumber,'_',suffix,'.zip'];
xmlname=[outputdirectory,embryonumber,'_',suffix,'.xml'];
file=fopen(xmlname,'w');

fprintf (file, '<?xml version=''1.0'' encoding=\''utf-8\''?>\n');
fprintf (file, '<embryo>\n');

if (newimage)
    %For 16-bit tiff usage
    %fprintf (file,'<useStack type="1"/>\n');
    %fprintf (file,['<image file="',embryodir,embryonumber,'_t',num2str(start_time),'.tif"/>\n']);
    fprintf (file,['<image file="',embryodir,'image/tif/',embryonumber,'-t',num2str(start_time,'%03d'),'-p01.tif"/>\n']);
else
    fprintf (file,['<image file="',embryodir,'image/tif/',embryonumber,'-t',num2str(start_time,'%03d'),'-p01.tif"/>\n']);
end
fprintf (file,['<nuclei file="',zipname,'"/>\n']);

fprintf (file,'<end index="475"/>\n');
fprintf(file,['<resolution xyRes="',num2str(xyres),'" zRes="',num2str(zres),'" planeEnd="',num2str(slices),'"/> <exprCorr type="blot"/>\n']); 
fprintf (file,'<polar size="15"/>\n');
fprintf (file,'</embryo>');
fclose(file);

master_esequence = {};

saveGreedyNucleiFiles(master_esequence,0,[outputdirectory,suffix,embryonumber,'/nuclei'],anisotropy,ROIxmin,ROIymin);
zip(zipname,[outputdirectory,suffix,embryonumber]);

%Fire up AceTree with the config file that was just generated
AceTree.main(xmlname)
pause(10)
%Get handles to the AceTree and NucleiManager instances that AT has created
%Java objects pulled from NucleiManager are not copies but the instances
%themselves and edits will be reflected to NucleiManager immediately.
AT = AceTree.getAceTree([]);
NM = AT.getNucleiMgr;

run_time = end_time;
persistance=0;

for init = 1:run_time %Force pipeline to run one image at a time
    ttotal = tic;
    %Check to see if the master record has already been initialized, if
    %yes, read in NM data and merge into the master record if any changes
    %have been made in AT
    
    if init == 1
        firsttimestep = true;
    end
    esequence={};
    tlist = init;
    %pause(30) %Process sequence now includes image waiter, comment this out when running live
    start_time = init;
    end_time = init;
    tseg = tic;
    processSequence; %start segmentation for the timepoint as defined by init
    segmentation_time(init) = toc(tseg);
    if(~isempty(esequence{init}.finalpoints))
        %note off by 2 correction here ROI is first included pixel so
        %pos=pos+roi-1  Acetree coordinate system is 0 origin rather than 1
        %which is the origin of second subtraction
        esequence{init}.finalpoints(:,1)=  esequence{init}.finalpoints(:,1)+ROIxmin-2;
        esequence{init}.finalpoints(:,2)=  esequence{init}.finalpoints(:,2)+ROIymin-2;
    end
    try 
        a = esequence{init}.pred;
    catch
        for j = 1:length(esequence{init}.diams)
            esequence{init}.pred(j) = -1;
            esequence{init}.suc(j,1) = -1;
            esequence{init}.suc(j,2) = -1;
        end
    end
    master_esequence{init} = esequence{init}; %add the results of segmentation into the appropriate position of the master structure
    
    if exist('master_esequence','var')
        while true
            %attempt to lock NM
            success = AT.SNLockNucleiMgr(true);
            if success
                %if successful (ie. not already locked by AT) break and
                %proceed
                break
            else
                %otherwise, pause for 100ms and loop around
                pause(0.1)
            end
        end
        if init > 1    
            master_esequence = DiffAndMergeNuclei(master_esequence, NM);
        end
    end
    
    if init>1
        'detection completed, beginning tracking'
        skipbifurcation = 1; %for pairwise analysis, ignore bifurcation analysis
        trackingparameters.trainingmode=false;
        trackingparameters.recordanswers=false;
        evalforced=false;
        trackingparameters.starttime=1;
        trackingparameters.endtime=2; %since we're only going to pass two timepoints, the endtime is 2
        trackingparameters.anisotropyvector=[1,1,anisotropy];
        parameters.anisotropyvector=[1,1,anisotropy];
        esequence = cell(1,2);
        esequence{1} = master_esequence{init-1}; %pull out the last two timepoints to pass on for tracking
        esequence{2} = master_esequence{init};
        ttracking = tic;
        tracking_driver_new_classifier_based_version
        tracking_time(init) = toc(ttracking);
        [m,~] = size(esequence{1}.suc_time);
        %since we process pairs here, the sucessor only is calculated for
        %init - 1 while the predecessor can only be determined for init and
        %if one of either exists, we know it's either the previous time or
        %the current time
        for i = 1:m
            if esequence{1}.suc_time(i,1) ~= -1
                esequence{1}.suc_time(i,1) = init; %since we're only looking at pairs, the sucessor, if it exists is always where we are
            end
            if esequence{1}.suc_time(i,2) ~= -1
                esequence{1}.suc_time(i,2) = init;
            end
        end
        m = length(esequence{2}.pred_time);
        for i = 1:m
            if esequence{2}.pred_time(i) ~= -1
                esequence{2}.pred_time(i) = init-1; %since this chain only processes pairs, the predecessor will always be only 1 timepoint prior
            end
        end
        master_esequence{init-1} = esequence{1}; %merge the results back into the master structure
        master_esequence{init} = esequence{2};
    end
    
    if (init > 10)% && ~rem(init,5)) || init == run_time
        skipbifurcation = 0;
        'starting windowed tracking'
        trackingparameters.trainingmode=false;
        trackingparameters.recordanswers=false;
        evalforced=false;
        trackingparameters.starttime=1;
        trackingparameters.endtime=10; %since we're only going to pass 10 timepoints, the endtime is 10
        trackingparameters.anisotropyvector=[1,1,anisotropy];
        parameters.anisotropyvector=[1,1,anisotropy];
        parfor k = init-10:init
            esequence{k} = master_esequence{k};
        end
        ttracking = tic;
        tracking_driver_new_classifier_based_version
        tracking_time(init) = toc(ttracking);
        %To adjust time here, we just slide the scale based on where we are.
        %We don't have to worry about adjustments that were made
        %sequentially since those just get overwritten with the current
        %1:10 timescale
        parfor j = 1:10
            [m,~] = size(esequence{j}.suc_time);
            for i = 1:m
                if esequence{j}.suc_time(i,1) ~= -1
                    esequence{j}.suc_time(i,1) = esequence{j}.suc_time(i,1) + (init - 10);
                end
                if esequence{j}.suc_time(i,2) ~= -1
                    esequence{j}.suc_time(i,2) = esequence{j}.suc_time(i,2) + (init - 10);
                end
            end
        end
        parfor j = 1:10
            m = length(esequence{j}.pred_time);
            for i = 1:m
                if esequence{j}.pred_time(i) ~= -1
                    esequence{j}.pred_time(i) = esequence{j}.pred_time(i) + (init - 10);
                end
            end
        end
        
        %The window breaks the pred / suc link between the first timepoint
        %in the window and the last timepoint before the window, re-link 
        %them here
        
        if init>11;
            [m,n] = size(master_esequence{init-10}.suc);
            for i = 1:m
                for j = 1:n
                    if master_esequence{init-10}.suc(i,j) ~= -1
                        ind = master_esequence{init-10}.suc(i,j);
                        esequence{1}.pred(ind) = i;
                        esequence{1}.pred_time(ind) = init-10;
                    end
                end
            end
        end
        
        for i = 1:10
            master_esequence{init - (10-i)} = esequence{i};
        end
        

    end
    
    %Merge the new tracking data into NM
    if init>1
        master_esequence = DiffAndMergeNuclei(master_esequence, NM);
        %Rebuild the tree here since this is the only place where data can flow
        %from SN -> AT
        AT.clearTree();
        AT.buildTree(true);
        AT.buildTree(true); %sometimes acetree loses names, rebuild twice to fix this
        %Now that NM has been updated, unlock it
    end
    

    if init>10
        fid = fopen('C:\Users\shahp2\Desktop\MicroPointTarget\TargetCell.txt','r');
        if fid ~= -1
            TargetCell = fscanf(fid,'%s');
            if ~strcmp(TargetCell,'none')
                fprintf(1,'Target cell defined\n');
                TargetCell
                [Exists,Pos] = TargetCellExists(NM,TargetCell,init);
                if Exists
                    persistance = persistance + 1;
                    if persistance == 5 || persistance == 7
                        fprintf(1,'Targeted cell found\n');
                        fprintf(1,'Persisted for %f timepoints',persistance);
                        fprintf(1,'Targeted cell found, starting killing');
                        fid = fopen('C:\Users\shahp2\Desktop\MicroPointTarget\coords.ini','w');
                        fprintf(fid,strcat('[main]\r\nstatus=Ready\r\nX=',num2str(512-Pos(1)),'\r\nY=',num2str(Pos(2)),'\r\nZ=',num2str(Pos(3)),'\r\n'));
                        fclose(fid);
                    end
                else
                    fprintf(1,'Targeted cell not found\n');
                end
            end
        else
            fprintf(1,'Target cell not identified\n');
        end
        fclose(fid);
    end


    success = AT.SNLockNucleiMgr(false);
    total_time(init) = toc(ttotal);

end

esequence = master_esequence;
save('times.mat','esequence','time_to_read','time_to_link_hard','total_time','tracking_time','segmentation_time','time_to_filter','time_to_castrays','time_to_merge_nuclei','time_to_link_easy');
%Write out a backup
%saveGreedyNucleiFiles(master_esequence,init,[outputdirectory,suffix,embryonumber,'/nuclei'],anisotropy,ROIxmin,ROIymin);
%zip(zipname,[outputdirectory,suffix,embryonumber,'/nuclei']);
%rmdir([outputdirectory,suffix,embryonumber],'s');

'lineaging completed'