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

function detect_track_driver_allmatlab_REALTIME(parameterfile,embryodir,embryonumber,suffix,outputdirectory,polygon_points,isred,TargetCell)
%Import AceTree classes
import org.rhwlab.acetree.AceTreeLauncher
import org.rhwlab.acetree.AceTree
import org.rhwlab.snight.NucleiMgr

runexpression=true;
newimage=true;

'beginning lineaging'
tic

TargetCell = '';

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
    %old AT format fprintf (file,['<image file="',embryodir,'image/tif/',embryonumber,'-t',num2str(start_time,'%03d'),'-p01.tif"/>\n']);
    fprintf (file,['<image file="',embryodir,embryonumber,'_t',num2str(start_time),'.tif"/>\n']);
else
    %old AT format
    fprintf (file,['<image file="',embryodir,embryonumber,'_t',num2str(start_time),'.tif"/>\n']);
end
fprintf (file,['<nuclei file="',zipname,'"/>\n']);

fprintf (file,'<end index="475"/>\n');
fprintf(file,['<resolution xyRes="',num2str(xyres),'" zRes="',num2str(zres),'" planeEnd="',num2str(slices),'"/> <exprCorr type="blot"/>\n']); 
fprintf (file,'<polar size="15"/>\n');
fprintf (file,'<useStack type="1"/>\n');
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

for init = 1:run_time %Force pipeline to run one image at a time
    %Check to see if the master record has already been initialized, if
    %yes, read in NM data and merge into the master record if any changes
    %have been made in AT
    
    if init == 1
        firsttimestep = true;
    end
    esequence={};
    tlist = init;
    %pause(30) %Process sequence now includes image waiter, comment this out when running live
    start_time = 1;
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
 
    master_esequence{init} = esequence{init}; %add the results of segmentation into the appropriate position of the master structure
    
    if init>1
        'detection completed, beginning tracking'
        if init>6
            skipbifurcation = 0; 
        else
            skipbifurcation = 1;%for pairwise analysis, ignore bifurcation analysis
        end
        trackingparameters.trainingmode=false;
        trackingparameters.recordanswers=false;
        evalforced=false;
        trackingparameters.starttime=1;
        trackingparameters.endtime=init; %since we're only going to pass two timepoints, the endtime is 2
        trackingparameters.anisotropyvector=[1,1,anisotropy];
        parameters.anisotropyvector=[1,1,anisotropy];
        esequence = master_esequence;
        ttracking = tic;
        tracking_driver_new_classifier_based_version_REALTIME
        master_esequence = esequence;
        tracking_time(init) = toc(ttracking);
    end
    
    %Merge the new tracking data into NM
    if exist('master_esequence','var')
        while true
            %attempt to lock NM
            success = AT.SNLockNucleiMgr(true);
            if success
                %if successful (ie. not already locked by AT) break and
                %proceed
                break
            else
                %otherwise, pause for 10ms and loop around
                pause(0.01)
            end
        end
        if init>1
            mergeT = tic;
            master_esequence = DiffAndMergeNuclei(master_esequence, NM);
            TimeSpentMerging(init) = toc(mergeT);
            %Rebuild the tree here since this is the only place where data can flow
            %from SN -> AT
            try
            AT.clearTree();
            AT.buildTree(true);
            catch err
            end
            try
            AT.buildTree(true);
            catch err
            end
        end
    end
    
    if init>5
        fid = fopen('C:\Users\shahp2\Desktop\MicroPointTarget\TargetCell.txt','r');
        if fid ~= -1
            text = fscanf(fid,'%s');
            fclose(fid);
            TargetCells = strsplit(text,',');
            NumTargets = length(TargetCells);
            if ~exist('persistance','var') || (length(persistance)~=NumTargets)
                persistance = zeros(1,NumTargets);
                killed = zeros(1,NumTargets);
                %thresh sets the number of timepoints before ablating
                thresh = 7.*ones(1,NumTargets);
            end
            used = 0;
            for k = 1:NumTargets
                if ~strcmp(TargetCell,'none')
                    TargetCell = TargetCells{k};
                    fprintf(1,strcat('\nTarget cell defined\n',TargetCell,'\n'));
                    check = 1;
                    while check
                        try
                            [Exists,Pos] = TargetCellExists(NM,TargetCell,init);
                            check = 0;
                        catch err
                            pause(0.5);
                        end
                    end
                    if Exists
                        fid2 = fopen('C:\Users\shahp2\Desktop\MicroPointTarget\log.txt','a');
                        persistance(k) = persistance(k) + 1;
                        fprintf(1,strcat('\nTargeted cell found\n',TargetCell,'\n'));
                        fprintf(1,'Persisted for %1.1f timepoints\r\n',persistance(k));
                        fprintf(fid2,strcat('\nTimepoint = \n',num2str(init),'\n'));
                        fprintf(fid2,strcat('\nTargeted cell found\n',TargetCell,'\n'));
                        fprintf(fid2,'Persisted for %1.1f timepoints\r\n',persistance(k));
                        %thresh is the number of timepoints in a row a cell
                        %has to be seen to trigger killing, killed is the
                        %number of times to ablate a cell. If ablating many
                        %times, the if statement sets the gap before the
                        %cell is ablated again
                        if persistance(k) >= thresh(k) && killed(k) < 1  && used == 0
                            fprintf(1,'Targeted cell found, starting killing\r\n');
                            fprintf(fid2,'Targeted cell found, starting killing\r\n');
                            killed(k)=killed(k)+1;
                            if killed(k) == 1
                                thresh(k) = thresh(k) + 5;
                            end
                            fid = fopen('C:\Users\shahp2\Desktop\MicroPointTarget\coords.ini','w');
                            fprintf(fid,strcat('[main]\r\nstatus=Ready\r\nX=',num2str(512-Pos(1)),'\r\nY=',num2str(Pos(2)+3),'\r\nZ=',num2str(Pos(3)),'\r\n'));
                            fprintf(fid2,strcat('[main]\r\nstatus=Ready\r\nX=',num2str(512-Pos(1)),'\r\nY=',num2str(Pos(2)+3),'\r\nZ=',num2str(Pos(3)),'\r\n'));
                            fclose(fid);
                            fclose(fid2);
                            used=1;
                        end
                    else
                        fprintf(1,'Targeted cell not found\r\n');
                    end
                end
            end
        else
            fprintf(1,'Target cell not identified\n');
        end
    end

    success = AT.SNLockNucleiMgr(false);
    total_time(init) = toc(ttotal);
    fprintf(1,'Total time elapse = %f s\n',total_time(init));
    
    if init>5 && ~rem(init,5)
        save('times.mat','esequence','time_to_read','time_to_link_hard','total_time','tracking_time','segmentation_time','time_to_filter','TimeSpentMerging');
    end
    
end

esequence = master_esequence;
save('times.mat','esequence','time_to_read','time_to_link_hard','total_time','tracking_time','segmentation_time','time_to_filter','TimeSpentMerging');
%Write out a backup
%saveGreedyNucleiFiles(master_esequence,init,[outputdirectory,suffix,embryonumber,'/nuclei'],anisotropy,ROIxmin,ROIymin);
%zip(zipname,[outputdirectory,suffix,embryonumber,'/nuclei']);
%rmdir([outputdirectory,suffix,embryonumber],'s');

'lineaging completed'