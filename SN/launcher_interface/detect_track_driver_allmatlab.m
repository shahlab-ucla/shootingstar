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

function detect_track_driver_allmatlab(parameterfile,embryodir,embryonumber,suffix,outputdirectory,polygon_points,isred);%,lineageParameterFileName)
%Import AceTree classes
import org.rhwlab.acetree.AceTreeLauncher
import org.rhwlab.acetree.AceTree
import org.rhwlab.snight.NucleiMgr
%Fire up AceTree with a config file pointing to a blank directory
AceTree.main('C:/Users/shahp2/Documents/MATLAB/SN/launcher_interface/BlankConfig.xml')
%Get handles to the AceTree and NucleiManager instances that AT has created
%Java objects pulled from NucleiManager are not copies but the instances
%themselves and edits will be reflected to NucleiManager immediately.
%TODO:
%1. Write a conversion function to put esequence in an AT vector format
%Since this is currently done recursively, one might have to just build up
%the whole representation using the recursive functions directly
%2. Implement a function (mergeNucleiRecord) to diff and merge a passed Nuclei instance and the
%AT formatted representation of an Esequence entry
%3. At the beginning of each process sequence, diff every entry in the AT
%NM with the SN esequence, update preds and sucs, incorporate new nuclei
%4. To incorporate new nuclei, a function will have to be written that
%checks for the existence of preds and sucs, if one of the two exists, copy
%tracking params from that one, if both exist, average the tracking params,
%if neither exists, compute the mean tracking params of all cells at that
%timepoint
%5. At the end of each process sequence, merge new results into the AT NM
%from esequence
%6. Implement a semaphore flag in the NM class that AT will check before
%modifying and that SN will check before modifying
%7. Change the update nuclei function to just trigger a rebuild which will
%pull new entries from NM that SN has loaded into the GUI tree
AT = AceTree.getAceTree([])
NM = AT.getNucleiMgr

runexpression=true;
newimage=false;

'beginning lineaging'
tic
global parameters;


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
    fprintf (file,'<useStack type="1"/>\n');
    fprintf (file,['<image file="',embryodir,embryonumber,'_t',num2str(start_time),'.tif"/>\n']);
else
    fprintf (file,['<image file="',embryodir,'image/tif/',embryonumber,'-t',num2str(start_time,'%03d'),'-p01.tif"/>\n']);
end
fprintf (file,['<nuclei file="',zipname,'"/>\n']);

fprintf (file,'<end index="475"/>\n');
fprintf(file,['<resolution xyRes="',num2str(xyres),'" zRes="',num2str(zres),'" planeEnd="',num2str(slices),'"/> <exprCorr type="blot"/>\n']); 
fprintf (file,'<polar size="15"/>\n');
fprintf (file,'</embryo>');
fclose(file);

for init = 1:run_time %Force pipeline to run one image at a time
    if init == 1
        firsttimestep = true;
    end
    esequence={};
    tlist = init;
    start_time = init;
    end_time = init;
    processSequence; %start segmentation for the timepoint as defined by init
    master_esequence{init} = esequence{init}; %add the results of segmentation into the appropriate position of the master structure
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
        tracking_driver_new_classifier_based_version;
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
    
    if init > 5 && ~rem(init,5)
        skipbifurcation = 0;
        'starting windowed tracking'
        trackingparameters.trainingmode=false;
        trackingparameters.recordanswers=false;
        evalforced=false;
        trackingparameters.starttime=1;
        trackingparameters.endtime=10; %since we're only going to pass two timepoints, the endtime is 2
        trackingparameters.anisotropyvector=[1,1,anisotropy];
        parameters.anisotropyvector=[1,1,anisotropy];
        esequence = cell(1,10);
        for i = 1:10
            esequence{i} = master_esequence{init - (10-i)};
        end
        tracking_driver_new_classifier_based_version
        %To adjust time here, we just slide the scale based on where we are.
        %We don't have to worry about adjustments that were made
        %sequentially since those just get overwritten with the current
        %1:10 timescale
        for j = 1:10
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
        for j = 1:10
            m = length(esequence{j}.pred_time);
            for i = 1:m
                if esequence{j}.pred_time(i) ~= -1
                    esequence{j}.pred_time(i) = esequence{j}.pred_time(i) + (init - 10);
                end
            end
        end
        %windowed = batch('tracking_driver_new_classifier_based_version','AdditionalPaths',genpath('C:\Users\Pavak\Documents\MATLAB\Starrynite')); %this job can be allocated workers using matlabpool
%       jstart = init;

        for i = 1:10
            master_esequence{init - (10-i)} = esequence{i};
        end
        
        %writes out all of the nuclei files and zips them up, this
        %overwrites the zip but should be OK since AT won't save out the
        %edited tree until the very end
        saveGreedyNucleiFiles(master_esequence,init,[outputdirectory,suffix,embryonumber,'/nuclei'],anisotropy,ROIxmin,ROIymin);
        check = 1;
        while check ==1
            try
                zip(zipname,[outputdirectory,suffix,embryonumber,'/nuclei']);
                check = 0;
            catch exception
                pause(0.1);
            end
        end

    end

%    if exist('windowed','var')
%        if strcmp(windowed.State,'finished')
%            S = windowed.load('esequence');
            
%            windowed.delete()
%        end
%    end
end

esequence = master_esequence;

%if (end_time-start_time>0)
    %parameterConfiguration;
    %parameterfile=lineageParameterFileName;
    %readParameters;
%{
    trackingparameters.trainingmode=false;
    trackingparameters.recordanswers=false;
    evalforced=false;
    endtime=end_time;
    trackingparameters.endtime=endtime;
    trackingparameters.anisotropyvector=[1,1,anisotropy];
    parameters.anisotropyvector=[1,1,anisotropy];
    
    tracking_driver_new_classifier_based_version;
%}
    %{
if(start_time~=1)
        tempesequence=esequence;
        esequence=cell(end_time,1);
        for i=1:length(tempesequence)
            esequence{start_time+i-1}=tempesequence{i};
        end
    end
%}  
    
    saveGreedyNucleiFiles(esequence,run_time,[outputdirectory,suffix,embryonumber,'/nuclei'],anisotropy,ROIxmin,ROIymin);
%else
      %now done in detection 
      %{
    if(start_time~=1)
        tempesequence=esequence;
        esequence=cell(end_time,1);
        for i=1:length(tempesequence)
            esequence{start_time+i-1}=tempesequence{i};
        end
    end
      %}
    %base=[outputdirectory,suffix,embryonumber,'/nuclei/'];
    %output_unlineaged_acetree;
%end

zip(zipname,[outputdirectory,suffix,embryonumber,'/nuclei']);
%remove temp directory
rmdir([outputdirectory,suffix,embryonumber],'s');

%then need to output xml file
%and run naming program

%currentdir=pwd;
%cd ('l:/bin/starryniteII');
%system([' java -Xmx500m -cp acebatch2.jar Measure1 ',xmlname]);
%cd (currentdir);

'lineaging completed'

'running expression'

%{
if(runexpression)
if(~isred)
     system(['java -cp acebatch2.jar SixteenBitGreenExtractor1 ',xmlname,' 400']);  
  
else
     system(['java -cp acebatch2.jar SixteenBitRedExtractor1 ',xmlname,' 400']);    
end
end

%lineage processing
time=end_time;
startt=start_time;
movieInfo=[];

global parameters;

for i=startt:time
    if ~isempty(esequence{i-startt+1}.finalpoints)
        [integratedGFP,area]=integrateGFP(esequence{i-startt+1},parameters);
        
        %amp=esequence{i}.finalmaximas;
        amp=integratedGFP;
        %   pos=esequence{i}.finalaveragepoints(1:length(amp),:);
        pos=esequence{i-startt+1}.finalpoints(1:length(amp),:);
        pos(:,3)=pos(:,3)*anisotropy;
        
        %filter nuclei based on polygon bounding box
        filtervalues=inpolygon(pos(:,1)+ROIxmin,pos(:,2)+ROIymin,polygon_points(:,1),polygon_points(:,2));
        pos=pos(filtervalues,:);
        amp=amp(filtervalues);
        
        movieInfo(i).xCoord=[pos(:,1)+ROIxmin,zeros(size(pos(:,1)))];
        movieInfo(i).yCoord=[pos(:,2)+ROIymin,zeros(size(pos(:,2)))];
        movieInfo(i).zCoord=[pos(:,3),zeros(size(pos(:,3)))];
        movieInfo(i).amp=[amp,zeros(size(amp))];
        
        %also filter esequence by roi for single point output
       %
       esequence{i-startt+1}.finalpoints=esequence{i-startt+1}.finalpoints(filtervalues,:);
       esequence{i-startt+1}.finaldiams=esequence{i-startt+1}.finaldiams(filtervalues);
        esequence{i-startt+1}.finalmaximas=esequence{i-startt+1}.finalmaximas(filtervalues);
           esequence{i-startt+1}.finalpoints(:,1)=esequence{i-startt+1}.finalpoints(:,1)+ROIxmin;
              esequence{i-startt+1}.finalpoints(:,2)=esequence{i-startt+1}.finalpoints(:,2)+ROIymin;
    end
end

%remove big data structure from memory before tracking...
%clear esequence;
clear backup_esequence;


%parameter configuration




%% Cost functions

%compiler note to include these files that are not explicitly called
%#function costMatLinearMotionLink costMatLinearMotionCloseGaps kalmanResMemLM kalmanInitLinearMotion kalmanGainLinearMotion

%Frame-to-frame linking
costMatrices(1).funcName = 'costMatLinearMotionLink';

%Gap closing, merging and splitting
costMatrices(2).funcName = 'costMatLinearMotionCloseGaps';

%--------------------------------------------------------------------------

%% Kalman filter functions

%Memory reservation
kalmanFunctions.reserveMem = 'kalmanResMemLM';

%Filter initialization
kalmanFunctions.initialize = 'kalmanInitLinearMotion';

%Gain calculation based on linking history
kalmanFunctions.calcGain = 'kalmanGainLinearMotion';

%--------------------------------------------------------------------------

%% General tracking parameters

%Gap closing time window
gapCloseParam.timeWindow = 5;%10;5 me

%Flag for merging and splitting
gapCloseParam.mergeSplit = 1;

%Minimum track segment length used in the gap closing, merging and
%splitting step
gapCloseParam.minTrackLen = 2;

%--------------------------------------------------------------------------

%% Cost function specific parameters: Frame-to-frame linking

%Flag for linear motion
parameters.linearMotion = 0;

%Search radius lower limit
parameters.minSearchRadius = 2;%2

%Search radius upper limit
parameters.maxSearchRadius = 50;%100;%5

%Standard deviation multiplication factor
parameters.brownStdMult =frame2frameStd;% 2.4;%5;%8

%Flag for using local density in search radius estimation
parameters.useLocalDensity = 1;

%Number of past frames used in nearest neighbor calculation
parameters.nnWindow = gapCloseParam.timeWindow;

%Store parameters for function call
costMatrices(1).parameters = parameters;
clear parameters

%--------------------------------------------------------------------------

%% Cost cunction specific parameters: Gap closing, merging and splitting

%Same parameters as for the frame-to-frame linking cost function
parameters.linearMotion = costMatrices(1).parameters.linearMotion;
parameters.useLocalDensity = costMatrices(1).parameters.useLocalDensity;
parameters.maxSearchRadius = costMatrices(1).parameters.maxSearchRadius;
parameters.minSearchRadius = costMatrices(1).parameters.minSearchRadius;

%*** I changed thiis brownian std /2 before 1 is default .75 for test 1 .5
%conservative %2
parameters.brownStdMult = gapStd*costMatrices(1).parameters.brownStdMult*ones(gapCloseParam.timeWindow,1);
parameters.nnWindow = costMatrices(1).parameters.nnWindow;

%Gap length (frames) at which f(gap) (in search radius definition) reaches its
%plateau
parameters.timeReachConfB = 2;

%Amplitude ratio lower and upper limits
parameters.ampRatioLimit = [.25, 4];%[0.1 8];% [0.5 4];

%Minimum length (frames) for track segment analysis
parameters.lenForClassify = 4;%5

%Standard deviation multiplication factor along preferred direction of
%motion
parameters.linStdMult = gapStd*ones(gapCloseParam.timeWindow,1);
%parameters.linStdMult = 3*ones(gapCloseParam.timeWindow,1);
%Gap length (frames) at which f'(gap) (in definition of search radius
%parallel to preferred direction of motion) reaches its plateau
parameters.timeReachConfL = gapCloseParam.timeWindow;

%Maximum angle between the directions of motion of two linear track
%segments that are allowed to get linked
parameters.maxAngleVV = 90;%45;

%Store parameters for function call
costMatrices(2).parameters = parameters;
clear parameters







scriptTrackGeneral_noparam;

mkdir([outputdirectory,suffix,embryonumber,'/nuclei']);

%mkdir([outputdirectory,'nuclei']);

%hack for fact that save was written assuming esequence starts at t 1 for
%case where start time was not 1
if(start_time~=1)
tempesequence=esequence;
esequence=cell(end_time,1);
for i=1:length(tempesequence)
    esequence{start_time+i-1}=tempesequence{i};
end
end

if (end_time-start_time>0)
saveNucleiFiles(tracksFinal,esequence,time,[outputdirectory,suffix,embryonumber,'/nuclei'],anisotropy)
else
 base=[outputdirectory,suffix,embryonumber,'/nuclei/'];
 output_unlineaged_acetree;
end
zipname=[outputdirectory,embryonumber,'_',suffix,'.zip'];
zip(zipname,[outputdirectory,suffix,embryonumber,'/nuclei']);
%remove temp directory
rmdir([outputdirectory,suffix,embryonumber],'s');



if(start_time~=1)
tempesequence=esequence;
esequence=cell(end_time,1);
for i=1:length(tempesequence)
    esequence{start_time+i-1}=tempesequence{i};
end
end

if (end_time-start_time>0)
saveNucleiFiles(tracksFinal,esequence,time,[outputdirectory,suffix,embryonumber,'/nuclei'],anisotropy)
else
 base=[outputdirectory,suffix,embryonumber,'/nuclei/'];
 output_unlineaged_acetree;
end
zipname=[outputdirectory,embryonumber,'_',suffix,'.zip'];
zip(zipname,[outputdirectory,suffix,embryonumber,'/nuclei']);
%remove temp directory
rmdir([outputdirectory,suffix,embryonumber],'s');


%then need to output xml file
%and run naming program
xmlname=[outputdirectory,embryonumber,'_',suffix,'.xml'];
file=fopen(xmlname,'w');

fprintf (file, '<?xml version=''1.0'' encoding=\''utf-8\''?>\n');
fprintf (file, '<embryo>\n');

if (newimage)
    fprintf (file,'<useStack type="1"/>\n');
    fprintf (file,['<image file="',embryodir,embryonumber,'_t1.tif"/>\n']);
else
    fprintf (file,['<image file="',embryodir,'image/tif/',embryonumber,'-t001-p01.tif"/>\n']);
end
fprintf (file,['<nuclei file="',zipname,'"/>\n']);

fprintf (file,'<end index="475"/>\n');
fprintf(file,['<resolution xyRes="',num2str(xyres),'" zRes="',num2str(zres),'" planeEnd="',num2str(slices),'"/> <exprCorr type="blot"/>\n']); 
fprintf (file,'<polar size="15"/>\n');
fprintf (file,'</embryo>');
fclose(file);

%edited duplicate 
zipnameedited=[outputdirectory,embryonumber,'_',suffix,'_edited.zip'];
xmlname=[outputdirectory,embryonumber,'_',suffix,'_edited.xml'];
file=fopen(xmlname,'w');

fprintf (file, '<?xml version=''1.0'' encoding=\''utf-8\''?>\n');
fprintf (file, '<embryo>\n');
if (newimage)
    %if looking at new images delete old images after use
   % rmdir([embryodir,'image/'],'s');
    fprintf (file,'<useStack type="1"/>\n');
    fprintf (file,['<image file="',embryodir,embryonumber,'_t1.tif"/>\n']);
else
    fprintf (file,['<image file="',embryodir,'image/tif/',embryonumber,'-t001-p01.tif"/>\n']);
end
fprintf (file,['<nuclei file="',zipnameedited,'"/>\n']);

fprintf (file,'<end index="475"/>\n');
fprintf(file,['<resolution xyRes="',num2str(xyres),'" zRes="',num2str(zres),'" planeEnd="',num2str(slices),'"/> <exprCorr type="blot"/>\n']); 
fprintf (file,'<polar size="15"/>\n');
fprintf (file,'</embryo>');
fclose(file);
%edited duplicate 
copyfile (zipname,zipnameedited);

%currentdir=pwd;
%cd ('l:/bin/starryniteII');
system([' java -Xmx500m -cp acebatch2.jar Measure1 ',xmlname]);
%cd (currentdir);

'lineaging completed'

'running expression'
if(isred)
     system(['java -cp acebatch2.jar SixteenBitGreenExtractor1 ',xmlname,' 400']);  
  
else
     system(['java -cp acebatch2.jar SixteenBitRedExtractor1 ',xmlname,' 400']);    
end
%}