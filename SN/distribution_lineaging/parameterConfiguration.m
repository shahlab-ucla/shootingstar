%tracking parameter file used for model retraining, these should match
%final parameters used for tracking 
%revised 2/26/2019 only in that it loads a model with current classifier
%type (not necessary for correctness but I deleted all the old no longer
%working model files from the codebase)


%GENERAL IMAGING PARAMS
% trained division and bifurcation classifier model
%load 'clean_red_model.mat';
%load 'redmodel_newnondiv_4nn_recursive_wlinkconfidence_allnewmeasures_dimreduction.mat';
%load 'redmodel_newnondiv_4nn_recursive_wlinkconfidence_leaveout.mat';
%load 'redmodel_newnondiv_4nn_recursive_wlinkconfidence_leaveout.mat';
%load 'red_5_confidence_trainingdata_allfeatures_stable_leaveout_halfdivbad.mat';
%load 'red_5_confidence_trainingdata_allfeatures_stable_leaveout_halfdivbad_differentnondivbinning.mat';
load '2019TrackingModelv2.mat';
%load 'clean_red_singlemodel_red_kernel.mat';
%load 'clean_red_singlemodel_red_normal.mat';
%load 'clean_red_multimodel_red_normal.mat';
%load 'hobert_greenmodel_multimodel_normal.mat'
%load 'hobert_greenmodel_v2.mat'
%load '10min_yalezebrafish_singlemodel.mat'

%start time for analysis
trackingparameters.starttime=1;

% Whether to force all tentative bifurcations and perform bifurcation
% resolution
%if true it creates bifurcations up to trackingparameters.maxdivscore
%(defined below) and returns result as final lineage. If false then creates
% all possible bifurcations and runs tentative bifurcation resolution 
skipbifurcation=false;

% time interval (in minutes) between imageframes
trackingparameters.interval=1;

% whether to use the average of segmented z planes as cell centroid or the
% location of the 3D intensity maxima
% maxima appears to be more reliable in our data
trackingparameters.useAveragePoint=false;

%EASY STAGE
%max multiple of avg nn distance for timestep used as max sanity check on match distance
trackingparameters.candidateCutoff=1.2;

%UNLIKELY TO CHANGE
%safe filter easy cases
trackingparameters.safefilter=true;
trackingparameters.safefactor=2; %2 = normal safe, can be more conservative
%filter out cases with a conflicting claim
trackingparameters.conflictfilter=true;
%candidate selection
% # of back nn compiled into candidate list used in div and noneasy 1-1
trackingparameters.nnnumber=2;
%filter 2 nn back candidates based on nn forward ranking
%prevent more than n forward nn to be on candidate list
trackingparameters.forwardnnnumber=4; 

%TENTATIVE BIFURCATION CREATION

%current non div cost function (normalized distance)
trackingparameters.nonDivCostFunction=@distanceCostFunction;%
%trackingparameters.nonDivCostFunction=@nondivScoreModelCostFunction;%
%start and end distances (in avg nn units) of optimization of distance for
%noneasy 1-1 cases
trackingparameters.minnondivscore=.125;% 
trackingparameters.nondivscorestep=.125;
trackingparameters.maxnondivscore=.875;

%likely ranges for if using 1:1 liklihood score instead of distance
%trackingparameters.minnondivscore=-20;
%trackingparameters.nondivscorestep=1;
%trackingparameters.maxnondivscore=45;

%cost function used for division scoring/ feeding into bifurcaton model
trackingparameters.DivCostFunction=@divScoreModelCostFunction;
%threshold and increment for sliding threshold optimization of divisions
%log product of pdf score for division models
trackingparameters.mindivscore=-20;
trackingparameters.divscorestep=2;
trackingparameters.maxdivscore=5;


%parameters for deleting fragments so isolated they cannot form bifurcations
trackingparameters.deleteisolated=true; %whether to do any deletion of isolated points
%delete if
%smaller than thresh
trackingparameters.FPsizethresh=2;
%and less than # of cells
trackingparameters.earlythresh=250;
%or if smaller than or = to second thresh
trackingparameters.FPsizethreshsmall=1;
%note that these are currently redundant with eachother
%and with third non parameterized check that discards anything isolated of length 1
%but might be useful to be able to set them independently higher in some
%situations

%TENATATIVE BIFURCATION RESOLUTION
%above this size things are considered unlikely to be FP so their FN forward gap properties are not %taken into account when judging them
trackingparameters.smallcutoff=4; 
%max temporal gap that can be closed during bifurcation processing stage
trackingparameters.temporalcutoff=6;

%starting negative cutoff for FN search (2=FN only, 1=FN or calssifier 1:1)
%though it can be changed anytime only makes sense to use it in the same
%mode the model has been trained in
trackingparameters.temporalcutoffstart=2;
%note trackingparameters.candidateCutoff is used in this processing step also

%final confidence calculation
%temporal size of wider window in estimating properties for temporal neighborhood
%around link
trackingparameters.wideWindow=4;
