function [trackingparameters,esequence]=initializeTrackingStructures_REALTIME(esequence,trackingparameters,parameters,init)
%initialize data structures and perform precomputations for tracking 

if init == 2
    [trackingparameters,esequence]=initializeTrackingStructures_REALTIME(esequence,trackingparameters,parameters,1);
end

%initialize structures for sucessors and deletion
t = init;
i = init;
time = init;

if(isempty(esequence{t}))
    esequence{t}.finalpoints=[];
end
esequence{t}.suc=-1*ones(size(esequence{t}.finalpoints,1),2);
esequence{t}.suc_time=-1*ones(size(esequence{t}.finalpoints,1),2);
esequence{t}.pred=-1*ones(size(esequence{t}.finalpoints,1),1);
esequence{t}.pred_time=-1*ones(size(esequence{t}.finalpoints,1),1);
esequence{t}.delete=zeros(size(esequence{t}.finalpoints,1),1);

esequence{t}.mergedlogoddssum=zeros(size(esequence{t}.finalpoints,1),1);


if(~isempty(esequence{i}.finalpoints))
%initialize avg nn distance
distances=distance_anisotropic(esequence{i}.finalpoints',esequence{i}.finalpoints',trackingparameters.anisotropyvector);
for j=1:max(size(distances))
    distances(j,j)=Inf;
end
mindistances=min(distances);
diam=mean(mindistances);

 trackingparameters.forwardcutoff(i)=diam*trackingparameters.candidateCutoff;%1.5;%candidate cutoff (distance)
 esequence{i}.selfdistance=mindistances;


%integrate GFP
[integratedGFP,area]=integrateGFP(esequence{i},parameters);
esequence{i}.totalGFP=integratedGFP;
esequence{i}.avgGFP=integratedGFP./area;
else
    trackingparameters.forwardcutoff(i)=-1;%1.5;%candidate cutoff (distance)
    esequence{i}.selfdistance=[];
    esequence{i}.totalGFP=[];
    esequence{i}.avgGFP=[];
end


if (isfield(trackingparameters,'abscutoff')&&trackingparameters.abscutoff)
    'working'
    trackingparameters.forwardcutoff=trackingparameters.candidateCutoff*ones(size(trackingparameters.forwardcutoff));
end


%calculate confidences
if(~isempty(esequence{time}.finalpoints))

    confidencedata=calculateConfidenceVector(esequence{time},parameters);
    esequence{time}.confidencevector=confidencedata;
    tempconfidences=mvnpdf(confidencedata, trackingparameters.model.con_goodmean,trackingparameters.model.con_goodstd)...
        ./mvnpdf(confidencedata, trackingparameters.model.con_badmean,trackingparameters.model.con_badstd);
    tempconfidences=1./tempconfidences; %pdf ratio of bad
    %use raw confidences
    tempconfidences(isnan(tempconfidences))=Inf;%...
    esequence{time}.confidences=tempconfidences;
else
    esequence{time}.confidencevector=[];
        esequence{time}.confidences=[];
end



