%lineage driver
%itializes structures
%does easy, gathers candidates, noneasy 1-1, divisions, division classify
%{
if (~isfield(trackingparameters, 'useAveragePoint'))
    trackingparameters.useAveragePoint=false;
end
if(trackingparameters.useAveragePoint)
    for t=1:length(esequence)
        esequence{t}.finalpoints(:,3)=esequence{t}.finalaveragepoints(:,3);
    end
end
%}
%#function distanceCostFunction divScoreModelCostFunction
%#function NaiveBayes
%#function NaiveBayes\predict
%#function predict
%#function wibble


[trackingparameters,esequence]=initializeTrackingStructures_REALTIME(esequence,trackingparameters,parameters,init);

easy = tic;
%linking
esequence=linkEasyCases_REALTIME(esequence,trackingparameters,init);
esequence=gatherEndCandidates_REALTIME(esequence,trackingparameters,init);
time_to_link_easy(init) = toc(easy);

%do non gap
%these arent really parameters but internal configuration of endmatch in
%next block
hard = tic;
trackingparameters.dotrack=true;
trackingparameters.trackdiv=false;
trackingparameters.tracknondiv=true;
trackingparameters.gapthresh=0;
global answers
answers=[];
wrong=[];
numanswered=[];
result=[];
firstloop=true;
global alldatanondiv;

for offset=trackingparameters.minnondivscore:trackingparameters.nondivscorestep:trackingparameters.maxnondivscore%.875%1 %.1:.0125:.9
    
    nondivthresh=(offset);%.01*1.175^offset;
    alldatadiv=[];
    alldatanondiv=[];
    trackingparameters.endscorethresh_nondiv=nondivthresh;
    esequence=greedyEndScore_REALTIME(esequence,trackingparameters,parameters,init,1);
   %{
    if(firstloop)
        firstloop=false;
        if(trackingparameters.trainingmode)
            allbifurcationinfo(lin).initialnondivdistribution=alldatanondiv(:,3);
        end
        %trackingparameters.endnodivthreshold=prctile(alldatanondiv(:,3),91);
        %nondivthress(lin)=trackingparameters.maxnondivscore;
    end
    %}
    %compute correct, wrong if in recordanswer mode with answer key
    if(~isempty(answers))
        isdivision=answers(:,15)~=-1;%is chosen answer division
        wrong=[wrong;length(find(~answers(:,3)&isdivision))...
            ,length(find(~answers(:,4)&~isdivision))];
        numanswered=[numanswered;length(find(isdivision)),length(find(~isdivision))];                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             [length(find(isdivision)),length(find(~isdivision)),length(find(~answers(:,3)&isdivision))...
            ,length(find(~answers(:,4)&~isdivision))]
        result=[result;[length(find(isdivision)),length(find(~isdivision)),length(find(~answers(:,3)&isdivision))...
            ,length(find(~answers(:,4)&~isdivision))]];
    end
    
end

if(isfield(trackingparameters,'hysteresis')&&trackingparameters.hysteresis)
    esequence=cleanUnlinkedHysteresis(esequence,trackingparameters);  
end

%division closing
trackingparameters.do_enddiv=true;
trackingparameters.dotrack=true;
trackingparameters.trackdiv=true;
trackingparameters.tracknondiv=false;
trackingparameters.gapthresh=0;
answers=[];
%do reasonableish divisions in score order
for offset=trackingparameters.mindivscore:trackingparameters.divscorestep:trackingparameters.maxdivscore %1.5 %.1:.0125:.9
    alldatadiv=[];
    alldatanondiv=[];
    divthresh=offset;%.01*1.175^offset;
      trackingparameters.endscorethresh_div=divthresh;
    esequence=greedyEndScore_REALTIME(esequence,trackingparameters,parameters,init,1);
    %min score div, min score nondiv,
    %divright,nondivright,answerpresent,answer, answertime, i,t
    %min division is not min
    %isdivision=answers(:,7)~=-1;%is real answer division
    if(~isempty(answers))
        isdivision=answers(:,15)~=-1;%is chosen answer division
        wrong=[wrong;length(find(~answers(:,3)&isdivision))...
            ,length(find(~answers(:,4)&~isdivision))];
        numanswered=[numanswered;length(find(isdivision)),length(find(~isdivision))];
        
        [length(find(isdivision)),length(find(~isdivision)),length(find(~answers(:,3)&isdivision))...
            ,length(find(~answers(:,4)&~isdivision))]
    end
end


if(~exist('skipbifurcation')||~skipbifurcation)

    %do all remaining possible divisions
trackingparameters.endscorethresh_div=inf;
esequence=greedyEndScore_REALTIME(esequence,trackingparameters,parameters,init,0);

%time adjustment back for answer key mode
if(trackingparameters.trainingmode==true)
    trackingparameters.recordanswers=true;
    endtime=endtime-11;
    trackingparameters.endtime=trackingparameters.endtime-11;
else
    trackingparameters.recordanswers=false;%true;
end

%initialize data structures for storing bifurcation results
global bifurcationPoints;
bifurcationPoints=[];
global FNtype;
FNtype=[];
global simpleFNcorrect;
simpleFNcorrect=[];
global removed
global removedi;
global BifurcationMeasures;
global bestCandidateInfo;
global endIndicies;
global bestEndCandidateInfo;
global confidenceData;
global splitFNMatchScore;
global  computedclassificationvector;
global refclassificationvector;
global classround;
classround=[];
refclassificationvector=[];
computedclassificationvector=[];
splitFNMatchScore.forwardd1xy=[];
splitFNMatchScore.forwardd2xy=[];
splitFNMatchScore.forwardd1z=[];
splitFNMatchScore.forwardd2z=[];
splitFNMatchScore.backxy=[];
splitFNMatchScore.backz=[];
splitFNMatchScore.forwardd1gapsize=[];
splitFNMatchScore.forwardd2gapsize=[];
splitFNMatchScore.backgapsize=[];

confidenceData=[];
confidenceData.bifcon=[];
confidenceData.bifconvector=[];
confidenceData.bestbackconv=[];
confidenceData.bestforward1conv=[];
confidenceData.bestforward2conv=[];
confidenceData.rforwardd1_conv=[];
confidenceData.rforwardd1_consum=[];
confidenceData.rforwardd1_flength=[];
confidenceData.rforwardd1_solidlength=[];
confidenceData.rforwardd2_conv=[];
confidenceData.rforwardd2_consum=[];
confidenceData.rforwardd2_flength=[];
confidenceData.rforwardd2_solidlength=[];
endIndicies=[];
bestEndCandidateInfo=[];
bestCandidateInfo=[];
removedi=[];
removed=[];
BifurcationMeasures=[];



esequence=greedydeleteFPbranches_REALTIME(esequence,trackingparameters,parameters,init);


time_to_link_hard(init) = toc(hard);

end


