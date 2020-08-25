%{
%remember you cano nly run this for the embryo that is currently in
% workspace and need to reload it for others to get points, dumbass

emb=lin;
cases=allbifurcationinfo(emb).refclassificationvector==0;
cases=allbifurcationinfo(emb).refclassificationvector~=allbifurcationinfo(emb).computedclassificationvector;
cases=allbifurcationinfo(emb).refclassificationvector==0&allbifurcationinfo(emb).classround'==2;

info=[];
for i=1:length(cases)
    if(cases(i))
        info=[info;allbifurcationinfo(emb).removed(i,1),...
            esequence{allbifurcationinfo(emb).removed(i,1)}.finalpoints(...
        allbifurcationinfo(emb).removed(i,2),:)];
    end
end
if(ROI)
info(:,2)=info(:,2)+ROIxmin;
info(:,3)=info(:,3)+ROIymin
end
  info=[info,allbifurcationinfo(emb).refclassificationvector(cases)',allbifurcationinfo(emb).computedclassificationvector(cases)'];

cmats={};
for i=1:length(allbifurcationinfo)
cmats{i}=confusionmat(allbifurcationinfo(i).refclassificationvector,allbifurcationinfo(i).computedclassificationvector);
end
%}
%evaluate error over set of lineages for marker 1
basedir='L:\duz\project\Imaging\';

%{
lineages={...
    'ZD_RW10348_WT_20110126_2_s1_emb2','ZD_RW10348_WT_20110126_2_s1_emb3',...
    'ZD_RW10348_WT_20110126_2_s2_emb1','ZD_RW10348_WT_20110126_2_s2_emb2',...
    'ZD_RW10348_WT_20110427_3_s1_emb1','ZD_RW10348_WT_20110427_3_s1_emb2',...
    'ZD_RW10348_WT_20110427_3_s2_emb1','ZD_RW10348_WT_20110427_3_s2_emb2',...
    'ZD_RW10348_WT_20110427_3_s4_emb1','ZD_RW10348_WT_20110427_3_s4_emb2'...
'ZD_RW10434_20110126_2_s3_emb1',...
'ZD_RW10434_20110126_2_s3_emb2',...
'ZD_RW10434_20110126_2_s4_emb1',...
'ZD_RW10434_20110126_2_s4_emb2',...
'ZD_RW10434_20110126_2_s4_emb3',...
'ZD_RW10434_WT_20110429_2_s1_emb1',...
'ZD_RW10434_WT_20110429_2_s3_emb1',...
'ZD_RW10434_WT_20110502_1_s1_emb1',...
'ZD_RW10434_WT_20110502_1_s3_emb1',...
'ZD_RW10434_WT_20110502_1_s3_emb2',...
'ZD_BV82_WT_20110419_3_s1_emb1',...
'ZD_BV82_WT_20110419_3_s1_emb2',...
'ZD_BV82_WT_20110419_3_s2_emb1',...
'ZD_BV82_WT_20110419_3_s2_emb2',...
'ZD_BV82_WT_20110419_3_s2_emb3',...
'ZD_BV82_WT_20110426_1_s1_emb2',...
'ZD_BV82_WT_20110426_1_s1_emb3',...
'ZD_BV82_WT_20110426_1_s1_emb4',...
'ZD_BV82_WT_20110426_1_s2_emb2',...
'ZD_BV82_WT_20110725_2_s1_emb1',...
'ZD_RW10425_WT_20100412_2_s1_emb1',...
'ZD_RW10425_WT_20100412_2_s1_emb2',...
'ZD_RW10425_WT_20100412_2_s1_emb3',...
'ZD_RW10425_WT_20100412_2_s2_emb1',...
'ZD_RW10425_WT_20100412_2_s2_emb2',...
'ZD_RW10425_WT_20110428_3_s1_emb1',...
'ZD_RW10425_WT_20110428_3_s1_emb3',...
'ZD_RW10425_WT_20110428_3_s2_emb1',...
'ZD_RW10425_WT_20110428_3_s2_emb3',...
'ZD_RW10425_WT_20110428_3_s3_emb1',...
'ZD_RW10714_WT_20110724_1_s1_emb1',...
'ZD_RW10714_WT_20110724_1_s1_emb2',...
'ZD_RW10714_WT_20110724_1_s1_emb3',...
'ZD_RW10714_WT_20110724_1_s2_emb1',...
'ZD_RW10714_WT_20110724_1_s2_emb2',...
'ZD_RW10714_WT_20110724_1_s2_emb3',...
'ZD_RW10714_WT_20110724_1_s3_emb1',...
'ZD_RW10714_WT_20110724_1_s3_emb2',...
'ZD_RW10714_WT_20110724_1_s3_emb3',...
'ZD_RW10714_WT_20110726_1_s1_emb1',...
};
%}
lineages={'ZD_BV82_SKN-1_20110422_1_s2_emb1','ZD_BV82_SKN-1_20110422_1_s2_emb2','ZD_BV82_SKN-1_20120802_1_s1_emb3',...
'ZD_BV82_SKN-1_20120801_2_s1_emb3','ZD_BV82_MEX-3_20110507_2_s3_emb1','ZD_BV82_MEX-3_20110507_2_s3_emb2',...
'ZD_BV82_MEX-3_20110527_1_s1_emb1','ZD_BV82_MEX-3_20110328_2_s2_emb1','ZD_BV82_MEX-3_20110328_2_s2_emb2',...
'ZD_BV82_MEX-3_20110328_2_s2_emb3','ZD_BV82_GLP-1_20110318_2_s2_emb2','ZD_BV82_GLP-1_20110318_2_s2_emb3',...
'ZD_BV82_PAL-1_20110318_2_s2_emb1','ZD_BV82_PAL-1_20110318_2_s2_emb2'};

lineageimage=lineages;
for i=1:length(lineages)
    lineageimage{i}=lineageimage{i}(1:end-5)
end
lineagedir=lineages;
for i=1:length(lineages)
    lineagedir{i}=lineagedir{i}(1:end-8)
end


edittimes=400*ones(size(lineagedir));

eightcelltimes=30*ones(size(lineagedir));
%eightcelltimes=[35,35,25,25,35,25,35,35,35,35];
eightcelltimes=eightcelltimes+20;%go to AB8 instead of P08
nondivthress=[];
%lineages={'ZD_RW10348_HDA-1_20110331_1_s3_emb3'};
%edittimes=[300];
%red 10 training set, when these are used filtering needs to be on in
%answer key generation

%11-14 use this
%edittimes=[190,200,210,200,190,   180,190,180,190,180];
edittimes=200.*ones(1,14);
%{
    lineages={'ZD_BV82_APX-1_20110415_1_s2_emb1','ZD_BV82_CUL-1_20110329_1_s1_emb1',...
    'ZD_RW10425_ELT-1V_20110916_1_s1_emb2','ZD_BV82_LIT-1_20110419_1_s1_emb2',...
    'ZD_RW10434_LIT-1_20110419_1_s4_emb1',...
    'ZD_BV82_WT_20110419_3_s2_emb2','ZD_BV82_WT_20110426_1_s1_emb3',...
    'ZD_BV82_WT_20110426_1_s1_emb4', 'ZD_BV82_WT_20110426_1_s2_emb2',...
    'ZD_RW10425_WT_20100412_2_s1_emb1'};
%}
% basedir='L:\santella\unzipped_lineages\training\'


%}
%green
%edittimes=280;
%lineages={'ZD_BV24_JournalV_1_s1_emb1'};



%red 5 training better for doing confidence training because lack 'wrong'
%polar links
%when these are used for training filtering should be off
%{
lineages={'ZD_BV82_WT_20100809_2_s1_emb_linbiftest_lowthresh1','ZD_RW10425_WT_20100412_2_s1_emb1',...
    'ZD_RW10425_WT_20100412_2_s1_emb2','ZD_RW10425_WT_20100412_2_s1_emb3','ZD_RW10434_WT_20110429_2_s1_emb_linbiftest1'};

lineagedir={'ZD_BV82_WT_20100809_2','ZD_RW10425_WT_20100412_2',...
    'ZD_RW10425_WT_20100412_2','ZD_RW10425_WT_20100412_2','ZD_RW10434_WT_20110429_2'};
lineageimage={'ZD_BV82_WT_20100809_2_s1','ZD_RW10425_WT_20100412_2_s1',...
    'ZD_RW10425_WT_20100412_2_s1','ZD_RW10425_WT_20100412_2_s1','ZD_RW10434_WT_20110429_2_s1'};
edittimes=[200,192,201,185,185]
eightcelltimes=[20,20,30,20,20];
basedir='L:\santella\unzipped_lineages\test_data\'
%}




%lineages={'ZD_RW10714_WT_20101012_3_s3_emb1','ZD_RW10425_WT_20110428_3_s1_emb1','ZD_BV82_PIE-1_20110315_1_s2_emb1',...
%    'ZD_RW10348_PAL-1_20110318_2_s3_emb1','ZD_RW10348_HDA-1_20110331_1_s3_emb3'};
%edittimes=[170,160,180,165,200];

%edittimes=edittimes-10;%avoid area where cant recover fn at end


%'green mode'
%parameterConfigurationGreen
'red mode'
parameterConfiguration

%trainingmode refers to bif. classifier training not confidence
%uses oracle
trackingparameters.trainingmode=false;
trackingparameters.recordanswers=false;
answerkey=false;
outputtrimmed=true;
errors=cell(1,length(lineages));
replacemodel=false;
allexpectedchange={};
evalforced=false;
evalfinal=false;
allbifurcationinfo=[];

%for lin=1:length(lineages)
    for lin=12:length(lineages)
    endtime=edittimes(lin);
    trackingparameters.endtime=endtime;
    %'loader for files in zhuos directory'
     load([basedir,lineagedir{lin},'/',lineages{lin},'_fullmatlabresult.mat']);
    'loader for training files'
    % load([basedir,lineages{lin},'_fullmatlabresult.mat']);
    trackingparameters.anisotropyvector=[1,1,anisotropy];
    parameters.anisotropyvector=[1,1,anisotropy];
    
    embryonumbers = {};
    nucleidir=basedir;
    embryonumbers_c={[lineages{lin},'_edited\nuclei\']};
    %endtime=edittimes(lin)+11;
    outputdirectory=[lineages{lin},'/nuclei/'];
    if (answerkey)
        train_tracking_statistics_function;
        train_confidence_function;
    end
    endtime=endtime+10;
    trackingparameters.endtime=trackingparameters.endtime+10;
    
    %for live run reset to all timepoints
%   'process all time points for live trim'
%    
%   'process 250 for top down'  
% trackingparameters.endtime=min(250,length(esequence)); 
% trackingparameters.endtime=length(esequence);
 %       endtime=trackingparameters.endtime;
    
  
    
    trackingparameters.recordanswers=false;
    
    tic
    tracking_driver_new_classifier_based_version;
    toc
    
    %time adjustment back for final evaluation
    %endtime=endtime-9;
    %trackingparameters.endtime=trackingparameters.endtime-9;
    %reset back for final evaluation
    if(trackingparameters.trainingmode==false)
        endtime=endtime-11;
        trackingparameters.endtime=trackingparameters.endtime-11;
    end
    
    
    %section for storing training info
    ncells=[];
    for i=1:size(removed,1)
        % ncells=[ncells,length(esequence{removed(i,1)}.FP)];
        ncells=[ncells,mean(esequence{removed(i,1)}.selfdistance)./mean(esequence{removed(i,1)}.finaldiams)];
    end
    if (trackingparameters.trainingmode) %training mode reall means training for bifurcation model
        allbifurcationinfo(lin).Divdata=Divdata;
        allbifurcationinfo(lin).Tripledata=Tripledata;
        allbifurcationinfo(lin).NoDivdata=NoDivdata;
        allbifurcationinfo(lin).ncells=ncells;
        allbifurcationinfo(lin).removed=removed;
        allbifurcationinfo(lin).confidenceData=confidenceData;
        allbifurcationinfo(lin).splitFNMatchScore=splitFNMatchScore;
        allbifurcationinfo(lin).BifurcationMeasures=BifurcationMeasures;
        allbifurcationinfo(lin).classround=classround;
        allbifurcationinfo(lin).computedclassificationvector=computedclassificationvector;
        allbifurcationinfo(lin).refclassificationvector= refclassificationvector;
    end
    if(~outputtrimmed)
        [ linkconfidencedata] ...
            = extractTrainingConfidenceDataVectors( esequence,trackingparameters,embryonumbers_c,nucleidir,ROI,ROIxmin,ROIymin  );
        allbifurcationinfo(lin).linkconfidencedata=linkconfidencedata;
    end
    if(evalfinal)
        
        mkdir(outputdirectory);
        saveGreedyNucleiFiles(esequence,endtime,outputdirectory,anisotropy,ROIxmin,ROIymin);
        zipname=[lineages{lin},'/',embryonumber,'_',suffix,'.zip'];
        zip(zipname,[outputdirectory,'']);
        
        uneddir=outputdirectory;
        eddir=[basedir,lineages{lin},'_edited\nuclei\'];
        %    uneddir=[basedir,lineages{i},'\nuclei\'];
        test=evaluate_lineage_error(eddir,uneddir,endtime,anisotropy);
        test{37}=test{10}+test{21}+test{32};
        test{38}=test{12}+test{23}++test{34};
        test{39}=sum(test{37}+test{38});
        test{40}=sum(test{37});
        test{41}=sum(test{38});
        
        errors{lin}=test;
    end
    if(trackingparameters.trainingmode)
        allexpectedchange{lin}=expected_corrections;
    end
    
    
    if(outputtrimmed)
        esequence_con= scoreLinkConfidence(esequence,trackingparameters);
        %{ 
        %cutting trimming version
        for t=1:trackingparameters.endtime-1
            for i=1:length(esequence_con{t}.suc)
                if(esequence_con{t}.linkconfidences(i)<.8)
                    suct=esequence_con{t}.suc_time(i,:);
                    suc=esequence_con{t}.suc(i,:);
                     if(suc(1)~=-1)
                    %wipe out pred of suc
                    esequence_con{suct(1)}.pred(suc(1))=-1;
                    esequence_con{suct(1)}.pred_time(suc(1))=-1;
                     end
                    if(suc(2)~=-1)
                        esequence_con{suct(2)}.pred(suc(2))=-1;
                        esequence_con{suct(2)}.pred_time(suc(2))=-1;
                    end
                    esequence_con{t}.suc(i,:)=[-1,-1];
                    esequence_con{t}.suc_time(i,:)=[-1,-1];
                end
            end
        end
        %}
        %clipped trimming version
      %{
        %initialize
        for t=1:trackingparameters.endtime
            esequence_con{t}.path_confidence=zeros(size(esequence{t}.delete));
        end
        %assign path confidence
        startt=eightcelltimes(lin);
        for i=1:size(esequence_con{startt}.finalpoints,1)
            esequence_con=recursiveComputePathConfidence_min(esequence_con,startt,i,1);
        end
        %delete below .8 previously
        for t=startt:trackingparameters.endtime
            esequence_con{t}.delete(esequence_con{t}.path_confidence<.8)=true;
        end
      %}
        %output
         % finaloutputdirectory='l:/santella/lineage_automerge/80version_liberalmodel_ab8_newversiontrim/';
        finaloutputdirectory='L:\santella\lineage_automerge\80version_conservativemodel_ab8_wholescore\';
          finaloutputdirectory='L:\santella\lineage_automerge\mutantanalysis\';
    
        trimmed_outputdirectory=['nuclei/'];
        mkdir(trimmed_outputdirectory);
        %note I'm not surw where code using the anlge and aleft is, or
        %whats up with it
        saveGreedyNucleiFilesAndConfidence(esequence_con,trackingparameters.endtime, trimmed_outputdirectory,anisotropy,ROIxmin,ROIymin,0,false);
        zipname=['trim_',embryonumber,'_',suffix,'.zip'];
        
        zip(zipname, trimmed_outputdirectory);
        rmdir( trimmed_outputdirectory,'s');
        
      
                          
        movefile (zipname, finaloutputdirectory)
        outputXMLfile([finaloutputdirectory,lineages{lin},'_edited.xml'],...
            xyres,zres,slices,[finaloutputdirectory,zipname],...
            ['l:/duz/project/imaging/',lineagedir{lin},'/'],...
            lineageimage{lin},false);
        %}
          %save result as mat file
     %  save([finaloutputdirectory,lineages{lin},'.mat']);
    end
     %finaloutputdirectory='L:\santella\lineage_automerge\80version_conservativemodel_ab8_wholescore\';
        save([finaloutputdirectory,lineages{lin},'.mat']);
   
    
end
%{
if(trackingparameters.trainingmode)
    multiple_embryo_train;
else
    multiple_embryo_train_confidence;
end

%}
    
           