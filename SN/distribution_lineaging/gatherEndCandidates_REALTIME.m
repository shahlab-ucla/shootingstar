function esequence=gatherEndCandidates_REALTIME(esequence,trackingparameters,init)
%gathers candidates for frame to frame matching of non easy cases and
%stores them in esequence data structure

%these are included here and not in param file because they are legacy 
%parameters whose setting reflects current use of this function (gathering
%candidates for non gap div and 1-1 unlike the possible use of it to gather
%gap canididates also
trackingparameters.nnnumber_gap=2;
trackingparameters.temporalcutoff=1;%dont do gaps
trackingparameters.forwardnnnumber_gap=7; %separate threshold for t>t+1

rocdata=[];
allresults={};
counter=1;


%init cand cells for all nuclei
for t=init-1:init
    esequence{t}.backcandidates=cell(1,size(esequence{t}.finalpoints,1));
    esequence{t}.forwardcandidates=cell(1,size(esequence{t}.finalpoints,1));
end

rc=0;
for t=init-1:init
    for i=1:size(esequence{t}.finalpoints,1)
        %if this is a start
        if(esequence{t}.pred(i,1)==-1)
            [candidates,candidates_t]=findBackwardCandidatesTime_NN(esequence,i,t,trackingparameters);
 
     if(~isempty(candidates))
       rc=rc+1;  
        esequence{t}.backcandidates{i}=[esequence{t}.backcandidates{i};[candidates,candidates_t]];
        for j=1:length(candidates)
            esequence{candidates_t(j)}.forwardcandidates{candidates(j)}=[esequence{candidates_t(j)}.forwardcandidates{candidates(j)};[i,t]];
        end
     end
        end
    end
end


for t=init:init-1
    for i=1:size(esequence{t}.finalpoints,1)
        if (~isempty(esequence{t}.forwardcandidates{i}))
            %filter to have no more than set # at each timepoint
            %remove from forward candidate list, and remove itself from backward
            %list of each removed forward candidate
            candidates=esequence{t}.forwardcandidates{i}(:,1);
            candidates_t=esequence{t}.forwardcandidates{i}(:,2); 

             
            for time=min(candidates_t):min(candidates_t)+1%max(candidates_t)              
                  
                if(time==t+1)
                    cutoff=trackingparameters.forwardnnnumber;
                else
                    cutoff=trackingparameters.forwardnnnumber_gap;
                end
                if(time==t+1)
                    thistime=find(candidates_t==time);
                    thistime_logical=(candidates_t==time);
                else
                       thistime=find(candidates_t>=time);
                    thistime_logical=(candidates_t>=time); 
                end
                if(cutoff<length(thistime)) %more than allowed exist

                    distances=zeros(size(thistime));
                    for j=1:length(thistime)
                        distances(j)=distance_anisotropic(esequence{t}.finalpoints(i,:)',esequence{candidates_t(thistime(j))}.finalpoints(candidates(thistime(j)),:)',trackingparameters.anisotropyvector);
                    end
                    bad=ones(size(distances));
                    for j=1:cutoff
                        [val,imin]=min(distances);
                        bad(imin)=0;
                        distances(imin)=inf;
                    end
                    %now have the bad ones
                    discards=zeros(size(candidates));
                    discards(thistime(logical(bad)))=1;
                    %remove them from forward candidates
                    esequence{t}.forwardcandidates{i}=[candidates(~discards),candidates_t(~discards)];
                    %remove this from back candidates
                    for j=1:length(discards)
                            if(discards(j))
                        succandidates=esequence{candidates_t(j)}.backcandidates{candidates(j)};
                        equalentry=succandidates(:,1)==candidates(j)&...
                            succandidates(:,2)==candidates_t(j);
                        esequence{candidates_t(j)}.backcandidates{candidates(j)}=...
                            succandidates(~equalentry,:);
                            end
                    end
                    %update working list to reflect removal
                    candidates=esequence{t}.forwardcandidates{i}(:,1);
                    candidates_t=esequence{t}.forwardcandidates{i}(:,2); 
       
                end
            end%fore each time
        end %are candidates
    end %nucleus loop
end%time loop
