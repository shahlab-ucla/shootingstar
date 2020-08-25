function [convector,consum,flength, solidlength]=computeConfidenceRecursive(esequence,t,i,trackingparameters)
%compute confidence vector length and ammount spent in solid parts along a forward branch
%called on best match

convector=[];
consum=0;
flength=0;
solidlength=0;
torig=t;

next=esequence{t}.suc(i,1);
nextt=esequence{t}.suc_time(i,1);
while(next~=-1)
    convector=[convector;esequence{t}.confidencevector(i,:)];
    % consum=consum+esequence{t}.confidences(i,:);
    t=nextt;
    i=next;
    nextt=esequence{t}.suc_time(i,1);
    next=esequence{t}.suc(i,1);
    solidlength=solidlength+1/trackingparameters.interval;
end
convector=[convector;esequence{t}.confidencevector(i,:)];
consum=consum+esequence{t}.confidences(i,:);
solidlength=solidlength+1/trackingparameters.interval;
%compute best match forward and recursively compute forward solidity on it

flength=(t-torig+1)/trackingparameters.interval;
consum=1;

%I added a check to minimize recursion which I think is the problem in
%longer periods and why tracking is taking way longer per frame  over
%longer time periods once you hit this length, it has little diagnostic
%value
if(solidlength>=10)
    
    return
else
    
    [reversecand,reversecandt] =findForwardCandidatesTime(esequence,i,t,trackingparameters);
    if(isempty(reversecand))
        return
    else
        [scores_forwardfromback,~,~,~]=...
            gapScore(esequence,i,t,reversecand,reversecandt,trackingparameters);
        [~,minirf]=min(scores_forwardfromback);
        
        
        [convectorf,consumf,flengthf, solidlengthf]=...
            computeConfidenceRecursive(esequence,reversecandt(minirf),reversecand(minirf),trackingparameters);
        
        
        %end
        
        
        %compute min size of conflict leading to best candidate and follow that
        %set t, i to min daughter or best score
        convector=[convector;convectorf];
        consum=consumf+consum;
        flength=flength+flengthf+(reversecandt(minirf)-t)/trackingparameters.interval;
        solidlength=solidlength+solidlengthf;
    end
end


end

