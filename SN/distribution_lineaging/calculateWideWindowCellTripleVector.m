function [ data ] = calculateWideWindowCellTripleVector...
    (esequence,t,i,tj,j,tk,k,anisotropyvector,wideWindow)


%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%dvector1=calculateCellPairVectorNondivision_wdiam(esequence{t},i,esequence{tj},j,anisotropyvector);
%dvector2=calculateCellPairVectorNondivision_wdiam(esequence{t},i,esequence{tk},k,anisotropyvector);

beforevectors=[];
curr=i;
curr_t=t;
pred=esequence{t}.pred(i);
pred_t=esequence{t}.pred_time(i);
count=1;
%while track doesnt end,within window and no division 
while (pred~=-1&&count<=wideWindow&&esequence{pred_t}.suc(pred,2)==-1)
    beforevectors=[beforevectors;...
        calculateCellPairVectorNondivision_wdiam(...
        esequence{pred_t},pred,esequence{curr_t},curr,anisotropyvector)...
        ];
    curr=pred;
    curr_t=pred_t;
    bkpred=pred;
    pred=esequence{pred_t}.pred(pred);
    pred_t=esequence{pred_t}.pred_time(bkpred);
    count=count+1;    
end

aftervectors=[];
curr=j;
curr_t=tj;
suc=esequence{tj}.suc(j,1);
suc_t=esequence{tj}.suc_time(j,1);
count=1;
%while track doesnt end,within window and no division 
while (suc~=-1&&count<=wideWindow&&esequence{curr_t}.suc(curr,2)==-1)
    aftervectors=[aftervectors;...
        calculateCellPairVectorNondivision_wdiam(...
        esequence{curr_t},curr,esequence{suc_t},suc,anisotropyvector)...
        ];
    curr=suc;
    curr_t=suc_t;
    suc=esequence{suc_t}.suc(curr,1);
    suc_t=esequence{suc_t}.suc_time(curr,1);
    count=count+1;    
end


aftervectors2=[];
curr=k;
curr_t=tk;
suc=esequence{tk}.suc(k,1);
suc_t=esequence{tk}.suc_time(k,1);
count=1;
%while track doesnt end,within window and no division 
while (suc~=-1&&count<=wideWindow&&esequence{curr_t}.suc(curr,2)==-1)
    aftervectors2=[aftervectors2;...
        calculateCellPairVectorNondivision_wdiam(...
        esequence{curr_t},curr,esequence{suc_t},suc,anisotropyvector)...
        ];
    curr=suc;
    curr_t=suc_t;
    suc=esequence{suc_t}.suc(curr,1);
    suc_t=esequence{suc_t}.suc_time(curr,1);
    count=count+1;    
end
if(size(beforevectors,1)>1)
    beforevectors=mean(beforevectors);
end
if(size(aftervectors,1)>1)
    aftervectors=mean(aftervectors);
end
if(size(aftervectors2,1)>1)
    aftervectors2=mean(aftervectors2);
end
if(isempty(beforevectors))
    beforevectors=ones(1,5);
end
if(isempty(aftervectors))
    aftervectors=ones(1,5);
end
if(isempty(aftervectors2))
    aftervectors2=ones(1,5);
end
data=[(beforevectors)./(aftervectors),...
    (beforevectors)./(aftervectors2),...
    (aftervectors)./(aftervectors2)];

data(isnan(data)|isinf(data))=0;
end

