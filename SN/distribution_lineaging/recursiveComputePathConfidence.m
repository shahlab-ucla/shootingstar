function esequence=recursiveComputePathConfidence(esequence,t,i,con)
if(isfield(esequence{t},'linkconfidences'))
newcon=con*esequence{t}.linkconfidences(i);
esequence{t}.path_confidence(i)=newcon;
if(esequence{t}.suc(i,1)~=-1)
    esequence=recursiveComputePathConfidence(...
        esequence,esequence{t}.suc_time(i,1),esequence{t}.suc(i,1),newcon);
    if(esequence{t}.suc(i,2)~=-1)
            esequence=recursiveComputePathConfidence(...
        esequence,esequence{t}.suc_time(i,2),esequence{t}.suc(i,2),newcon);
    end
end
end