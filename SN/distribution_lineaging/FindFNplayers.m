function [ startconflictplayers,startbacktrace,endconflictplayers,endforwardtrace ] = FindFNplayers( start_t,start_i,end_t,end_i,esequence )
%given putative FN endpoints pulls in other players via conflicting nn
%network
%if result of back/forward tracing from conflict points is equivalent to
%start then returns [] in backtrace variables
%players in conflict points (or just the input back again if there isnt a
%conflict
if(esequence{start_t}.fNN(start_i)~=-1)
    startconflictplayers=esequence{start_t+1}.predecessor_suitors{esequence{start_t}.fNN(start_i)};
    
      %filter out deteted suitors since they cannot be players
    startconflictplayers=startconflictplayers(~esequence{start_t}.delete(startconflictplayers));
    
    %ensure that in weird configurations the linked predecessor of the
    %conflict point is always a player even when it is not a nn suitor
    if(esequence{start_t+1}.pred(esequence{start_t}.fNN(start_i))~=-1)
        %the actual pred of collapse point should always be a player
        %if it is a gap link then the whole thing is a dirty case and we
        %give up finding players
        if(esequence{start_t+1}.pred_time(esequence{start_t}.fNN(start_i))~=start_t)
            startconflictplayers=[];
        else
            startconflictplayers=unique([startconflictplayers;...
            esequence{start_t+1}.pred(esequence{start_t}.fNN(start_i))]);
        end
    end
else
    startconflictplayers=[];
end
%if(esequence{end_t}.pred(end_i)~=-1)
%    endconflictplayers=esequence{end_t-1}.suc(esequence{end_t}.pred(end_i),:)';
%else
    %endconflictplayers=end_i; %if not the division nothing
    if(esequence{end_t}.bNN(end_i)~=-1)
      endconflictplayers=esequence{end_t-1}.sucessor_suitors{esequence{end_t}.bNN(end_i)};
    
     
      %filter out deteted suitors since they cannot be players
      endconflictplayers=endconflictplayers(~esequence{end_t}.delete(endconflictplayers));

    else
         endconflictplayers=[];
    end
%end
%endconflictplayers=unique([endconflictplayers;esequence{end_t-1}.sucessor_suitors{esequence{end_t}.bNN(end_i)}]);

%trace forward to purported FN 'beginning'
i=start_i;
i=esequence{start_t}.fNN(i);%step past link gap and follow suc rest  of way
%if the bnn is deleted then we dont even try to trace
if(i==-1||isempty(i)||esequence{start_t+1}.delete(i))
    i=-1;
end
if(i~=-1)
for t=start_t+1:end_t-1
    %not clean if: fn or if there is a division and this is not the
    %bifurcation we're currently examining
    
      %not sure why had to add these 5/15/2014
    if (i==-1)
        break
    end
    if((esequence{t}.suc_time(i,1)~=t+1||...
            (~(t==end_t-1&i==esequence{end_t}.pred(end_i))&esequence{t}.suc(i,2)~=-1) ))
        i=-1;
        break
    else
        %i=esequence{t}.fNN(i);
        i=esequence{t}.suc(i,1);
    end
end
end
endforwardtrace=i;

%trace back to purported FN 'start'
i=end_i;
i=esequence{end_t}.bNN(i);

%if the bnn is deleted then we dont even try to trace
if(i==-1||isempty(i)||esequence{end_t-1}.delete(i))
    i=-1;
end
%trace backward from div to start not clean if have fn or division
%interveene
if(i~=-1)
for t=end_t-1:-1:start_t+1
    %not sure why had to add these 5/15/2014
    if(i==-1)
        break
    end
    %same check but since we know the only time w'ere in the tp of bif we
    %are in bif can use simpler rule
    if(esequence{t}.pred_time(i,1)~=t-1||(t~=end_t-1&esequence{t}.suc(i,2)~=-1))
        i=-1;
        break
    else
        i=esequence{t}.pred(i);
        %i=esequence{t}.bNN(i);
    end
end
end
startbacktrace=i;

%if traced forward/back are one of conflicting then wipe it out to signal
%this
if (max(startbacktrace==startconflictplayers))
    startbacktrace=[];
end

if (max(endforwardtrace==endconflictplayers))
    endforwardtrace=[];
end

end

