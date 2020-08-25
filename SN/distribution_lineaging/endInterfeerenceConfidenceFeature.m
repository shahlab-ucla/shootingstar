function [dist ] = endInterfeerenceConfidenceFeature( esequence,t,ti,i,direction,trackingparameters )
%compute at time t pos i the vicinity of a 
%end if direction=0
%start if direction =1
%normalize this by local nn distance at i within timepoint
if(t>0&ti>0&t<length(esequence)&ti<length(esequence))
if (direction==0)
    interesting=esequence{t}.suc(:,1)==-1&~esequence{t}.delete;
else
     interesting=esequence{t}.pred==-1&~esequence{t}.delete;
end
else
    interesting=[];
end

if(isempty(interesting)||isempty(esequence{t}.selfdistance(interesting)))
  dist=[10,10];%10xnn sufficient%  dist=max(esequence{t}.selfdistance);%if there are no ends use dummy value of most distant point
else
    [mdist,mini]=min(distance_anisotropic(esequence{ti}.finalpoints(i,:)',...
        esequence{t}.finalpoints(interesting,:)',trackingparameters.anisotropyvector)); 
    dist=[distance(esequence{ti}.finalpoints(i,1:2)',esequence{t}.finalpoints(mini,1:2)'),...
        distance(esequence{ti}.finalpoints(i,3),esequence{t}.finalpoints(mini,3))*trackingparameters.anisotropyvector(3)];
    dist=dist./esequence{ti}.selfdistance(i);
end
  
end

