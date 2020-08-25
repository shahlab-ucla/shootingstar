function diskSet=calculateSliceGFP(Xorig,diskSet)
global parameters;
diskSet.GFPsums=zeros(1,length(diskSet.xydetdiameters));
diskSet.diskArea=zeros(1,length(diskSet.xydetdiameters));
diskSet.diskMax=zeros(1,length(diskSet.xydetdiameters));
raysort2s=[1,13,5,9,3,12,7,15,2,14,6,10,4,11,8,16];
for i=1:length(diskSet.xydetdiameters)
    
    
    
    good=logical(diskSet.xpositions(:,i)~=0|diskSet.ypositions(:,i)~=0);
    
    sortedgood=good(raysort2s);
    sortedx=diskSet.xpositions(raysort2s,i);
    sortedy=diskSet.ypositions(raysort2s,i);
    
    
    
    %need to order and filter points
    %
    pointsmask=poly2mask(sortedx,sortedy,size(Xorig,1),size(Xorig,2));
    points=find(pointsmask);
    diskSet.diskArea(i)=length(points);
    slice=Xorig(:,:,diskSet.xymax(i,3));
    diskSet.GFPsums(i)=sum(slice(points));
    %compute median value of non b
    if(~isfield(parameters,'polarBackgroundTrheshold'))
       % points=find(pointsmask&slice>2200);
       points=find(pointsmask);
    else
        points=find(pointsmask&slice>parameters.polarBackgroundThreshold);
    end
    if(isempty(points))
        diskSet.diskMax(i)=nan;
    else
    diskSet.diskMax(i)=prctile(slice(points),95);
    end
    %median(slice(points));%
end
diskSet;