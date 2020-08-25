function displayData(matches,matchesr,points,pointsr,stack)

anchorval=mean(mean(median(stack)))*4;

FP=points(find(matches==-1),:);
FN=pointsr(find(matchesr==-1),:);
pmatches=points(find(matches~=-1),:);



%detected centers
%visualize for heck of it
data=zeros(size(stack));

if (length(FP)~=0)
data(sub2ind(size(data),FP(:,2),FP(:,1),FP(:,3)))=.75*anchorval;
end
if (length(FN)~=0)
    data(sub2ind(size(data),round(FN(:,2)),round(FN(:,1)),round(FN(:,3))))=1*anchorval;
end

pmatches=round(pmatches);
data(sub2ind(size(data),pmatches(:,2),pmatches(:,1),pmatches(:,3)))=.5*anchorval;
data=imdilate(data,strel('disk',2));

%all labeled data are in yellow small
data(sub2ind(size(data),round(pointsr(:,2)),round(pointsr(:,1)),round(pointsr(:,3))))=.6*anchorval;
data=imdilate(data,strel('disk',2));

datatest=max(data,stack);

%{
maxval=max(max(max(datatest)));
minval=min(min(min(datatest)));
Xsave=uint8(floor(255*(datatest-minval)/(maxval-minval)));
for i = 1:size(datatest,3)
imwrite(Xsave(:,:,i),'Test_labeled_output_looser.tif','tif', 'Compression', 'none', 'WriteMode', 'append');
end
%}

nii=make_nii(datatest);
view_nii(nii);
