
dir='L:\santella\NIH_data\neural_for_seg\20121111_DCR553_for_segmentation\RL_Decon\'
figure
view(3);
startt=100;
endt=120;
m=moviein(endt-startt+1);
for time=startt:endt
    
data=loadsimpleStackTIFF([dir,'Decon_',num2str(time),'.tif']);
threshold=20;

[f,v]=isosurface(data,threshold);
%[f2,v2,c2]=isocaps(data,threshold);
p = patch('Faces',f,'Vertices',v);
        set(p,'SpecularColorReflectance',.1,'FaceColor',[.5,.5,.5],'BackFaceLighting','reverselit','EdgeColor','none','FaceAlpha',.5);
s=(size(data));
axis([0,s(2),0,s(1),0,s(3)]);
       % camlight('headlight')
        lighting phong
        frame=getframe(gcf);
          m(time-startt+1)=frame;
clf

end
%m(1)=m(2);
diru='L:\santella\NIH_data\neural_for_seg\20121111_DCR553_for_segmentation\'

movieName='isotest4.avi';
movie2avi(m,[diru,movieName],'fps',6);



%ignore for now and j
dir='L:\santella\NIH_data\neural_for_seg\20121111_DCR553_for_segmentation\RL_Decon\'
time=420;
data=loadsimpleStackTIFF([dir,'Decon_',num2str(time),'_manualsegv2.tif']);

hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(data), hy, 'replicate');
Ix = imfilter(double(data), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);
baseimage=gradmag;
baseimage=double(data);
baseimage=1/baseimage;%inverted main image
baseimage(isinf(baseimage))=1;
datatemp=smooth3(data,'box',[5,5,5]);

%datatemp(datatemp<150)=0;
foreground=data;%foreground markers >200
foreground(datatemp>=275)=1;
foreground(datatemp<275)=0;

background=data;%foreground markers >200
background(datatemp>=30)=1;
background(datatemp<30)=0;
D = bwdist(background);
DL = watershed(D);
bgm = DL == 0;
bothmarks=bgm|foreground;

%mergedimage=imimposemin(baseimage,bothmarks);
mergedimage=imimposemin(baseimage,foreground,6);

pdata=1/mergedimage;
pdata(pdata==-0)=max(max(max(pdata)));

threshold=75;
[f,v]=isosurface(pdata,threshold);
%[f2,v2,c2]=isocaps(data,threshold);
p = patch('Faces',f,'Vertices',v);
        set(p,'SpecularColorReflectance',.1,'FaceColor',[.5,.5,.5],'BackFaceLighting','reverselit','EdgeColor','none','FaceAlpha',.75);
s=(size(data));
axis([0,s(2),0,s(1),0,s(3)]);
       % camlight('headlight')
        lighting phong

%L=watershed(mergedimage);

%[P,J]=regionGrowing(data,[60,195,166],100);

figure;
point=[190,66,166;205,127,210;295,101,167;308,107,146];
color=[1,0,0;0,1,0;0,0,1;1,0,1];
threshold=[20,20,75,75];
obs=[];
for i=1:4
mask=zeros(size(data));
p=point(i,:);
mask(p(1),p(2),p(3))=1;
temp=baseimage;
if(i<3)
temp=smooth3(temp,'box',[3,3,3]);
end
mergedimage=imimposemin(temp,mask,6);
pdata=1/mergedimage;
pdata(pdata==-0)=max(max(max(pdata)));
%threshold=75;

[f,v]=isosurface(pdata,threshold(i));
%[f2,v2,c2]=isocaps(data,threshold);
obs(i) = patch('Faces',f,'Vertices',v);
        set(obs(i),'SpecularColorReflectance',.1,'FaceColor',color(i,:),'BackFaceLighting','reverselit','EdgeColor','none','FaceAlpha',.75);
s=(size(data));
axis([0,s(2),0,s(1),0,s(3)]);
       % camlight('headlight')
        lighting phong
end




datatemp=loadsimpleStackTIFF([dir,'Decon_',num2str(time),'_manualsegv4.tif']);

%for i=1:size(data,3)
%    datatemp(:,:,i)=medfilt2(datatemp(:,:,i),[21,21]);
%end

figure
test2=imresize(datatemp,.25);
for i=1:4
test2=ordfilt3d(test2,14);
end
%test2=imopen(test2,strel('disk',31));
nhood=ones(11,11,11);
center=size(nhood)./2;
for x=1:size(nhood,1)
    for y=1:size(nhood,2)
        for z=1:size(nhood,3)
            if(distance(center',[x,y,z]')>size(nhood,1)/2)
                nhood(x,y,z)=0;
            end
        end
    end
end
test2=imopen(test2,strel('arbitrary',nhood));

%{
test3(test2>2)=3;
nii=make_nii(test3);
view_nii(nii)
%}
test2(:,:,230:end)=0;
smooth3(test2,'box',[11,11,11]);
[f,v]=isosurface(test2,1);
v(:,1:2)=v(:,1:2)*4;
%[f2,v2,c2]=isocaps(data,threshold);
obs(5) = patch('Faces',f,'Vertices',v);
%set(obs(5),'SpecularColorReflectance',.1,'FaceColor',[0,0,0],'BackFaceLighting','reverselit','EdgeColor','none','FaceAlpha',.05);
set(obs(5),'SpecularColorReflectance',.1,'FaceColor',[.9,.9,.9],'BackFaceLighting','reverselit','EdgeColor','none','FaceAlpha',.3);
diam=2;
[xp,yp,zp]=sphere(20);
xp=xp*diam;
yp=yp*diam;
zp=zp*diam;
count=0;
hold on
while count<650
    %for x=1:20:size(data,1)
    %    for y=1:20:size(data,2)
    %        for z=1:20:size(data,3)
    x=max(1,floor(rand(1)*size(data,1)));
    y=max(1,floor(rand(1)*size(data,2)));
    z=max(1,floor(rand(1)*size(data,3)));
    xtest=max(1,floor(x/4));
    ytest=max(1,floor(y/4));
    if (test2(xtest,ytest,z)>=1)
        obs(end+1)= surface(yp+y,xp+x,zp+z,'FaceColor',[.5,.5,.5],'EdgeColor','none','FaceAlpha',1);
        count=count+1;
    end
    %        end
    %    end
end
%{
for i=7:length(obs)
set(obs(i),'SpecularColorReflectance',.1,'FaceColor',[.5,.5,.5]','BackFaceLighting','reverselit','EdgeColor','none','FaceAlpha',1);
end
%}
zdir = [0 1 0];
%rotate(h12,zdir,25)

[Z,Y]=meshgrid(0:100:500,0:100:500);
X=ones(size(X)).*160;
testpatch=surf(X,Y,Z);
set(testpatch,'FaceColor','white','BackFaceLighting','lit','AmbientStrength',1,'DiffuseStrength',0,'EdgeColor','none','FaceAlpha',.2);
%obs(end)=testpatch;


 view(-90,0);
 h=camlight('right')
m=moviein(360/4);
axis off
set(gcf, 'color', 'white');
c=1
for degree=1:4:360
    %{
    if(degree>180)
        rot=-1*(360-degree);
        view(90,rot);
    else
    %}
        rot=degree;
       % view(-90,rot);
    %end
    for j=1:length(obs)
    rotate(obs(j),zdir,4);
    end
    camlight(h,'right');
    frame=getframe(gcf);
     m(c)=frame;
     c=c+1
end

diru='L:\santella\NIH_data\neural_for_seg\20121111_DCR553_for_segmentation\'

movieName='isoteststill_withcells.avi';
movie2avi(m,[diru,movieName],'fps',6);



