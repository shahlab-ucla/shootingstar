function stack=writeSimpleStackTiff(name,image)



depth=size(image,3);
for j=1:depth
    imwrite(image(:,:,j),name,'tif','Compression','none','WriteMode','append');

end

