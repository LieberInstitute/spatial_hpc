fname = '/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/Images/VistoSeg/VSPG/V12D07-335_D1_DAPI.tif'; 
cd /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/Images/VistoSeg/VSPG/cellpose_test
img = imread(fname);

for i = 1:8 
for j= 1:7
A{i*j} = img(y(i):y(i+1),x(j):x(j+1));
end
end

for i = 1:56                         
imwrite(A{i},sprintf('test_%d.tif',i))
end