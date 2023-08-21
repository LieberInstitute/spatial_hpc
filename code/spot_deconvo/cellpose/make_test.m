myfiles = dir('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/images/VistoSeg/VSPG/*.mat');

for i=1:numel(myfiles)
    fname = fullfile(myfiles(i).folder,myfiles(i).name);
    load(fname)
    imwrite(DAPI,[fname(1:end-4),'_DAPI.tif'])
end


%% make test samples with smaller size
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

%% run cellpose interactively
%% move files to respective folders

movefile /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/Images/VistoSeg/VSPG/cellpose_test/models/* /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/spot_deconvo/cellpose/maddy_models/
movefile /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/Images/VistoSeg/VSPG/cellpose_test/* /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/spot_deconvo/cellpose/maddy_models/
movefile /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/Images/VistoSeg/VSPG/*.npy /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/spot_deconvo/cellpose/final_masks
movefile /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/Images/VistoSeg/VSPG/*.png /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/spot_deconvo/cellpose/final_masks
movefile /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/Images/VistoSeg/VSPG/*DAPI.tif /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/spot_deconvo/cellpose/final_masks

