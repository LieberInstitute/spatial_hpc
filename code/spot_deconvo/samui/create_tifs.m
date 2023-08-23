copyfile /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/Images/VistoSeg/VSPG/*1.tif /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/spot_deconvo/samui/
%scp /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/Images/VistoSeg/VSPG/*1.tif /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/spot_deconvo/samui/


myfiles = dir('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/spot_deconvo/cellpose/final_masks/*.png');

for i = 1:numel(myfiles)
fname=fullfile(myfiles(i).folder, myfiles(i).name);
BW = imread(fname);
BW = label2rgb(BW,'hsv','w','shuffle');
Or = '/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/samui/';
imwrite(BW,fullfile(Or, [myfiles(i).name(1:end-18),'.tif']),"WriteMode","append")
disp(i)
end

