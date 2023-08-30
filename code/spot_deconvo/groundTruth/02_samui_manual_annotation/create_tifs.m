copyfile /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/Images/VistoSeg/VSPG/*1.tif /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/spot_deconvo/groundTruth/02_samui_manual_annotation/
%scp /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/Images/VistoSeg/VSPG/*1.tif /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/spot_deconvo/groundTruth/02_samui_manual_annotation/


myfiles = dir('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/spot_deconvo/groundTruth/01_cellpose/final_masks/*.png');

for i = 1:numel(myfiles)
fname=fullfile(myfiles(i).folder, myfiles(i).name);
BW = imread(fname);
BW = label2rgb(BW,'hsv','w','shuffle');
Or = '/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/spot_deconvo/groundTruth/02_samui_manual_annotation/';
imwrite(BW,fullfile(Or, [myfiles(i).name(1:end-18),'.tif']),"WriteMode","append")
disp(i)
end

