copyfile /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/Images/VistoSeg/VSPG/*1.tif /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/spot_deconvo/groundTruth/02_samui_manual_annotation/
%scp /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/Images/VistoSeg/VSPG/*1.tif /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/spot_deconvo/groundTruth/02_samui_manual_annotation/


%% subtract autofluorescence from TMEM
myfiles = dir('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/Images/VistoSeg/VSPG/*1.mat');
for i = 1:numel(myfiles)
    fname=fullfile(myfiles(i).folder, myfiles(i).name);
    load(fname);
    Alexa_555 = Alexa_555-Autofluorescence;
    fname1=fullfile('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/spot_deconvo/groundTruth/02_samui_manual_annotation/', myfiles(i).name);
    imwrite(mat2gray(DAPI),[fname1(1:end-4),'.tif'])
    imwrite(mat2gray(Alexa_488),[fname1(1:end-4),'.tif'],'writemode', 'append')
    imwrite(mat2gray(Alexa_555),[fname1(1:end-4),'.tif'],'writemode', 'append')
    imwrite(mat2gray(Alexa_594),[fname1(1:end-4),'.tif'],'writemode', 'append')
    imwrite(mat2gray(Alexa_647),[fname1(1:end-4),'.tif'],'writemode', 'append')
    imwrite(mat2gray(Autofluorescence),[fname1(1:end-4),'.tif'],'writemode', 'append')
end

%% attach segmented image to raw image
myfiles = dir('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/spot_deconvo/groundTruth/01_cellpose/final_masks/*.png');

for i = 1:numel(myfiles)
fname=fullfile(myfiles(i).folder, myfiles(i).name);
BW = imread(fname);
BW = label2rgb(BW,'hsv','w','shuffle');
Or = '/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/spot_deconvo/groundTruth/02_samui_manual_annotation/';
imwrite(BW,fullfile(Or, [myfiles(i).name(1:end-18),'.tif']),"WriteMode","append")
disp(i)
end

