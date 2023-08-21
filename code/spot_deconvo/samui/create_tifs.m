
myfiles = dir('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/spot_deconvo/samui/*1.tif');


fname='';
BW = imread(fname);


imwrite(mat2gray(I1.(O{1})),[fname(1:end-4),'_A1.tif'])