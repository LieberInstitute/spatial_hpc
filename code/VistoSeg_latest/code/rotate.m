

fname = '/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/raw-data/Images/VSPG/V12D07-332_rerun/V12D07-332_D1.mat';
I4 = load(fname);
O{1} = 'DAPI'; O{2} = 'Alexa_488'; O{3} = 'Alexa_555'; O{4} = 'Alexa_594'; O{5} = 'Alexa_647'; O{6} = 'Autofluorescence'; 
O = O';
I4.DAPI = imrotate(I4.DAPI,180);
I4.Alexa_488 = imrotate(I4.Alexa_488,180);
I4.Alexa_555 = imrotate(I4.Alexa_555,180);
I4.Alexa_594 = imrotate(I4.Alexa_594,180);
I4.Alexa_647 = imrotate(I4.Alexa_647,180);
I4.Autofluorescence = imrotate(I4.Autofluorescence,180);

save(fname,'-struct','I4', '-v7.3');

imwrite(mat2gray(I4.(O{1})),[fname(1:end-4),'.tif'])
for i = 2:numel(O)
imwrite(mat2gray(I4.(O{i})),[fname(1:end-4),'.tif'],'writemode', 'append')
end
