 img = load('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/Images/VistoSeg/VSPG/V12D07-332_B1.mat');
[x,y] = size(img.DAPI);
 DAPI = cat(3,zeros(x,y), zeros(x,y), img.DAPI);
 neuron = cat(3,zeros(x,y), img.Alexa_488, zeros(x,y));
 TMEM = cat(3,img.Alexa_555, img.Alexa_555, zeros(x,y));
 GFAP = cat(3, img.Alexa_594, zeros(x,y), zeros(x,y));
 OLIG = cat(3,img.Alexa_647, zeros(x,y), img.Alexa_647);


temp = cat(3, img.Alexa_555+img.Alexa_594+img.Alexa_647, img.Alexa_488+img.Alexa_555, img.DAPI+img.Alexa_647);

imwrite(temp, 'temp.png')