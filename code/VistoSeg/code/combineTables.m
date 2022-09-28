disp('loading data')
tic
fname = '/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/Images/VistoSeg/Capture_areas/*.csv';

myfiles = dir(fname);
length(myfiles)

T = cell2table(cell(0,16), 'VariableNames', {'Area', 'Centroid_1', 'Centroid_2', 'BoundingBox_1', 'BoundingBox_2', 'BoundingBox_3', 'BoundingBox_4', 'MajorAxisLength', 'MinorAxisLength', 'Eccentricity', 'Circularity', 'Perimeter', 'WeightedCentroid_1', 'WeightedCentroid_2', 'MeanIntensity', 'name'});

for i = 1:length(myfiles)
mask_name = fullfile(myfiles(i).folder,myfiles(i).name);
T1 = readtable(mask_name);
T1.name = repmat(myfiles(i).name,height(T1),1);
T = [T;T1];
end

writetable(T,fullfile(myfiles(1).folder,'combined_refineVNS_metric_WS.csv'))

