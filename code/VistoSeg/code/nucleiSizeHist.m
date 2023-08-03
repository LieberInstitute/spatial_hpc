opts = detectImportOptions('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg/code/ALLSAMPLES.txt','Delimiter','\t', 'ReadVariableNames', false);
opts.VariableNames= {'filepath','M'};
t = readtable('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg/code/ALLSAMPLES.txt',opts);

for i=1:40
fname = t.filepath{i};	
fname = [fname(1:end-4),'_nuclei.mat'];
load(fname)
stats = regionprops(mask_dark_blue, 'Area');
figure('visible','off')
histogram(struct2array(stats),'BinEdges', [50:50:2000]);
%set(gca, 'xscale','log')
ylim([0, 10000]) 
xlim([0, 2000]) 
hold on 
xline(100, 'r')
hold off
saveas(gcf,[fname(1:end-4),'_nucleiHist.png'])
end
	
for i=1:40
	fname = t.filepath{i};	
	fname = [fname(1:end-4),'_nuclei.mat'];
	load(fname)
	stats = regionprops(mask_dark_blue, 'Area');
	
	temp = [stats.Area];
	c = unique(temp); % the unique values in the A (1,2,3,4,5)    
	 for i = 1:length(c)
	   counts(i,1) = sum(temp==c(i)); % number of times each unique value is repeated
	 end
	
	plot(c,counts)
    xlim([0 1000])
    hold on 
    xline(100, 'r')
    hold off	  
end
		   

 for i=1:40
		   fname = t.filepath{i};	
		   he = imread(fname);
		   he = he(12501:13000,12501:13000,:);
		   fname = [fname(1:end-4),'_nuclei.mat'];
		   load(fname)
		   
		   s = 50;
		   mask_dark_blue = bwareaopen(mask_dark_blue, s, 8);
		   BW = 250*(double(mask_dark_blue(12501:13000,12501:13000)));
		   BW = cat(3, BW, BW, BW);
		   
		   img = [he; ones(50,500,3); BW];
		   imwrite(img, [fname(1:end-4),'_sizeThresh.png'])
		   save([fname(1:end-4),'_sizeThresh.mat'], 'mask_dark_blue','-v7.3')	
 end		   
