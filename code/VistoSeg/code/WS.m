%function WS(fname)
	opts = detectImportOptions('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg/code/ALLSAMPLES.txt','Delimiter','\t', 'ReadVariableNames', false);
	opts.VariableNames= {'filepath','M'};
	t = readtable('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg/code/ALLSAMPLES.txt',opts);

	for i = 1:40
		
	   fname = t.filepath{i};	
	   he = imread(fname);
	   he = he(12501:13000,12501:13000,:);
	   load([fname(1:end-4),'_mask.mat'])
		   seg1 = mask{t.M(i)};
	   load([fname(1:end-4),'_nuclei.mat'])
		   seg2 = mask_dark_blue;	
	   load([fname(1:end-4),'_sizeThresh.mat'])
		   seg3 = mask_dark_blue;

	   D = -bwdist(~mask_dark_blue);
	   mask = imextendedmin(D,1);
	   D2 = imimposemin(D,mask);
	   Ld2 = watershed(D2);
	   bw3 = mask_dark_blue;
	   bw3(Ld2 == 0) = 0;
	   mask_dark_blue = bw3;

	   save([fname(1:end-4),'_nuclei_WS.mat'],'mask_dark_blue','-v7.3')
	
   seg1 = 250*(double(seg1(12501:13000,12501:13000)));
   seg1 = cat(3, seg1, seg1, seg1);
   seg2 = 250*(double(seg2(12501:13000,12501:13000)));
   seg2 = cat(3, seg2, seg2, seg2);
   seg3 = 250*(double(seg3(12501:13000,12501:13000)));
   seg3 = cat(3, seg3, seg3, seg3);
   seg4 = 250*(double(bw3(12501:13000,12501:13000)));
   seg4 = cat(3, seg4, seg4, seg4);
   
   img = [he; 250*ones(10,500,3); seg1; 250*ones(10,500,3); seg2; 250*ones(10,500,3); seg3; 250*ones(10,500,3); seg4];
   imwrite(img, [fname(1:end-4),'_segFinal.png'])
	
end

