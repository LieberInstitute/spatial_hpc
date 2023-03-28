def process_samples(xml_files):
    for xml_file in xml_files:
        os.chdir("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/visium_stitcher/")
        
        # Load transforms
        transforms = vs.transform_finder(xml_file)
        stitched_image = xml_file.replace(".xml", ".png")
        
        # Extract sample names from the transforms dictionary
        sample_names = [key.replace("_tissue_hires_image.png", "") for key in transforms.keys()]
        
        # Create an empty list to store the AnnData objects
        adatas = []
        
        # Loop through the sample names and process data
        for sample_name in sample_names:
            adata = sc.read_visium(sample_name, count_file="raw_feature_bc_matrix.h5")
            adata.var_names_make_unique()
            adata.obs_names = [sample_name + "-" + i for i in adata.obs_names]
            adata.obs['sample'] = sample_name
            adata.uns['transform'] = transforms[sample_name + "_tissue_hires_image.png"]
            
            # Filter out the cells that are not in the tissue before stitching
            adata = adata[adata.obs['in_tissue'] == 1]
            adatas.append(adata)
        
        # Stitch together all the AnnData objects
        adata = vs.stitch(adatas, image=stitched_image)
        
        # Filter out the cells that are not in the tissue
        adata = adata[adata.obs['in_tissue'] == 1]
        
        adata_nooverlap = adata[~adata.obs["overlap"]]
        plot_file = xml_file.replace(".xml", "_capturearea.png")
        sc.pl.spatial(adata_nooverlap, color="sample", alpha=0.25, save=plot_file)
        
        adata.write("adata_{}_overlap.h5ad".format(xml_file.split('.')[0]))
        adata.write("adata_{}_no_overlap.h5ad".format(xml_file.split('.')[0]))

xml_files = [
    "081_A1B1_stitch.xml",
    "081_C1D1_stitch.xml",
    "084_stitch.xml",
    "085_stitch.xml",
    "086_A1B1_stitch.xml",
    "086_C1D1_stitch.xml",
    "2720_stitch.xml",
    "297_335_stitch.xml",
    "333_stitch.xml",
    "335_ABC_stitch.xml",
    "336_stitch.xml",
]

process_samples(xml_files)