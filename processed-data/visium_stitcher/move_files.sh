#!/bin/bash

base_dir="/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/visium_stitcher"

for subdir in "${base_dir}"/[BV]*; do
    if [[ -d "${subdir}" ]]; then
        # Create the "spatial" directory if it doesn't exist
        mkdir -p "${subdir}/spatial"

        # Move all files except "raw_feature_bc_matrix.h5" to the "spatial" directory
        find "${subdir}" -maxdepth 1 -type f ! -name "raw_feature_bc_matrix.h5" -exec mv {} "${subdir}/spatial" \;
    fi
done
