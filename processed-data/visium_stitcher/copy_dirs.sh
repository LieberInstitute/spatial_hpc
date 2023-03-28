#!/bin/bash

base_path="/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/visium_stitcher"

# Get the list of subdirectories
subdirs=$(ls "${base_path}")

# Initialize a flag to indicate whether all required files are present
all_files_present=true

# Iterate through each subdirectory
for subdir in ${subdirs}; do
  dir_path="${base_path}/${subdir}"
  file_path="${dir_path}/raw_feature_bc_matrix.h5"
  
  # Check if the raw_feature_bc_matrix.h5 file exists in the subdirectory
  if [ -f "${file_path}" ]; then
    echo "Found: ${file_path}"
  else
    echo "Not found: ${file_path}"
    all_files_present=false
  fi
done

# Print the final result
if $all_files_present; then
  echo "All raw_feature_bc_matrix.h5 files are present in the subdirectories."
else
  echo "Some raw_feature_bc_matrix.h5 files are missing in the subdirectories."
fi
