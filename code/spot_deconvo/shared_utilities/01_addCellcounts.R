#-------------------------------------------------------------------------------
#   Add cell counts to HE spatial object
#-------------------------------------------------------------------------------
# HE_id_path <- here("processed-data", "spot_deconvo", "HE_ID_table.csv")
# HE_counts_path <- here("processed-data", "rerun_spaceranger", "{sample_id}", "outs", "spatial", "tissue_spot_counts.csv")

# id_table <- read.csv(HE_id_path)
# 
# spe_HE$count <- NA
# for (sample_id in unique(spe_HE$sample_id)) {
#   #   Correctly determine the path for the cell counts for this sample, then
#   #   read in
#   long_id <- id_table[match(sample_id, id_table$short_id), "long_id"]
#   this_path <- sub("{sample_id}", long_id, HE_counts_path, fixed = TRUE)
#   cell_counts <- read.csv(this_path)
#   
#   #   All spots in the object should have counts
#   stopifnot(
#     all(
#       colnames(spe_nonIF[, spe_nonIF$sample_id == sample_id]) %in%
#         cell_counts$barcode
#     )
#   )
#   
#   #   Line up the rows of 'cell_counts' with the sample-subsetted SPE object
#   cell_counts <- cell_counts[
#     match(
#       colnames(spe_nonIF[, spe_nonIF$sample_id == sample_id]),
#       cell_counts$barcode
#     ),
#   ]
#   
#   #   Add this sample's counts to the SPE object
#   spe_nonIF$count[spe_nonIF$sample_id == sample_id] <- cell_counts$Nmask_dark_blue
# }
# 
# #   Ensure counts were read in for all spots in the object
# if (any(is.na(spe_nonIF$count))) {
#   stop("Did not find cell counts for all non-IF spots.")
# }

