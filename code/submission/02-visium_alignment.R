library(tidyverse)
library(here)
library(readxl)
library(sessioninfo)

he_sr_params_path = here(
    'code', 'submission', 'spaceranger_parameters.txt'
)
if_sr_params_path = here(
    "code", "submission", "spaceranger_IF_parameters.txt"
)
pd_path = here('processed-data', 'nda-submission', 'pheno_data.csv')

out_dir = here('processed-data', 'nda-submission')

col_names = c(
    "subjectkey", "src_subject_id", "interview_date", "interview_age", "sex",
    "alignment_file", "loupe_version"
)

#   Read in H&E spaceranger parameters
he_sr_params = read_tsv(
        he_sr_params_path, col_names = FALSE, show_col_types = FALSE
    ) |>
    mutate(
        donor = str_extract(X1, 'Br[0-9]{4}'),
        src_subject_id = sprintf('%s_%s', X2, X3),
        loupe_version = "6.2.0"
    ) |>
    rename(alignment_file = X5) |>
    select(donor, src_subject_id, alignment_file, loupe_version)

#   Read in IF spaceranger parameters
if_sr_params = read_tsv(
        if_sr_params_path, col_names = FALSE, show_col_types = FALSE
    ) |>
    mutate(
        donor = str_extract(X1, 'Br[0-9]{4}'),
        src_subject_id = sprintf('%s_%s', X2, X3),
        loupe_version = "6.4.1"
    ) |>
    rename(alignment_file = X5) |>
    select(donor, src_subject_id, alignment_file, loupe_version)

#   Read in general donor-level info
pd = read_csv(pd_path, show_col_types = FALSE) |>
    rename(donor = src_subject_id)


#   Combine all info, grab the necessary columns, and export to CSV
rbind(he_sr_params, if_sr_params) |>
    left_join(pd, by = 'donor') |>
    select(all_of(col_names)) |>
    write_csv(file.path(out_dir, 'visium_alignment.csv'))

session_info()
