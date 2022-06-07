offsite_dir <- .get_config_path("LOCAL_PAPALEXI_2021_DATA_DIR")
processed_dir <- paste0(offsite_dir, "processed/")
raw_dir <- paste0(offsite_dir, "raw/")

# load the gRNA-to-cell assignments
assignments <- readRDS(paste0(raw_dir, "perturbation_assignments.rds")) |> as.list()


# load the cell barcodes and gRNA IDs
grna_expression_odm <- ondisc::read_odm(odm_fp = paste0(processed_dir, "grna_expression/count_matrix.odm"),
                                        metadata_fp = paste0(processed_dir, "grna_expression/metadata.rds"))
cell_barcodes <- grna_expression_odm |> ondisc::get_cell_barcodes()
gRNA_ids <- grna_expression_odm |> ondisc::get_feature_ids()
gRNA_tbl <- grna_expression_odm |>
  ondisc::get_feature_covariates() |>
  dplyr::select(target, target_type, known_protein_effect)

# create the sparse ODM
ondisc::convert_assign_list_to_sparse_odm(cell_barcodes = cell_barcodes,
                                          gRNA_ids = gRNA_ids,
                                          gRNA_assignment_list = assignments,
                                          odm_fp = paste0(processed_dir, "grna_assignment/assignment_matrix.odm"),
                                          metadata_fp = paste0(processed_dir, "grna_assignment/metadata.rds"),
                                          features_metadata_df = gRNA_tbl)
