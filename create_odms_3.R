offsite_dir <- .get_config_path("LOCAL_PAPALEXI_2021_DATA_DIR")
processed_dir <- paste0(offsite_dir, "processed/")

# load the Seurat data object; wrangle the most important metadata features
eccite_obj <- LoadData(ds = "thp1.eccite")
key_covariates <- eccite_obj@meta.data |>
  dplyr::select(MULTI_classification.global, Phase, percent.mito) |>
  dplyr::rename("batch" = "MULTI_classification.global",
                "phase" = "Phase") |>
  dplyr::mutate(p_mito = percent.mito/100, percent.mito = NULL) |>
  dplyr::mutate(batch = factor(x = batch, levels = c("rep1-tx", "rep2-tx", "rep3-tx"), labels = c("rep_1", "rep_2", "rep_3")))


##################
# 1. gene modality
##################
rna_modality <- eccite_obj[["RNA"]]
gene_counts <- as.matrix(rna_modality@counts)
gene_odm <- ondisc::create_ondisc_matrix_from_R_matrix(r_matrix = gene_counts,
                                                       barcodes = colnames(gene_counts),
                                                       features_df = data.frame(row.names(gene_counts)),
                                                       odm_fp = paste0(processed_dir, "gene/expression_matrix.odm"))
# add the covariates p mito, batch, and cell cycle
gene_odm_mod <- gene_odm |> ondisc::mutate_cell_covariates(key_covariates)
ondisc::save_odm(odm = gene_odm_mod, metadata_fp = paste0(processed_dir, "gene/metadata.rds"))


##################
# 2. gRNA modality
##################
gRNA_modality <- eccite_obj[["GDO"]]
gRNA_counts <- as.matrix(gRNA_modality@counts - 1)
gRNA_odm <- ondisc::create_ondisc_matrix_from_R_matrix(r_matrix = gRNA_counts,
                                                       barcodes = colnames(gRNA_counts),
                                                       features_df = data.frame(row.names(gRNA_counts)),
                                                       odm_fp = paste0(processed_dir, "gRNA/count_matrix.odm"))
# add the covariates batch, cell cycle
gRNA_odm_mod <- gRNA_odm |> ondisc::mutate_cell_covariates(phase = key_covariates$phase,
                                                           batch = key_covariates$batch)
ondisc::save_odm(odm = gRNA_odm_mod, metadata_fp = paste0(processed_dir, "gRNA/metadata.rds"))


##################
# protein modality
##################
protein_modality <- eccite_obj[["ADT"]]
protein_counts <- as.matrix(protein_modality@counts)
protein_odm <- ondisc::create_ondisc_matrix_from_R_matrix(r_matrix = protein_counts,
                                                          barcodes = colnames(protein_modality),
                                                          features_df = data.frame(row.names(protein_modality)),
                                                          odm_fp = paste0(processed_dir, "protein/count_matrix.odm"),
                                                          metadata_fp = paste0(processed_dir, "protein/metadata.rds"))
# add the covariates batch, cell cycle
protein_odm_mod <- protein_odm |> ondisc::mutate_cell_covariates(phase = key_covariates$phase,
                                                                 batch = key_covariates$batch)
ondisc::save_odm(odm = gRNA_odm_mod, metadata_fp = paste0(processed_dir, "gRNA/metadata.rds"))
