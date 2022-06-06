offsite_dir <- .get_config_path("LOCAL_PAPALEXI_2021_DATA_DIR")
processed_dir <- paste0(offsite_dir, "processed/")

# load the Seurat data object; wrangle the most important metadata features
SeuratData::InstallData("thp1.eccite")
eccite_obj <- SeuratData::LoadData(ds = "thp1.eccite")
key_covariates <- eccite_obj@meta.data |>
  dplyr::select(MULTI_classification.global, Phase, percent.mito, guide_ID) |>
  dplyr::rename("batch" = "MULTI_classification.global",
                "phase" = "Phase",
                "perturbation" = "guide_ID") |>
  dplyr::mutate(p_mito = percent.mito/100, percent.mito = NULL,
                batch = factor(x = batch, levels = c("rep1-tx", "rep2-tx", "rep3-tx"),
                               labels = c("rep_1", "rep_2", "rep_3")))


##################
# 1. gene modality
##################
rna_modality <- eccite_obj[["RNA"]]
gene_counts <- as.matrix(rna_modality@counts)
gene_odm <- ondisc::create_ondisc_matrix_from_R_matrix(r_matrix = gene_counts,
                                                       barcodes = colnames(gene_counts),
                                                       features_df = data.frame(row.names(gene_counts)),
                                                       odm_fp = paste0(processed_dir, "gene/expression_matrix.odm"))
# add the covariates p mito, batch, and cell cycle, and perturbation
gene_odm_mod <- gene_odm |> ondisc::mutate_cell_covariates(key_covariates)
ondisc::save_odm(odm = gene_odm_mod, metadata_fp = paste0(processed_dir, "gene/metadata.rds"))


##################
# 2. gRNA modality
##################
gRNA_modality <- eccite_obj[["GDO"]]
gRNA_counts <- as.matrix(gRNA_modality@counts - 1)[-1,]

# gRNA_odm <- ondisc::read_odm(odm_fp = paste0(processed_dir, "gRNA/count_matrix.odm"),
#                              metadata_fp = paste0(processed_dir, "gRNA/metadata.rds"))

gRNA_odm <- ondisc::create_ondisc_matrix_from_R_matrix(r_matrix = gRNA_counts,
                                                       barcodes = colnames(gRNA_counts),
                                                       features_df = data.frame(row.names(gRNA_counts)),
                                                       odm_fp = paste0(processed_dir, "grna/count_matrix.odm"))
# extract the target and target type of each gRNA
gRNA_names <- row.names(gRNA_counts)
gRNA_targets <- gsub(pattern = "g[0-9]+", replacement = "", x = gRNA_names)
# replace "NT" with non-targeting for consistency
gRNA_targets[gRNA_targets == "NT"] <- "non-targeting"
# replace "PDL1" with "CD278," as the gene name in the gene data is "CD278" rather than "PDL1"
gRNA_targets[gRNA_targets == "PDL1"] <- "CD274"
# include the protein effect
protein_effect <- dplyr::recode(gRNA_targets, CD86 = "CD86", CD274 = "PDL1", PDCD1LG2 = "PDL2")
protein_effect[!(protein_effect %in% c("CD86", "PDL1", "PDL2"))] <- NA
# finally, get the "target_type"
target_type <- ifelse(gRNA_targets == "non-targeting", "non-targeting", "gene")
# add the covariates batch, cell cycle
gRNA_odm_mod <- gRNA_odm |>
  ondisc::mutate_cell_covariates(key_covariates) |>
  ondisc::mutate_feature_covariates(target = gRNA_targets,
                                    target_type = target_type,
                                    known_protein_effect = protein_effect)
ondisc::save_odm(odm = gRNA_odm_mod,
                 metadata_fp = paste0(processed_dir, "grna/metadata.rds"))


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

# protein_odm <- ondisc::read_odm(odm_fp = paste0(processed_dir, "protein/count_matrix.odm"),
# metadata_fp = paste0(processed_dir, "protein/metadata.rds"))

# add the covariates batch, cell cycle; also, add the gene that each protein is encoded by
protein_odm_mod <- protein_odm |>
  ondisc::mutate_cell_covariates(key_covariates) |>
  ondisc::mutate_feature_covariates(encoded_by = c("CD86", "CD274", "PDCD1LG2", "HAVCR2"))
ondisc::save_odm(odm = protein_odm_mod, metadata_fp = paste0(processed_dir, "protein/metadata.rds"))

