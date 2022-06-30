offsite_dir <- .get_config_path("LOCAL_PAPALEXI_2021_DATA_DIR")
processed_dir <- paste0(offsite_dir, "processed/")
raw_dir <- paste0(offsite_dir, "raw/")

# load the Seurat data object; wrangle the most important metadata features
cat(sprintf("Loading raw data...\n"))
eccite_obj <- SeuratData::LoadData(ds = "thp1.eccite")
key_covariates <- eccite_obj@meta.data |>
  dplyr::select(orig.ident, MULTI_classification.global, Phase, percent.mito) |>
  dplyr::rename(bio_rep = MULTI_classification.global,
                phase = Phase,
                lane = orig.ident) |>
  dplyr::mutate(p_mito = percent.mito/100, percent.mito = NULL,
                bio_rep = factor(x = bio_rep, levels = c("rep1-tx", "rep3-tx", "rep4-tx"),
                               labels = c("rep_1", "rep_2", "rep_3")))

##################################
# 0. save perturbation assignments
##################################

perturbation_assignment <- eccite_obj$guide_ID |> as.character()
saveRDS(perturbation_assignment, paste0(raw_dir, "perturbation_assignments.rds"))

##################
# 1. gene modality
##################
cat(sprintf("Creating gene expression odm...\n"))
rna_modality <- eccite_obj[["RNA"]]
gene_counts <- as.matrix(rna_modality@counts)
gene_odm <- ondisc::create_ondisc_matrix_from_R_matrix(r_matrix = gene_counts,
                                                       barcodes = colnames(gene_counts),
                                                       features_df = data.frame(row.names(gene_counts)),
                                                       odm_fp = paste0(processed_dir, "gene/expression_matrix.odm"))
# add the covariates p mito, batch, and cell cycle, and perturbation
gene_odm_mod <- gene_odm |> ondisc::mutate_cell_covariates(key_covariates)
ondisc::save_odm(odm = gene_odm_mod, metadata_fp = paste0(processed_dir, "gene/metadata.rds"))

#############################
# 2. grna expression modality
#############################
cat(sprintf("Creating grna expression odm...\n"))
grna_modality <- eccite_obj[["GDO"]]
grna_counts <- as.matrix(grna_modality@counts - 1)[-1,]

grna_odm <- ondisc::create_ondisc_matrix_from_R_matrix(r_matrix = grna_counts,
                                                       barcodes = colnames(grna_counts),
                                                       features_df = data.frame(row.names(grna_counts)),
                                                       odm_fp = paste0(processed_dir, "grna_expression/count_matrix.odm"))
# extract the target and target type of each grna
grna_names <- row.names(grna_counts)
grna_targets <- gsub(pattern = "g[0-9]+", replacement = "", x = grna_names)
# replace "NT" with non-targeting for consistency
grna_targets[grna_targets == "NT"] <- "non-targeting"
# replace "PDL1" with "CD278," as the gene name in the gene data is "CD278" rather than "PDL1"
grna_targets[grna_targets == "PDL1"] <- "CD274"
# include the protein effect
protein_effect <- dplyr::recode(grna_targets, CD86 = "CD86", CD274 = "PDL1", PDCD1LG2 = "PDL2")
protein_effect[!(protein_effect %in% c("CD86", "PDL1", "PDL2"))] <- NA
# finally, get the "target_type"
target_type <- ifelse(grna_targets == "non-targeting", "non-targeting", "gene")
# add the covariates batch, cell cycle
grna_odm_mod <- grna_odm |>
  ondisc::mutate_cell_covariates(key_covariates) |>
  ondisc::mutate_feature_covariates(target = grna_targets,
                                    target_type = target_type,
                                    known_protein_effect = protein_effect)
ondisc::save_odm(odm = grna_odm_mod,
                 metadata_fp = paste0(processed_dir, "grna_expression/metadata.rds"))

#############################
# 3. grna assignment
#############################

cat(sprintf("Creating grna assignment odm...\n"))

# read information from grna expression odm
grna_expression_odm <- grna_odm_mod
cell_barcodes <- grna_expression_odm |> ondisc::get_cell_barcodes()
grna_ids <- grna_expression_odm |> ondisc::get_feature_ids()
grna_tbl <- grna_expression_odm |>
  ondisc::get_feature_covariates() |>
  dplyr::select(target, target_type, known_protein_effect)

# create the sparse gRNA assignment ODM
ondisc::convert_assign_list_to_sparse_odm(cell_barcodes = cell_barcodes,
                                          gRNA_ids = grna_ids,
                                          gRNA_assignment_list = perturbation_assignment,
                                          odm_fp = paste0(processed_dir, "grna_assignment/assignment_matrix.odm"),
                                          metadata_fp = paste0(processed_dir, "grna_assignment/metadata.rds"),
                                          features_metadata_df = grna_tbl)

#####################
# 4. protein modality
#####################
cat(sprintf("Creating protein expression odm...\n"))
protein_modality <- eccite_obj[["ADT"]]
protein_counts <- as.matrix(protein_modality@counts)
protein_odm <- ondisc::create_ondisc_matrix_from_R_matrix(r_matrix = protein_counts,
                                                          barcodes = colnames(protein_modality),
                                                          features_df = data.frame(row.names(protein_modality)),
                                                          odm_fp = paste0(processed_dir, "protein/count_matrix.odm"),
                                                          metadata_fp = paste0(processed_dir, "protein/metadata.rds"))

# add the covariates batch, cell cycle; also, add the gene that each protein is encoded by
protein_odm_mod <- protein_odm |>
  ondisc::mutate_cell_covariates(key_covariates) |>
  ondisc::mutate_feature_covariates(encoded_by = c("CD86", "CD274", "PDCD1LG2", "HAVCR2"))
ondisc::save_odm(odm = protein_odm_mod, metadata_fp = paste0(processed_dir, "protein/metadata.rds"))
