Papalexi (2021) Data Documentation
================
Tim Barry
April 6, 2022

# Overview

This repository contains code to import and process the Papalexi 2021
data. Papalexi 2021 conducted many experiments, including flow
cytommetry, bulk RNA-seq, CITE-seq, pilot ECCITE-seq, and at-scale
ECCITE-seq experiments. Their GEO repository contains data only from the
final three experiments. This repository processes data only from the
final of these experiments (i.e., the at-scale ECCITE-seq experiment).

The top-level offsite `papalexi_2021` data directory contains two
subdirectories: `processed` and `raw`. The `raw` subdirectory contains
the raw data, and the `processed` subdirectory contains the processed
data in ODM format. There are four subdirectories in the `processed`
subdirectory, each corresponding to a diffrent modality: gene, gRNA
(expressions), gRNA (assignments), and protein.

``` r
# top-level directory
papalexi_dir <-.get_config_path("LOCAL_PAPALEXI_2021_DATA_DIR")
# processed subdirectory
processed_dir <- paste0(papalexi_dir, "processed/")
# list files in processed directory
list.files(processed_dir)
```

    ## [1] "gene"            "grna_assignment" "grna_expression" "protein"

We examine these modalities one at a time. Data are available on 20,729
cells, all of which have been stimulated by the IFN-gamma chemical to
upregulate expression of the protein PD-L1.

# gRNA (expressions)

The data contain 110 gRNAs. (The data originally contained 111 gRNAs, as
reported in the paper, but one of the gRNAs exhibited zero expression
across all cells, so I removed it.) The gRNAs either target genes or are
non-targeting. Each gene is targeted by 3-4 gRNAs. Meanwhile, there are
9 non-targeting gRNAs. The non-targeting gRNAs are labeled NTg1 – NTg10,
with NTg6 apparently missing.

``` r
# load gRNA odm; save to gRNA_odm
gRNA_odm_exp <- ondisc::read_odm(odm_fp = paste0(processed_dir, "grna_expression/count_matrix.odm"),
                             metadata_fp = paste0(processed_dir, "grna_expression/metadata.rds"))
```

The `feature_covariates` matrix of `gRNA_odm_exp` contains columns
`target`, `target_type`, and `known_protein_effect`. `target` lists the
gene that each gRNA targets. `target_type`, meanwhile, gives the type of
target of each gRNA (either `gene` or `non-targeting`). Finally,
`known_protein_effect` gives the protein that is encoded by a given
targeted gene (for genes whose encoded proteins were sequenced; NA
otherwise).

``` r
# print several feature covariates of gRNA modality
gRNA_odm_exp |>
  ondisc::get_feature_covariates() |>
  head() |> dplyr::select(target, target_type, known_protein_effect)
```

    ##         target target_type known_protein_effect
    ## CUL3g1    CUL3        gene                 <NA>
    ## CUL3g2    CUL3        gene                 <NA>
    ## CUL3g3    CUL3        gene                 <NA>
    ## CMTM6g1  CMTM6        gene                 <NA>
    ## CMTM6g2  CMTM6        gene                 <NA>
    ## CMTM6g3  CMTM6        gene                 <NA>

``` r
# create a table of the gRNA targets
gRNA_odm_exp |>
  ondisc::get_feature_covariates() |>
  dplyr::pull(target) |>
  table()
```

    ## 
    ##          ATF2          BRD4          CAV1         CD274          CD86 
    ##             4             4             4             3             4 
    ##         CMTM6          CUL3          ETV7        IFNGR1        IFNGR2 
    ##             3             3             4             4             4 
    ##          IRF1          IRF7          JAK2        MARCH8           MYC 
    ##             4             4             4             4             4 
    ##        NFKBIA non-targeting      PDCD1LG2        POU2F2         SMAD4 
    ##             4             9             4             4             4 
    ##          SPI1         STAT1         STAT2         STAT3        STAT5A 
    ##             4             4             4             4             4 
    ##      TNFRSF14        UBE2L6 
    ##             4             4

``` r
# create a table of the target types
gRNA_odm_exp |>
  ondisc::get_feature_covariates() |>
  dplyr::pull(target_type) |>
  table()
```

    ## 
    ##          gene non-targeting 
    ##           101             9

``` r
# create a table of the known protein effects
gRNA_odm_exp |> ondisc::get_feature_covariates() |>
  dplyr::pull(known_protein_effect) |>
  stringr::str_replace_na() |>
  table()
```

    ## 
    ## CD86   NA PDL1 PDL2 
    ##    4   99    3    4

# Gene

The gene expression data contain 18,649 genes.

``` r
# load gene odm
gene_odm <- ondisc::read_odm(odm_fp = paste0(processed_dir, "gene/expression_matrix.odm"),
                             metadata_fp = paste0(processed_dir, "gene/metadata.rds"))
```

The gene modality is linked to the gRNA modality by the `target` column
of the `feature_covariates` matrix of `gRNA_odm_exp`, as described
above. The `target` column is a subset of the `feature_ids` of
`gene_odm`.

``` r
gRNA_target_genes <- gRNA_odm_exp |> ondisc::get_feature_covariates() |>
  dplyr::filter(target_type == "gene") |>
  dplyr::pull(target)

all(gRNA_target_genes %in% ondisc::get_feature_ids(gene_odm))
```

    ## [1] TRUE

# Protein

The protein data contain four proteins: CD86, PDL1, PDL2, and CD366.
These proteins are encoded by the genes *CD86*, *CD274*, *PDCD1LG2*, and
*HAVCR2* respectively. CD366 also goes by the names HAVCR2 and TIM-3.

``` r
# load protein odm
protein_odm <- ondisc::read_odm(odm_fp = paste0(processed_dir, "protein/count_matrix.odm"),
                                metadata_fp = paste0(processed_dir, "protein/metadata.rds"))
```

The protein data are linked both to the gene and gRNA data. First, the
`feature_covariates` matrix of `protein_odm` lists the gene that encodes
each protein in the `encoded_by` column; this column is a subset of
`feature_ids` of `gene_odm`.

``` r
# feature covariates of proteins
protein_odm |> ondisc::get_feature_covariates()
```

    ##       mean_expression coef_of_variation n_nonzero encoded_by
    ## CD86        104.42404         1.0157685     20729       CD86
    ## PDL1        193.51247         0.9997360     20728      CD274
    ## PDL2         61.34787         0.7152157     20728   PDCD1LG2
    ## CD366        43.70930         0.9490350     20724     HAVCR2

``` r
all((protein_odm |>
    ondisc::get_feature_covariates() |>
    dplyr::pull()) %in% ondisc::get_feature_ids(gene_odm))
```

    ## [1] TRUE

Second, as described above, the `feature_covariates` data frame of
`gRNA_odm_exp` contains a column called `known_protein_effect`, which,
for genes that encode one of the sequenced proteins, lists the encoded
protein (`NA` otherwise). The column `known_protein_effect` is a subset
of the `feature_ids` of the `protein_odm`.

``` r
all((gRNA_odm_exp |> ondisc::get_feature_covariates() |>
    na.omit() |>
  dplyr::pull(known_protein_effect)) %in% ondisc::get_feature_ids(protein_odm))
```

    ## [1] TRUE

The genes that encode proteins CD86, PDL1, and PDL2 were all perturbed.
The gene that encodes CD366 (i.e., *HAVCR2*), by contrast, was *not*
perturbed, as far as I can tell. The protein CD366 is not mentioned
anywhere in the manuscript, so it probably is safe to ignore this
protein.

# gRNA (assignments)

The gRNA assignment matrix is a binary version of the gRNA expression
matrix, where the cell-to-gRNA assignments were called by the original
authors.

``` r
gRNA_odm_assign <- ondisc::read_odm(odm_fp = paste0(processed_dir, "grna_assignment/assignment_matrix.odm"),
                                    metadata_fp = paste0(processed_dir, "grna_assignment/metadata.rds"))
gRNA_odm_assign[[80:90, 1:5]]
```

    ## 11 x 5 sparse Matrix of class "lgCMatrix"
    ##                
    ##  [1,] . . . . .
    ##  [2,] . . . . .
    ##  [3,] . . . . .
    ##  [4,] . . . . .
    ##  [5,] . . . . .
    ##  [6,] . . | . .
    ##  [7,] . . . . .
    ##  [8,] . . . . .
    ##  [9,] . . . . .
    ## [10,] | . . . .
    ## [11,] . . . . .

# Cell-specific covariates across modalities

The `gene_odm` has covariates `n_nonzero`, `n_umis`, `lane`, `bio_rep`,
`phase`, and `p_mito`.

``` r
gene_odm |> ondisc::get_cell_covariates() |> head()
```

    ##                     n_nonzero n_umis  lane bio_rep phase     p_mito
    ## l1_AAACCTGAGCCAGAAC      3942  17207 Lane1   rep_1    G1 0.02295577
    ## l1_AAACCTGAGTGGACGT      2948   9506 Lane1   rep_1    G1 0.04512939
    ## l1_AAACCTGCATGAGCGA      4258  15256 Lane1   rep_1    G1 0.04116413
    ## l1_AAACCTGTCTTGTCAT      1780   5135 Lane1   rep_1    G1 0.05491723
    ## l1_AAACGGGAGAACAACT      2671   9673 Lane1   rep_1    G1 0.03359868
    ## l1_AAACGGGAGACAGAGA      3918  14941 Lane1   rep_1    G1 0.03379961

The other modalities (gRNA, protein) have the same covariates.
