# required R packages: Seurat, SeuratData, ondisc, dplyr,

#$ -l m_mem_free=16G

hpcc=$(hostname | grep "hpcc" | wc -l)
if [[ hpcc ]]
then
  module load R/R-4.1.2
fi
Rscript create_dirs_1.R
Rscript download_raw_2.R
Rscript create_odms_3.R
