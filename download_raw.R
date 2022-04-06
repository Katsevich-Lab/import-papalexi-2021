offsite_dir <- .get_config_path("LOCAL_PAPALEXI_2021_DATA_DIR")
raw_dir <- paste0(offsite_dir, "raw/")
count_dir <- paste0(raw_dir, "count/")

geo_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE153056&format=file&file="
filenames <- c("GSE153056_CITE_metadata.tsv.gz",
               "GSE153056_ECCITE_Arrayed_metadata.tsv.gz",
               "GSE153056_ECCITE_metadata.tsv.gz",
               "GSE153056_RAW.tar")
sources <- paste0(geo_url, filenames)
sources[4] <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE153056&format=file"
dests <- paste0(raw_dir, filenames)

# download files to raw directory
for (i in 1:length(filenames)) {
  print(paste0("Downloading ", filenames[i]))
  download.file(url = sources[i], destfile = dests[i])
  download.file(url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE153056&format=file", destfile = dests[i])
}

# untar/unzip the count data
untar(dests[4], exdir = count_dir)
fs_to_unzip <- paste0(count_dir, list.files(count_dir))
for (f in fs_to_unzip) R.utils::gunzip(f)
R.utils::gunzip(paste0(count_dir, list.files(count_dir)))
file.remove(dests[4])

# unzip the other data
for (f in dests[1:3]) {
  R.utils::gunzip(f)
}

# finally, get the vector of stimulated cells for the SueratData object
library(Seurat)
library(SeuratData)
InstallData(ds = "thp1.eccite")
eccite <- LoadData(ds = "thp1.eccite")
stim_cells <- Cells(eccite)
saveRDS(stim_cells, paste0(raw_dir, "/count/stim_cells.rds"))
