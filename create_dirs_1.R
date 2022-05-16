offsite_dir <- .get_config_path("LOCAL_PAPALEXI_2021_DATA_DIR")

# create raw and processed subdirectories
raw_dir <- paste0(offsite_dir, "raw/")
processed_dir <- paste0(offsite_dir, "processed/")
processed_subdirs <- paste0(processed_dir, c("gene", "grna", "protein"), "/")
counts_dir <- paste0(raw_dir, "count")

if (!dir.exists(raw_dir)) dir.create(raw_dir)
if (!dir.exists(processed_dir)) dir.create(processed_dir)
if (!dir.exists(counts_dir)) dir.create(counts_dir)

for (dir in processed_subdirs) {
  if (!dir.exists(dir)) dir.create(dir)
}
