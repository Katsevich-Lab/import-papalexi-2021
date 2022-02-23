offsite_dir <- .get_config_path("LOCAL_PAPALEXI_2019_DATA_DIR")

# create raw and processed subdirectories
raw_dir <- paste0(offsite_dir, "raw/")
processed_dir <- paste0(offsite_dir, "processed/")
counts_dir <- paste0(raw_dir, "count")

if (!dir.exists(raw_dir)) dir.create(raw_dir)
if (!dir.exists(processed_dir)) dir.create(processed_dir)
if (!dir.exists(counts_dir)) dir.create(counts_dir)
