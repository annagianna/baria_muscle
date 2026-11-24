# Download MitoCarta3.0 files 
# Anna Giannakogeorgou

mito_dir <- "data/raw_data/MitoCarta"
dir.create(mito_dir, recursive = TRUE, showWarnings = FALSE)

base_url <- "https://sites.broadinstitute.org/files/mitocarta/files/documents/"
mito_files <- c(
  "Human.MitoCarta3.0.xls" = "human.mitocarta3.0.xls",
  "Human.MitoPathways3.0.gmx" = "human.mitopathways3.0.gmx_.txt"
)

for (local_name in names(mito_files)) {
  dest <- file.path(mito_dir, local_name)
  if (file.exists(dest)) {
    message("Already present: ", dest)
    next
  }
  message("Downloading: ", local_name)
  download.file(paste0(base_url, mito_files[[local_name]]), dest, mode = "wb")
  if (!file.exists(dest) || file.size(dest) == 0) stop("Download failed: ", local_name)
}