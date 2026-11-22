# Installation of Bioconductor packages (with BiocManager)
# Anna Giannakogeorgou

# Bioconductor annotation/data packages exist on bioconda only as stubs whose
# post-link download pixi does not run. Run via `pixi run install-bioc`.

bioc_packages <- c("org.Hs.eg.db")

# Stop if this is not the pixi R library (e.g. when run from a system R)
if (!any(grepl(".pixi", .libPaths(), fixed = TRUE))) {
  stop("R library is not inside .pixi/ - run via `pixi run install-bioc` ",
       "or switch Positron to the pixi R interpreter.")
}

missing_packages <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) == 0) {
  message("All Bioconductor packages already installed: ", paste(bioc_packages, collapse = ", "))
} else {
  message("Installing: ", paste(missing_packages, collapse = ", "))
  BiocManager::install(missing_packages, update = FALSE, ask = FALSE)

  # Stop if an install failed, so dependent pixi tasks don't run
  still_missing <- missing_packages[!vapply(missing_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still_missing) > 0) stop("Installation failed for: ", paste(still_missing, collapse = ", "))
}

# Print versions (not recorded in pixi.lock)
for (pkg in bioc_packages) message(pkg, " ", as.character(packageVersion(pkg)))