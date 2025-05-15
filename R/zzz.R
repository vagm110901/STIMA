.onLoad <- function(libname, pkgname) {
  if (packageVersion("Seurat") != "5.0.2" | packageVersion("SeuratObject") != "5.0.2") {
    remove.packages("Seurat")
    remove.packages("SeuratObject")
    remotes::install_version("Seurat", version = "5.0.2")
    remotes::install_version("SeuratObject", version = "5.0.2")
  }
}
