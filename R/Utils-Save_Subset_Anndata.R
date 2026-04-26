
# R/Utils-Save_Subset_Anndata.R

#' @export

Save_Subset_Anndata_Pango <- function(H5ad_path,
                                      subset_var,
                                      var_value,
                                      subset_path = NULL){

  on.exit({
    rm(python_file_path)
    gc()
  })

  if (is.null(subset_path)) {

    H5ad_file_na <- strsplit(H5ad_path,"/") %>%
      PangoICH::Subset_PangoICH(1) %>%
      PangoICH::Subset_PangoICH(last = TRUE) %>%
      strsplit("\\.") %>%
      PangoICH::Subset_PangoICH(1) %>%
      PangoICH::Subset_PangoICH(1)

    file_dir <- strsplit(H5ad_path,"/") %>%
      PangoICH::Subset_PangoICH(1) %>%
      PangoICH::Subset_PangoICH(last = TRUE,
                                except = TRUE) %>%
      paste(collapse = "/")

    subset_path <- paste(file_dir,
                         H5ad_file_na,
                         sep = "/") %>%
      paste(subset_var,
            sep = "_") %>%
      paste(var_value,
            sep = "_") %>%
      paste("h5ad",
            sep = ".") %>%
      paste("gzip",
            sep = ".")

  }

  python_file_path <- system.file("Python",
                                  "Utils_Subset_Anndata.py",
                                  package = "PangoICH")
  source_python(python_file_path)

  Subset_Anndata_From_H5ad_Py(H5ad_path = H5ad_path,
                              subset_var = subset_var,
                              var_value = var_value,
                              subset_path = subset_path)

  cat("Successfully Save Subset Anndata to",subset_path,sep = " ")

}
