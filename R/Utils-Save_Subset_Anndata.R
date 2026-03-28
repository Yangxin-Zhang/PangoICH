
# R/Utils-Save_Subset_Anndata.R

#' @export

Save_Subset_Anndata_Pango <- function(H5ad_path,
                                      subset_var,
                                      var_value,
                                      subset_path){

  on.exit({
    rm(python_file_path)
    gc()
  })

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
