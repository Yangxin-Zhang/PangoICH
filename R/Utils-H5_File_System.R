
# R/Utils_File_System.R

#' @export

Extract_X_From_H5ad_Pango <- function(H5ad_Path)
  {

  on.exit({
    rm(adata,
       X,
       var_names,
       obs_names)
    gc()
  })

  adata <- HDF5Array::H5ADMatrix(H5ad_Path)

  X <- as(adata,"dgCMatrix")

  var_names <- h5read(H5ad_Path, "/var/_index/values") %>%
    as.character()
  obs_names <- h5read(H5ad_Path, "/obs/_index/values") %>%
    as.character()

  rownames(X) <- var_names
  colnames(X) <- obs_names

  return(X)

}

#' @export

Extract_Spatial_From_H5ad_Pango <- function(H5ad_Path)
  {

  on.exit({
    rm(raw_spatial,
       spatial,
       spatial_dt)
    gc()
  })

  raw_spatial <- h5read(H5ad_Path, "/obsm/raw_spatial")
  spatial <- h5read(H5ad_Path,"/obsm/spatial")

  spatial_dt <- data.table(raw_x = raw_spatial[1,],
                           raw_y = raw_spatial[2,],
                           x = spatial[1,],
                           y = spatial[2,])

  return(spatial_dt)

}
