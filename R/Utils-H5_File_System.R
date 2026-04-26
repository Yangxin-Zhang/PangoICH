
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

#' @export

Extract_X_From_H5_10X_Pango <- function(H510X_Path)
  {

  on.exit({
    rm(X)
    gc()
  })

  X <- Read10X_h5(H510X_Path)

  return(X)

}

#' @export

Extract_Spatial_Matrix_Pango <- function(Spatial_Matrix_Path)
  {

  on.exit({
    rm(sp_mat,
       spatial_dt)
    gc()
  })

  sp_mat <- Read10X_Coordinates(filename = Spatial_Matrix_Path,
                                filter.matrix = TRUE)

  spatial_dt <- data.table(barcode = rownames(sp_mat),
                           raw_x = sp_mat[,"imagerow"],
                           raw_y = sp_mat[,"imagecol"],
                           x = sp_mat[,"row"],
                           y = sp_mat[,"col"])

  return(spatial_dt)

}
