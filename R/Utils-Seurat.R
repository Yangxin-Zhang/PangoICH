
# R/Utils-Seurat.R

#' @export

CreateSeurat_From_H5ad_Pango <- function(H5ad_Path){

  on.exit({
    rm(expr,
       seu_obj,
       spatial)
    gc()
    })

  expr <- PangoICH::Extract_X_From_H5ad_Pango(H5ad_Path)

  seu_obj <- CreateSeuratObject(expr)

  spatial <- PangoICH::Extract_Spatial_From_H5ad_Pango(H5ad_Path)

  seu_obj@meta.data <- bind_cols(seu_obj@meta.data,spatial)

  return(seu_obj)

}

#' @export

Add_X_To_Seurat_Meta <- function(Seu_Obj,
                                 Gene_ls)
  {

  on.exit({
    rm(X,
       gene_x)
    gc()
  })

  X <- Seu_Obj@assays$RNA$counts

  if (length(Gene_ls) == 1) {

    gene_x <- X[Gene_ls,] %>%
      as.matrix()
    colnames(gene_x) <- Gene_ls

  } else {

    gene_x <- X[Gene_ls,] %>%
      as.matrix() %>%
      t() %>%
      as.data.frame()

  }

  Seu_Obj@meta.data <- bind_cols(Seu_Obj@meta.data,gene_x)

  return(Seu_Obj)

}
