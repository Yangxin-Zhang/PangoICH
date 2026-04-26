
# R/Utils-Seurat.R

#' Origi_Data_Type H5Seurat H5ad H510X
#' @export

CreateSeurat_Pango <- function(Origi_Data_Path,
                               Spatial_Matrix_Path = NULL,
                               Origi_Data_Type = "H5Seurat") {

  if (Origi_Data_Type == "H5Seurat") {

    Seu_Obj <- CreateSeurat_From_H5_Seurat_Pango(H5Seurat_Path = Origi_Data_Path)

  }

  if (Origi_Data_Type == "H510X") {

    Seu_Obj <- CreateSeurat_From_H5_10X_Pango(H510X_Path = Origi_Data_Path,
                                              Spatial_Matrix_Path = Spatial_Matrix_Path)  %>%
      PangoICH::Standardization_Seu_Obj_Pango()

  }

  if (Origi_Data_Type == "H5ad") {

    Seu_Obj <- CreateSeurat_From_H5ad_Pango(H5ad_Path = Origi_Data_Path)

  }

  return(Seu_Obj)

}

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

  seu_obj[["percent.mt"]] <- PercentageFeatureSet(seu_obj,pattern = "^mt-")
  seu_obj[["percent.rb"]] <- PercentageFeatureSet(seu_obj,pattern = "^Rp|Rn")

  seu_obj[["log_nCount_RNA"]] <- log1p(seu_obj@meta.data$nCount_RNA) %>%
    scale()

  seu_obj[["log_nFeature_RNA"]] <- log1p(seu_obj@meta.data$nFeature_RNA) %>%
    scale()

  spatial <- PangoICH::Extract_Spatial_From_H5ad_Pango(H5ad_Path)

  seu_obj@meta.data <- bind_cols(seu_obj@meta.data,spatial)

  return(seu_obj)

}

#' @export

Add_X_To_Seurat_Meta <- function(Seu_Obj,
                                 Gene_ls,
                                 layer = "counts")
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

  for (i in 1:length(Gene_ls)) {

    gene_na <- paste(Gene_ls[i],"counts",sep = "_")

    Seu_Obj[[gene_na]] <- gene_x[,Gene_ls[i]]

  }

  if (layer == "scale.data") {

    if (length(Gene_ls) < 2) {

      fea_ls <- c(Gene_ls,rownames(Seu_Obj)[1:10]) %>%
        unique()

    } else {

      fea_ls <- Gene_ls

    }

    Seu_Obj <- ScaleData(object = Seu_Obj,
                         features = fea_ls,
                         verbose = FALSE)

    X <- Seu_Obj@assays$RNA$scale.data

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

    for (i in 1:length(Gene_ls)) {

      gene_na <- paste(Gene_ls[i],"scale.data",sep = "_")

      Seu_Obj[[gene_na]] <- gene_x[,Gene_ls[i]]

    }
  }

  return(Seu_Obj)

}

#' @export

CreateSeurat_From_H5_10X_Pango <- function(H510X_Path,
                                           Spatial_Matrix_Path)
{

  on.exit({
    rm(expr,
       seu_obj,
       spatial,
       meta_data)
    gc()
  })

  expr <- PangoICH::Extract_X_From_H5_10X_Pango(H510X_Path)

  seu_obj <- CreateSeuratObject(expr)

  seu_obj[["percent.mt"]] <- PercentageFeatureSet(seu_obj,pattern = "^mt-")
  seu_obj[["percent.rb"]] <- PercentageFeatureSet(seu_obj,pattern = "^Rp|Rn")

  seu_obj[["log_nCount_RNA"]] <- log1p(seu_obj@meta.data$nCount_RNA) %>%
    scale()

  seu_obj[["log_nFeature_RNA"]] <- log1p(seu_obj@meta.data$nFeature_RNA) %>%
    scale()

  spatial <- PangoICH::Extract_Spatial_Matrix_Pango(Spatial_Matrix_Path)

  seu_obj[["barcode"]] <- rownames(seu_obj@meta.data)
  meta_data <- as.data.table(seu_obj@meta.data)

  meta_data <- meta_data[spatial,on = "barcode"] %>%
    as.data.frame()

  rownames(meta_data) <- meta_data$barcode

  seu_obj@meta.data <- meta_data

  return(seu_obj)

}

#' @export
CreateSeurat_From_H5_Seurat_Pango <- function(H5Seurat_Path)
{

  on.exit({
    rm(Seu_Obj)
    gc()
  })

  Seu_Obj <- LoadH5Seurat(H5Seurat_Path)

  return(Seu_Obj)

}

#' @export

Save_Seu_Obj_As_H5Seurat_Pango <- function(Seu_Obj,
                                           File_Name,
                                           Path = NULL){

  if (is.null(Path)) {

    Path <- getwd() %>%
      paste(File_Name,
            sep = "/")

  }

  SaveH5Seurat(Seu_Obj,Path)

  cat("Successfully Save to",Path,sep = " ")

}
