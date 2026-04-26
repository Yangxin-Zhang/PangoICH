
# R/Utils_Seurat_Analysis_Pipeline.R

#' @export

Standardization_Seu_Obj_Pango <- function(Seu_Obj,
                                          normalization.method = "LogNormalize",
                                          scale.factor = 10000)
  {

  on.exit({
    rm()
    gc()
  })

  Seu_Obj <- subset(Seu_Obj,subset = (nFeature_RNA > 5 &
                                        percent.mt <= 20)) %>%
    NormalizeData(normalization.method = normalization.method,
                  scale.factor = scale.factor,
                  verbose = FALSE)

  return(Seu_Obj)

}
