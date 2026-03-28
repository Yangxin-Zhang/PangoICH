
# R/Spatial_Graph.R

#' @export

Gene_Distribution_Graph_Pango <- function(Seu_Obj,
                                          Gene)
  {

  on.exit({
    rm(seu_obj,
       P_X,
       gene_x,
       graph_plot)
    gc()
  })

  seu_obj <- Add_X_To_Seurat_Meta(Seu_Obj = Seu_Obj,
                                  Gene_ls = Gene)

  P_X <- seu_obj@meta.data %>%
    as.data.table()

  gene_x <- P_X[,..Gene]
  P_X[,plotting_gene := gene_x]

  graph_plot <- ggplot() +
    geom_point(data = P_X[plotting_gene == 0],
               mapping = aes(x=-x,
                             y=-y),
               size = 0.01,
               color = "grey90") +
    geom_point(data = P_X[plotting_gene != 0],
               mapping = aes(x=-x,
                             y=-y,
                             colour = plotting_gene),
               size = 0.01) +
    scale_colour_gradientn(colors = c("#FEF4E8", "#FED9A6", "#FEB24C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026"),
                           values = scales::rescale(seq(from = min(P_X[plotting_gene != 0,plotting_gene]),
                                                        to = max(P_X[plotting_gene != 0,plotting_gene]),
                                                        length.out = 7),
                                                    to = c(0,1)),
                           limits = c(min(P_X[plotting_gene != 0,plotting_gene]),
                                      max(P_X[plotting_gene != 0,plotting_gene]))) +
    theme(panel.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          aspect.ratio = 0.8,
          legend.location = "none")

  return(graph_plot)

}
