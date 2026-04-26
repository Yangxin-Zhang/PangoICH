
# R/Spatial_Graph.R

#' @export

Gene_Distribution_Graph_Pango <- function(Seu_Obj,
                                          Gene,
                                          layer = "counts",
                                          Dispose_Main = FALSE)
  {

  on.exit({
    rm(seu_obj,
       P_X,
       gene_x,
       graph_plot)
    gc()
  })

  plot_title <- Gene
  seu_obj <- Add_X_To_Seurat_Meta(Seu_Obj = Seu_Obj,
                                  Gene_ls = Gene,
                                  layer = layer) %>%
    .add_background_genes(Background_Gene = c("Hbb-bt","Hbb-bs","Hba-a2"))

  P_X <- seu_obj@meta.data %>%
    as.data.table()

  Gene_counts <- paste(Gene,"counts",sep = "_")

  if (layer == "counts") {

    Gene <- paste(Gene,"counts",sep = "_")

  }

  if (layer == "scale.data") {

    Gene <- paste(Gene,"scale.data",sep = "_")

  }

  gene_x <- P_X[,..Gene]
  P_X[,plotting_gene := gene_x]

  gene_x_counts <- P_X[,..Gene_counts]
  P_X[,plotting_gene_counts := gene_x_counts]

  if (!Dispose_Main) {

    graph_plot <- ggplot() +
      geom_point(data = P_X[plotting_gene_counts == 0],
                 mapping = aes(x=-x,
                               y=-y),
                 size = 0.01,
                 color = "grey90") +
      geom_point(data = P_X[background_genes != 0],
                 mapping = aes(x=-x,
                               y=-y),
                 size = 0.01,
                 color = "white")+
      geom_point(data = P_X[plotting_gene_counts != 0],
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
      labs(title = plot_title) +
      theme(panel.background = element_rect(colour = "black",
                                            fill = "white",
                                            linetype = "solid"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            aspect.ratio = 1,
            legend.position = "none",
            legend.title = element_blank(),
            plot.margin = margin(0,0,0,0),
            plot.title = element_text(family = "Arial",
                                      size = 12,
                                      colour = "black",
                                      vjust = 0.5,
                                      hjust = 0.5,
                                      face = "bold",
                                      margin = margin(5,0,3,0)),
            plot.background = element_rect(fill = "white",
                                           colour = "white"))
  } else {

    graph_plot_text <- ggplot(data = P_X,
                              mapping = aes(x=-x,
                                            y=-y)) +
            labs(title = plot_title) +
      theme(panel.background = element_rect(colour = "white",
                                            fill = "white",
                                            linetype = "solid"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            aspect.ratio = 1,
            legend.position = "none",
            legend.title = element_blank(),
            plot.margin = margin(0,0,0,0),
            plot.title = element_text(family = "Arial",
                                      size = 12,
                                      colour = "black",
                                      vjust = 0.5,
                                      hjust = 0.5,
                                      face = "bold",
                                      margin = margin(0,0,0,0)),
            plot.background = element_rect(fill = "white",
                                           colour = "white"))

    graph_plot_main <- ggplot() +
      scale_x_continuous(expand = c(0.01, 0.01)) +
      scale_y_continuous(expand = c(0.01, 0.01)) +
      geom_point(data = P_X[plotting_gene_counts == 0],
                 mapping = aes(x=-x,
                               y=-y),
                 size = 0.01,
                 color = "grey90") +
      geom_point(data = P_X[background_genes != 0],
                 mapping = aes(x=-x,
                               y=-y),
                 size = 0.01,
                 color = "white")+
      geom_point(data = P_X[plotting_gene_counts != 0],
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
      labs(title = plot_title) +
      theme(panel.background = element_rect(colour = "black",
                                            fill = "white",
                                            linetype = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            aspect.ratio = 1,
            legend.position = "none",
            legend.title = element_blank(),
            plot.margin = margin(0,0,0,0),
            plot.title = element_blank(),
            plot.background = element_rect(fill = "white",
                                           colour = "white"))

    graph_plot <- list(Main = PangoICH::Convert_ggplot_As_RasterGrob(graph_plot_main),
                       Other = ggplotGrob(graph_plot_text))


  }


  return(graph_plot)

}

#' @export

Meta_Data_Graph_Pango <- function(Seu_Obj,
                                  Aim_Var,
                                  GreyScale_Image = TRUE){

  on.exit({
    rm(P_X,
       meta_data,
       graph_plot)
    gc()
  })

  P_X <- Seu_Obj@meta.data %>%
    as.data.table()

  meta_data <- P_X[,..Aim_Var]
  P_X[,plotting_data := meta_data]

  if (GreyScale_Image) {

    scale_color <- scale_colour_gradientn(colours = gray(seq(0, 1, length.out = 10)),
                                          values = scales::rescale(seq(from = min(P_X[plotting_data != 0,plotting_data]),
                                                                       to = max(P_X[plotting_data != 0,plotting_data]),
                                                                       length.out = 7),
                                                                   to = c(0,1)),
                                          limits = c(min(P_X[plotting_data != 0,plotting_data]),
                                                     max(P_X[plotting_data != 0,plotting_data])))

  } else {

    scale_color <- scale_colour_gradientn(colors = c("#FEF4E8", "#FED9A6", "#FEB24C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026"),
                                          values = scales::rescale(seq(from = min(P_X[plotting_data != 0,plotting_data]),
                                                                       to = max(P_X[plotting_data != 0,plotting_data]),
                                                                       length.out = 7),
                                                                   to = c(0,1)),
                                          limits = c(min(P_X[plotting_data != 0,plotting_data]),
                                                     max(P_X[plotting_data != 0,plotting_data])))

  }

  graph_plot <- ggplot() +
    geom_point(data = P_X,
               mapping = aes(x=-x,
                             y=-y,
                             colour = plotting_data),
               size = 0.01) +
    scale_color +
    theme(panel.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          aspect.ratio = 1,
          legend.title = element_blank(),
          legend.position = "none")

  return(graph_plot)

}

#'

.add_background_genes <- function(Seu_Obj,
                                  Background_Gene)
  {

  X <- Seu_Obj@assays$RNA$counts

  if (length(Background_Gene) == 1) {

    gene_expr <- X[Background_Gene,] %>%
      as.matrix()
    colnames(gene_expr) <- "background_genes"

    Seu_Obj[["background_genes"]] <- gene_expr[,"background_genes"]

  }
  if (length(Background_Gene) > 1) {

    gene_expr <- X[Background_Gene,] %>%
      as.matrix() %>%
      t() %>%
      as.data.frame() %>%
      Matrix::rowSums() %>%
      as.matrix()

    colnames(gene_expr) <- "background_genes"

    Seu_Obj[["background_genes"]] <- gene_expr[,"background_genes"]

  }
  return(Seu_Obj)
}






