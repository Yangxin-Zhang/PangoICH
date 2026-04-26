
# R/Utils-Picture_Related_Function.R

#' @export

Convert_ggplot_As_RasterGrob <- function(ggplot_Obj,
                                         width = 8,
                                         height = 8,
                                         units = "in",
                                         dpi = 600,
                                         device = "png"){

  tmp_file <- tempfile(fileext = ".png")

  ggsave(tmp_file,
         plot = ggplot_Obj,
         device = device,
         width = width,
         height = height,
         units = units,
         dpi = dpi)

  img <- image_read(tmp_file) %>%
    rasterGrob()

  return(img)

}

#' @export

Annotation_Label <- function(label,
                             size = 12,
                             angle = 0,
                             fontcolour = "black",
                             fontface = "bold",
                             background_fill = "white",
                             background_colour = "white"){

  ann_lab <- ggplot() +
    annotate(geom = "text",
             label = label,
             x = 0.5,
             y = 0.5,
             size = size*0.3527778,
             fontface = fontface,
             colour = fontcolour,
             angle = angle,)+
    theme_void() +
    theme(panel.background = element_rect(fill = background_fill,
                                          colour = background_colour),
          plot.background = element_rect(fill = "white",
                                         colour = "white"))

  return(ggplotGrob(ann_lab))

}
