
# R/Utils-Shiny_PangoICH.R

#' @export

Gene_Distribution_Map_Shiny_PangoICH <- function(H5_Database_Path,
                                                 Spatial_Matrix_Path = NULL){

  gene_distribution_map_shiny <- shinyApp(ui = .gene_distribution_map_shiny_UI(H5_Database_Path = H5_Database_Path),
                                          server = .gene_distribution_map_shiny_Server(H5_Database_Path = H5_Database_Path,
                                                                                       Spatial_Matrix_Path = Spatial_Matrix_Path))

  return(gene_distribution_map_shiny)

}
#' shiny UI

.gene_distribution_map_shiny_UI <- function(H5_Database_Path) {

  ui <- fluidPage(sidebarLayout(sidebarPanel(width = 3,
                                             selectInput(inputId = "H5_Database",
                                                         label = "H5 Database",
                                                         choices = c("",names(H5_Database_Path)),
                                                         multiple = FALSE),
                                             selectInput(inputId = "H5_File_Type",
                                                         label = "H5 File Type",
                                                         choices = c("",
                                                                     "H5ad",
                                                                     "H510X",
                                                                     "H5Seurat"),
                                                         multiple = FALSE),
                                             uiOutput("Spatial_Matrix_Path"),
                                             selectInput(inputId = "Graph_Type",
                                                         label = "Graph Type",
                                                         choices = c("",
                                                                     "gene_distribution",
                                                                     "meta_data"),
                                                         multiple = FALSE),
                                             uiOutput("Plotting_Data_Layer"),
                                             uiOutput("Gene_ID"),
                                             uiOutput("Meta_Data_Var"),
                                             uiOutput("Save_Plot"),
                                             uiOutput("Gene_Distribution_Graph"),
                                             uiOutput("Save_As_RData")
  ),
  mainPanel(width = 9,
            plotOutput("Gene_Distribution_Map",
                       height = "600px"))))

  return(ui)

}

#' shiny server

.gene_distribution_map_shiny_Server <- function(H5_Database_Path,
                                                Spatial_Matrix_Path) {

  server <- function(input,output,session)
  {

    output$Spatial_Matrix_Path <- renderUI({

      req(input$H5_File_Type)

      if (input$H5_File_Type == "H510X") {

        selectInput(inputId = "Spatial_Matrix",
                    label = "Spatial Matrix",
                    choices = c("",names(Spatial_Matrix_Path)),
                    multiple = FALSE)

      }
    })

    output$Plotting_Data_Layer <- renderUI({

      req(input$Graph_Type)

      if (input$Graph_Type == "gene_distribution") {

        selectInput(inputId = "Plotting_Data_Layer",
                    label = "Plotting Data Layer",
                    choices = c("",
                                "counts",
                                "scale.data"),
                    multiple = FALSE)

      }
    })

    output$Gene_ID <- renderUI({

      req(input$Graph_Type)

      if (input$Graph_Type == "gene_distribution") {

        selectInput(inputId = "Gene_ID",
                    label = "Gene ID",
                    choices = NULL,
                    selected = NULL,
                    multiple = FALSE)

      }
    })

    output$Gene_Distribution_Graph <- renderUI({

      req(input$Graph_Type)

      if (input$Graph_Type == "gene_distribution") {

        selectInput(inputId = "Gene_Distribution_Graph",
                    label = "Gene Distribution Graph",
                    choices = NULL,
                    selected = NULL,
                    multiple = TRUE)

      }
    })

    output$Meta_Data_Var <- renderUI({

      req(input$Graph_Type)

      if (input$Graph_Type == "meta_data") {

        selectInput(inputId = "Meta_Data_Var",
                    label = "Meta Data Var",
                    choices = c("",Meta_Data()),
                    multiple = FALSE)

      }
    })



    observe({

      req(Gene_IDs())

      if (input$Graph_Type == "gene_distribution") {

        updateSelectizeInput(
          session,
          "Gene_ID",
          choices = Gene_IDs(),
          server = TRUE)

        updateSelectizeInput(
          session,
          "Gene_Distribution_Graph",
          choices = Gene_IDs(),
          server = TRUE)

      }
    })

    observeEvent(input$Save_Plot, {

      if (length(input$Gene_ID) !=0 ) {

        file_na <- paste(input$Gene_ID,"png",sep = ".")

      } else {

        file_na <- "test.png"

      }

      ggsave(plot = Plotting_Graph(),
             filename = file_na,
             device = "png",
             path = "/home/youngxin/Documents/PangoICH/Graph",
             width = 200,
             height = 200,
             units = "mm",
             dpi = 600)

    })

    observeEvent(input$Save_As_RData, {

      req(input$Gene_Distribution_Graph)

      genes <- input$Gene_Distribution_Graph

      for (i in 1:length(genes)) {

        file_na_main <- paste(genes[i],"main",sep = "_") %>%
          paste("RData",sep = ".")

        file_na_other <- paste(genes[i],"other",sep = "_") %>%
          paste("RData",sep = ".")

        graph <- PangoICH::Gene_Distribution_Graph_Pango(Seu_Obj = Seu_Obj(),
                                                         Gene = genes[i],
                                                         layer = input$Plotting_Data_Layer,
                                                         Dispose_Main = TRUE)

        saveRDS(graph$Main,file = paste("/home/youngxin/Documents/PangoICH/Graph",file_na_main,sep = "/"))
        saveRDS(graph$Other,file = paste("/home/youngxin/Documents/PangoICH/Graph",file_na_other,sep = "/"))

        print("Save Successfully!")

      }
    })

    output$Gene_Distribution_Map <- renderPlot({

      Plotting_Graph()

    })

    output$Save_Plot <- renderUI({

      actionButton(inputId = "Save_Plot",
                   label = "Save Plot")

    })

    output$Save_As_RData <- renderUI({

      req(input$Graph_Type)
      if (input$Graph_Type == "gene_distribution") {

        actionButton(inputId = "Save_As_RData",
                     label = "Save As RData")

      }
    })

    # Reactive Object Module
    ##

    Seu_Obj <- reactive({

      req(input$H5_File_Type)

      if (input$H5_File_Type == "H510X") {

        req(input$Spatial_Matrix)

        seu_obj <- PangoICH::CreateSeurat_Pango(H5_Database_Path[input$H5_Database],
                                                Origi_Data_Type = input$H5_File_Type,
                                                Spatial_Matrix_Path = Spatial_Matrix_Path[input$Spatial_Matrix])

      } else  {

        seu_obj <- PangoICH::CreateSeurat_Pango(H5_Database_Path[input$H5_Database],
                                                Origi_Data_Type = input$H5_File_Type)

      }

      return(seu_obj)

    })

    Meta_Data <- reactive({

      req(input$Graph_Type,Seu_Obj())

      if (input$Graph_Type == "meta_data") {

        meta_data_vars <- Seu_Obj()@meta.data %>%
          colnames()

        return(meta_data_vars)

      }

    })

    Gene_IDs <- reactive({

      req(input$Graph_Type,Seu_Obj())

      if (input$Graph_Type == "gene_distribution") {

        genes <- rownames(Seu_Obj())

        return(genes)

      }
    })

    Plotting_Graph <- reactive({

      req(input$Graph_Type)

      if (input$Graph_Type == "gene_distribution") {

        graph <- PangoICH::Gene_Distribution_Graph_Pango(Seu_Obj = Seu_Obj(),
                                                         Gene = input$Gene_ID,
                                                         layer = input$Plotting_Data_Layer) +
          theme(plot.margin = margin(40,0,0,40),
                plot.background = element_rect(fill = "white",
                                               colour = "white"))

      }

      if (input$Graph_Type == "meta_data") {

        graph <- PangoICH::Meta_Data_Graph_Pango(Seu_Obj = Seu_Obj(),
                                                 Aim_Var = input$Meta_Data_Var)

      }

      return(graph)
    })

    ##

  }
  return(server)
}

#' @export

Combine_Plots_Shiny_PangoICH <- function(Plots,
                                         Area_List = NULL)
{

  combine_plots_shiny <- shinyApp(ui = .combine_plots_shiny_UI(Plots = Plots),
                                  server = .combine_plots_shiny_Server(Plots = Plots,
                                                                       Area_List = Area_List))

  return(combine_plots_shiny)

}

.combine_plots_shiny_UI <- function(Plots){

  ui <- fluidPage(sidebarLayout(sidebarPanel(width = 3,
                                             selectInput(inputId = "Chosed_Plots",
                                                         label = "Choose Plots",
                                                         choices = c(names(Plots)),
                                                         multiple = TRUE),
                                             textInput(inputId = "Paper_Size",
                                                       label = "Paper_Size",
                                                       value = "",
                                                       updateOn = "blur",
                                                       placeholder = "Width:  Height:  (mm)"),
                                             uiOutput("Direction"),
                                             uiOutput("Plots_Location"),
                                             uiOutput("Save_Plot")
  ),
  mainPanel(width = 9,
            imageOutput("Combined_Plot"))))

  return(ui)

}

.combine_plots_shiny_Server <- function(Plots,
                                        Area_List){

  Tmp_File <- tempfile(fileext = ".png")

  server <- function(input,output,session)
  {

    output$Plots_Location <- renderUI(
      {
        req(input$Chosed_Plots,
            input$Paper_Size,
            input$Direction)

        if(input$Direction == "Horizontal"){

          input_label <- paste(Paper_Width(),Paper_Height(),sep = "X") %>%
            paste("t,l,b,r",sep = ": ")

        } else {

          input_label <- paste(Paper_Height(),Paper_Width(),sep = "X") %>%
            paste("t,l,b,r",sep = ": ")

        }

        n_plots <- length(input$Chosed_Plots)
        text_input_ls <- tagList()
        for (i in 1:n_plots) {

          if (input$Chosed_Plots[i] %in% names(Area_List) == TRUE) {

            text_input_ls <- tagList(text_input_ls,
                                     textInput(inputId = input$Chosed_Plots[i],
                                               label = input$Chosed_Plots[i],
                                               value = "",
                                               updateOn = "blur",
                                               placeholder = paste(as.character(Area_List[[input$Chosed_Plots[i]]]),collapse = ",")))

          } else {

            text_input_ls <- tagList(text_input_ls,
                                     textInput(inputId = input$Chosed_Plots[i],
                                               label = input$Chosed_Plots[i],
                                               value = "",
                                               updateOn = "blur",
                                               placeholder = input_label))


          }
        }
        return(text_input_ls)
      }
    )

    Spacer_Area <- reactive({
      req(input$Direction)

      if (input$Direction == "Horizontal") {

        return(list(plot_spacer = area(1,1,Paper_Width(),Paper_Height())))

      }
      if (input$Direction == "Vertical") {

        return(list(plot_spacer = area(1,1,Paper_Height(),Paper_Width())))

      }

    })

    Plot_Area <- reactive(
      {

        n_plots <- length(input$Chosed_Plots)
        area_ls <- vector("list",length = n_plots)
        names(area_ls) <- input$Chosed_Plots
        for (i in 1:n_plots) {

          area_na <- input$Chosed_Plots[i]

          if (area_na %in% names(Area_List) == TRUE && input[[area_na]] == "") {

            area_ls[[area_na]] <- Area_List[[area_na]]

          } else {

            req(input[[area_na]])
            area_num <- input[[area_na]] %>%
              strsplit(split = ",")
            area_num <- area_num[[1]] %>%
              as.numeric()
            if (length(area_num) == 4) {

              area_plot <- area(area_num[1],area_num[2],area_num[3],area_num[4])
              area_ls[[area_na]] <- area_plot

            }
          }
        }
        return(area_ls)
      })

    Paper_Width <- reactive({

      req(input$Paper_Size)

      paper_size <- strsplit(input$Paper_Size,split = ",")
      paper_size <- paper_size[[1]]

      return(as.integer(paper_size[1]))

    })

    Paper_Height <- reactive({

      req(input$Paper_Size)

      paper_size <- strsplit(input$Paper_Size,split = ",")
      paper_size <- paper_size[[1]]

      return(as.integer(paper_size[2]))

    })

    Combined_Plot <- reactive(
      {

        req(Plot_Area(),Spacer_Area(),input$Direction)

        plots_location <- append(Spacer_Area(),Plot_Area())
        plots_location <- do.call(c,plots_location)

        plots_ls <- append(list(plot_spacer()),Plots[input$Chosed_Plots])
        combined_plot <- wrap_plots(plots_ls,design = plots_location)

        Combined_Plot_with_backgroud <- combined_plot +
          plot_annotation(theme = theme(plot.background = element_rect(colour = "black")))
        Combined_Plot_with_backgroud <- patchworkGrob(Combined_Plot_with_backgroud) %>%
          as.ggplot()

        if (input$Direction == "Vertical") {

          ggsave(Tmp_File,Combined_Plot_with_backgroud,width = Paper_Height(),height = Paper_Width(),units = "mm",dpi = 600)

        }
        if (input$Direction == "Horizontal") {

          ggsave(Tmp_File,Combined_Plot_with_backgroud,width = Paper_Width(),height = Paper_Height(),units = "mm",dpi = 600)

        }

        return(combined_plot)

      })

    output$Combined_Plot <- renderImage({

      req(Combined_Plot())

      list(src = Tmp_File,
           contentType = "image/png",
           width = "50%",
           height = "auto")
    },deleteFile = FALSE)

    output$Direction <- renderUI({

      req(input$Paper_Size)
      selectInput(inputId = "Direction",
                  label = "Direction",
                  choices = c("",
                              "Horizontal",
                              "Vertical"),
                  multiple = FALSE)

    })

    output$Save_Plot <- renderUI({

      actionButton(inputId = "Save_Plot",
                   label = "Save Plot")

    })

    observeEvent(input$Save_Plot, {

      req(Combined_Plot,input$Direction)

      if (input$Direction == "Vertical") {

        ggsave(plot = Combined_Plot(),
               filename = "vertical.png",
               device = "png",
               path = "/home/youngxin/Documents/PangoICH/Graph",
               width = Paper_Height(),
               height = Paper_Width(),
               units = "mm",
               dpi = 600)

      }
      if (input$Direction == "Horizontal") {

        ggsave(plot = Combined_Plot(),
               filename = "horizontal.png",
               device = "png",
               path = "/home/youngxin/Documents/PangoICH/Graph",
               width = Paper_Width(),
               height = Paper_Height(),
               units = "mm",
               dpi = 600)

      }
    })

  }
  return(server)
}

