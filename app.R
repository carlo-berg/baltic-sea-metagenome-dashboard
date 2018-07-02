## dashboard.app.R
# https://github.com/carlo-berg/baltic-sea-metagenome-dashboard

library(shiny)
library(shinydashboard)
library(reshape)
library(ggplot2)
library(plotly)
library(d3heatmap)
library(DT)
library(RColorBrewer)
library(tidyverse)
library(magrittr)
library(randomForest)

# setwd("~/Transfer/shinydashboard/app/")

# setting colors
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

# loading data ----

KEGG.tpm <- read.delim(file = "../data/lmo2012_transect2014_redox2014.KEGG-pathway-module.tpm.tsv")
eggNOG.tpm <- read.delim(file = "../data/lmo2012_transect2014_redox2014.eggNOG.tpm.tsv")

envdata <- read.delim(file = "../data/lmo2012_transect2014_redox2014.env-data.tsv")
envdata_t <- envdata %>% 
  t() 

colnames(envdata_t) <- envdata_t[1,]
envdata_t <- envdata_t[-1,]
envdata_t <- as.data.frame(envdata_t) %>% 
  mutate(samples=rownames(envdata_t))

envdata_t[, "Sal"] <- as.numeric(as.character(envdata_t[, "Sal"]))
envdata_t[, "Depth"] <- as.numeric(as.character(envdata_t[, "Depth"]))
envdata_t[, "DOC"] <- as.numeric(as.character(envdata_t[, "DOC"]))
envdata_t[, "O2"] <- as.numeric(as.character(envdata_t[, "O2"]))
envdata_t[, "Temp"] <- as.numeric(as.character(envdata_t[, "Temp"]))
envdata_t[, "NH4"] <- as.numeric(as.character(envdata_t[, "NH4"]))
envdata_t[, "NO3"] <- as.numeric(as.character(envdata_t[, "NO3"]))
envdata_t[, "PO4"] <- as.numeric(as.character(envdata_t[, "PO4"]))
envdata_t[, "Chla"] <- as.numeric(as.character(envdata_t[, "Chla"]))

envdata_t[, "Date"] <- as.Date(envdata_t[, "Date"], format="%d/%m/%y")

# sample groups

lmo2012 <- colnames(envdata[c(2:34)])
transect2014 <- colnames(envdata[c(35:64,73:76)])
redoxgradient2014 <- colnames(envdata[c(65:72, 77, 78)])

# preparing data ----

names <- KEGG.tpm[,1]
KEGG.tpm <- KEGG.tpm[,-1]
KEGG.tpm <- apply(KEGG.tpm, 2, as.numeric)
KEGG.tpm <- as.matrix(KEGG.tpm)
rownames(KEGG.tpm) <- names

# list of KEGG modules or eggNOGS
modules <- rownames(KEGG.tpm)
eggNOGs <- as.vector(eggNOG.tpm[c(1:200), "X"])

# BEGIN UI #############################################################################################

ui <- dashboardPage(
  skin = "green",
  
  # header content ----
  dashboardHeader(title = "Baltic Sea metagenome",
                  titleWidth = 300),
  
  # sidebar content ----
  dashboardSidebar(
    width = 300,
    sidebarMenu(
      # icons from https://getbootstrap.com/docs/3.3/components/
      menuItem("1. Settings", tabName = "settings", icon = icon("cog", lib = "glyphicon")),
      
      menuItem("2. Heatmap", tabName = "dashboard", icon = icon("equalizer", lib = "glyphicon")),
      
      menuItem(
        "3. Environmental data",
        tabName = "contextual",
        icon = icon("bar-chart")
      ),
      menuItem(
        "4. View and export count data",
        tabName = "description",
        icon = icon("save", lib = "glyphicon")
      ),
      
      menuItem(
        "5. Prediction of env. parameters",
        tabName = "predict",
        icon = icon("tree-deciduous", lib = "glyphicon")
      )
      
    )
  ),
  
  # body content ----
  dashboardBody(
    tabItems(
      
      # First tab content

      tabItem(tabName = "settings",
              h2("1. Settings"),
              
              fluidRow(
                column(
                  width = 12,
                  box(title = "Information",
                      width = NULL,
                      solidHeader = TRUE,
                      collapsible = TRUE,
                      "Choose here samples after specific categories or upload an external sample file to compare to the metagenomic data."
                  )
                )),
              
              
              fluidRow(
                column(
                  width = 4,
                  box(
                    title = "A. Select samples",
                    width = NULL,
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    "Three datasets are available of which samples can be chosen here (A). In the following steps (C, D) you can filter the samples by ranges in environmental parameters or date.",
                    tags$br(),tags$br(),
                    selectInput('lmo_dataset_list', 'Timepoints LMO 2012', lmo2012, multiple=TRUE, selectize=FALSE, selected = lmo2012[c(1:2)]),
                    
                    selectInput('transect_dataset_list', 'Stations Transect 2014', transect2014, multiple=TRUE, selectize=FALSE, selected = transect2014[c(1:2)]),
                    
                    selectInput('redox_dataset_list', 'Depths Redox gradient 2014', redoxgradient2014, multiple=TRUE, selectize=FALSE,selected = redoxgradient2014[c(1:2)])
                    
                  )
                  
                  
                  
                ),
                
                column(
                  width = 4,
          
                  box(
                    title = "B. Filter samples by parameter range",
                    "All criteria are applied together, samples with NA values in one parameter will be kept but may be excluded based on the filtering of the other parameters. By default, the sliders display the range in the data to include all samples.",
                    tags$br(),tags$br(),
                    width = NULL,
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    sliderInput("depth", "Depth:", min = 0, max = 439, value = c(floor(min(na.omit(envdata_t[, "Depth"]))), ceiling(max(na.omit(envdata_t[, "Depth"]))))),
                    sliderInput("salinity", "Salinity:", min = 0, max = 40, value = c(floor(min(na.omit(envdata_t[, "Sal"]))), ceiling(max(na.omit(envdata_t[, "Sal"]))))),
                    sliderInput("oxygen", "Oxygen:", min = 0, max = 450, value = c(floor(min(na.omit(envdata_t[, "O2"]))), ceiling(max(na.omit(envdata_t[, "O2"]))))),
                    sliderInput("temp", "Temperature:", min = 0, max = 30, value = c(floor(min(na.omit(envdata_t[, "Temp"]))), ceiling(max(na.omit(envdata_t[, "Temp"]))))),
                    sliderInput("nh4", "Ammonium:", min = 0, max = 20, value = c(floor(min(na.omit(envdata_t[, "NH4"]))), ceiling(max(na.omit(envdata_t[, "NH4"]))))),
                    sliderInput("no3", "Nitrate:", min = 0, max = 15, value = c(floor(min(na.omit(envdata_t[, "NO3"]))), ceiling(max(na.omit(envdata_t[, "NO3"]))))),
                    sliderInput("po4", "Phosphate:", min = 0, max = 15, value = c(floor(min(na.omit(envdata_t[, "PO4"]))), ceiling(max(na.omit(envdata_t[, "PO4"]))))),
                    sliderInput("chla", "Chlorophyll a:", min = 0, max = 10, value = c(floor(min(na.omit(envdata_t[, "Chla"]))), ceiling(max(na.omit(envdata_t[, "Chla"])))))
                    
                    
                  )
                  
                  
                  
                  
                ),
                
                column(
                  width = 4,
                  
                  box(
                    
                    title = "C. Filter samples by date range",
                    width = NULL,
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    
                    dateRangeInput(
                      "dates",
                      label = h4("Date range"),
                      start = "2012-01-01",
                      end = "2014-12-31",
                      min = "2012-01-01",
                      max = "2014-12-31"
                    )
                  ),
                  box(
                    title = "D. Select functional annotation",
                    width = NULL,
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    
                    radioButtons("annotation_data", "TPM-normalized counts of:", c("KEGG metabolic pathway modules" = "KEGG", "eggNOGs" = "eggNOG"))
                    
                    
                  ),
                  
                  box(title = "E. Upload external sample data",
                      width = NULL,
                      solidHeader = TRUE,
                      collapsible = TRUE,
                      
                      "Tab-separated textfile. The first column must be the functional annotation category, every other column is one sample.",
                      tags$br(), tags$br(),
                      fileInput("external_data_file", "Choose data File",
                                accept = c(
                                  "text/tab-separated-values,text/plain",
                                  ".tsv")
                                )
                      
                  )
                  
                )
                
                
              )),
      
      
      
      
      
      
      
      tabItem(tabName = "dashboard",
              h2("2. Heatmap"),
              fluidRow(
                column(
                  width = 9,
                  box(
                    title = "Information",
                    width = NULL,
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    "Click and draw a rectangel to select specific samples in the heatmap. Re-loading may take a while."
                  ),
                  
                  box(
                    d3heatmapOutput("heatmap"),
                    title = "Heatmap plot",
                    width = NULL,
                    solidHeader = TRUE,
                    collapsible = TRUE
                  )
                  
                ),
                
                column(
                  width = 3,
                  box(
                    title = "A. Settings",
                    width = NULL,
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    selectInput('heatmap_modules', 'Select modules:', modules, multiple=TRUE, selectize=FALSE, selected = modules[c(1:20)]),
                    tags$b("Perform clustering of:"),
                    checkboxInput("cluster_samples", "Samples", FALSE),
                    checkboxInput("cluster_modules", "Modules", TRUE),
                    tags$b("External uploaded sample:"),
                    checkboxInput("use_external_data_file", "include", FALSE)
                    
                  )
                )
                
                
              )),
      
      
      
      
      
      tabItem(tabName = "contextual",
              h2("3. Environmental data"),
              fluidRow(
                column(
                  width = 8,
                  
                  box(
                    "Environmental data for the samples selected in the settings tab is displayed here. Note that values may be missing for some of the samples. Chose the environmental parameter from the dropdown menu on the right.",
                    width = NULL,
                    title = "Information",
                    solidHeader = TRUE,
                    collapsible = TRUE
                  ), 
                  box(
                    plotlyOutput("contextual"),
                    "Below the same plot with flipped coordinates (suited for the redox gradient samples).",
                    plotlyOutput("contextual_redox"),
                    width = NULL,
                    title = "Environmental data plot",
                    solidHeader = TRUE,
                    collapsible = TRUE
                  )
                  
                ),
                
                column(
                  width = 4,
                  box(
                    title = "A. Select environmental parameter",
                    width = NULL,
                    solidHeader = TRUE,
                    selectInput(
                      inputId = "nutrients",
                      label = "Variable:",
                      choices = c(
                        "Temp",
                        "DOC",
                        "NH4",
                        "O2",
                        "Sal",
                        "Chla"
                      ),
                      selected = "Temperature"
                    )
                  )
                )
              )),
      
      
      
      tabItem(tabName = "description",
              h2("4. Data table"),
              fluidRow( 
              column(
                width = 12,
                box(
                  title = "Description",
                  solidHeader = TRUE,
                  width = NULL,
                  "This table displays the data which you selected on the left side and which is displayed in the heatmap."
                ),
                box(
                  title = "View Data",
                  solidHeader = TRUE,
                  width = NULL,
                  DTOutput('tbl')
                )
              )
              )),
      
      
      tabItem(tabName = "predict",
              h2("5. Prediction of environmental parameters from metagenomic data"),
              fluidRow(
                column(
                  width = 12,
                  box(title = "Information",
                      width = NULL,
                      solidHeader = TRUE,
                      collapsible = TRUE,
                      "If you uploaded an external data file, environmental parameters will be predicted using a random forest model that is trained with the Baltic Sea metagenomic data. Select the parameters that you want values to be predicted for from the dropdown menu on the right." 
                  )
                )),
              
              fluidRow(
                column(
                  width = 8,
                  
                  box(
                    title = "Table of predicted values",
                    solidHeader = TRUE,
                    width = NULL,
                    DTOutput('tbl_pred')
                  )
                  
                ),
                
                column(
                  width = 4,
                  box(
                    title = "Select environmental parameter",
                    width = NULL,
                    solidHeader = TRUE,
                    selectInput(
                      inputId = "env_param",
                      label = "Environmental parameter:",
                      choices = c(
                        "Temp",
                        "NH4",
                        "Sal",
                        "DOC",
                        "Chla"
                      ),
                      selected = "Temp"
                    )
                  )
                )
              ))
      
      
      
      
      
      
    ),
    

    
    
    
    
    
    
    
    fluidRow(
      column(
        width = 12,
        box(width = NULL,
            solidHeader = FALSE,
            collapsible = FALSE,
            div(
              class = "footer",
              HTML("Baltic Sea metagenome dashboard. Carlo Berg (<a href=http://edu.cberg.de>http://edu.cberg.de</a> | <a href=mailto:cb@edu.cberg.de>cb@edu.cberg.de</a>), Anders F. Andersson and the <a href=https://blueprint-project.org/>BONUS BLUEPRINT</a> project.<br>The BONUS BLUEPRINT project has received funding from BONUS (Art 185), funded jointly by the EU and the national funding institutions of Denmark, Sweden, Germany, Finland, and Estonia.")
            ),
            img(src='img/eu.png', align = "left"),
            "  ",
            img(src='img/bonus.png', align = "left"),
            "  ",
            img(src='img/scilifelab.png', align = "left") 
        )
      ))
    
    
    
  )
  
)


# BEGIN SERVER ###################################################
server <- function(input, output) {
  
# clustering yes/no
  
  output$cluster_samples <- reactive({
  
    input$cluster_samples  
  
  })
    
 
  
  # subsetting data interactively
  
  filtered_samples <-reactive({

  filtered_samples <- envdata_t %>% 
    filter( is.na(Sal) | (Sal >= input$salinity[1] & Sal <= input$salinity[2]))  %>% 
    filter( is.na(O2) | (O2 >= input$oxygen[1] & O2 <= input$oxygen[2]))  %>%
    filter( is.na(Depth) | (Depth >= input$depth[1] & Depth <= input$depth[2]))  %>%
    filter( is.na(Temp) | (Temp >= input$temp[1] & Temp <= input$temp[2]))  %>%
    filter( is.na(NH4) | (NH4 >= input$nh4[1] & NH4 <= input$nh4[2]))  %>%
    filter( is.na(NO3) | (NO3 >= input$no3[1] & NO3 <= input$no3[2]))  %>%
    filter( is.na(PO4) | (PO4 >= input$po4[1] & PO4 <= input$po4[2]))  %>%
    filter( is.na(Chla) | (Chla >= input$chla[1] & Chla <= input$chla[2]))  %>%
    filter(Date >= input$dates[1] & Date <= input$dates[2]) %>% 
    select(samples)
  
  }) 
  
  
  
  
  
  
  selectedData <- reactive({
    
    if (input$annotation_data == "KEGG") {
      selectedData <- KEGG.tpm[c(input$heatmap_modules), intersect(filtered_samples()[,1], c(input$lmo_dataset_list, input$transect_dataset_list, input$redox_dataset_list))]
    
      if (input$use_external_data_file == TRUE) {
        ext_data <- read.delim(input$external_data_file$datapath#,
                         # header = TRUE,
                         # sep = "\t"
                         )
        
        print(input$external_data_file)
        row.names(ext_data) <- ext_data[, 1]
        ext_data <- as.matrix(ext_data[-1])
        selectedData <- merge(selectedData, ext_data, by="row.names", all.x = TRUE)
        row.names(selectedData) <- selectedData[, 1]
        selectedData <- selectedData %>% select(-"Row.names")
      } else {
        selectedData <- KEGG.tpm[c(input$heatmap_modules), intersect(filtered_samples()[,1], c(input$lmo_dataset_list, input$transect_dataset_list, input$redox_dataset_list))]
      }
      
      } else  
      
      if (input$annotation_data == "eggNOG") {
        selectedData <- eggNOG.tpm[c(1:input$normCount), c(input$lmo_dataset_list, input$transect_dataset_list, input$redox_dataset_list)]
      }
   
    
  })
  
  
  
  # module list to select
  output$module_list <- reactive({
    
    if (input$annotation_data == "KEGG") {
      module_list <- modules
    } else  
      
      if (input$annotation_data == "eggNOG") {
        module_list <- eggNOGs
      }
    
    
  })
  
  
  #datatable
  
  {
    output$tbl = DT::renderDT(
      selectedData(), 
      extensions = c('Buttons', 'FixedColumns','Scroller'),
      options = list(scrollX = TRUE,
                     deferRender = TRUE,
                     scrollY = 500,
                     scroller = TRUE,
                     fixedColumns = list(leftColumns = 1),
                     lengthChange = FALSE,
                     dom = 'Bfrtip',
                     buttons = I(c('colvis', 'copy', 'csv', 'excel')))
    )
  }
  
  
 # environmental data plots
  output$contextual <- renderPlotly({
    
    envdata_t_plot <- envdata_t %>% 
      filter(samples %in% c(input$lmo_dataset_list, input$transect_dataset_list, input$redox_dataset_list))
    
    contextualplot <- ggplot(envdata_t_plot) +
      geom_line(aes(x = samples, y = envdata_t_plot[, input$nutrients])) +
      geom_point(aes(x = samples, y = envdata_t_plot[, input$nutrients])) +
      # scale_x_date(limits = c(input$dates[1], input$dates[2])) + 
      ylab("Concentration or level") +
      xlab("Sample") +
      ggtitle(input$nutrients) 
    ggplotly(contextualplot)
  })
  
  output$contextual_redox <- renderPlotly({
    
    envdata_t_plot <- envdata_t %>% 
      filter(samples %in% c(input$lmo_dataset_list, input$transect_dataset_list, input$redox_dataset_list))
    
    contextual_redoxplot <- ggplot(envdata_t_plot) +
      geom_line(aes(x = Depth, y = envdata_t_plot[, input$nutrients])) +
      geom_point(aes(x = Depth, y = envdata_t_plot[, input$nutrients])) +
      # scale_x_date(limits = c(input$dates[1], input$dates[2])) + 
      ylab("Concentration or level") +
      xlab("Depth") +
      coord_flip() +
      scale_x_reverse() + 
      ggtitle(input$nutrients) 
    ggplotly(contextual_redoxplot)
  })
  
  
  
  # reading rf object and corresponding feature list
  # correlation plot      
  pred_corr <- reactive({
    
    
    if (input$env_param == "Temp") { rf = readRDS("../data/predict/rf.kegg.temperature.rds") }
    if (input$env_param == "NH4") { rf = readRDS("../data/predict/rf.kegg.nh4.rds") }
    if (input$env_param == "Sal") { rf = readRDS("../data/predict/rf.kegg.salinity.rds") }
    if (input$env_param == "DOC") { rf = readRDS("../data/predict/rf.kegg.doc.rds") }
    if (input$env_param == "Chla") { rf = readRDS("../data/predict/rf.kegg.chla.rds") }
    
    features = readRDS("../data/predict/feature-list.kegg.rds")
    
    # reading count data
    
    # tab <- read.delim("../data/predict/Transect2014_EggNOG.tpm.annotated.tsv")
    tab <- read.delim(input$external_data_file$datapath)
    id <- as.character(tab[,1])
    counts = tab[,2:ncol(tab)]
    
    # sorting count data by the feature list
    ix = match(features, id)
    id = id[ix]
    counts = counts[ix,]
    
    # do the predictions
    pred = c()
    #predict(rf, t(counts[,1]), type="response")
    for (i in 1:ncol(counts)) {
      pred[i] = round(predict(rf, t(counts[,i]), type="response"), digits = 3)
    }
    # pred
    temp_df <- data.frame(colnames(tab[-1]), pred)
    colnames(temp_df) <- c("sample", "predicted_value")
    
    pred_corr <- temp_df
    
  })  
  
  output$tbl_pred = DT::renderDT(
    pred_corr(), 
    extensions = c('Buttons', 'FixedColumns','Scroller'),
    options = list(scrollX = TRUE,
                   deferRender = TRUE,
                   scrollY = 500,
                   scroller = TRUE,
                   fixedColumns = list(leftColumns = 1),
                   lengthChange = FALSE,
                   dom = 'Bfrtip',
                   buttons = I(c('copy', 'csv', 'excel')))
  )
  
  output$heatmap <- renderD3heatmap({
    d3heatmap(
      scale(selectedData()), Colv=input$cluster_samples, Rowv = input$cluster_modules
    )
  })
  
  


  
}

shinyApp(ui, server)