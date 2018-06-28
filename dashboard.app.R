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
envdata_t[, "O2"] <- as.numeric(as.character(envdata_t[, "O2"]))
envdata_t[, "Temp"] <- as.numeric(as.character(envdata_t[, "Temp"]))




# sample groups

lmo2012 <- colnames(envdata[c(2:34)])
transect2014 <- colnames(envdata[c(35:78)])
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
      
      menuItem("1. Settings", tabName = "settings", icon = icon("bar-chart")),
      
      menuItem("2. Heatmap", tabName = "dashboard", icon = icon("bar-chart")),
      
      menuItem(
        "3. Environmental data",
        tabName = "contextual",
        icon = icon("database")
      ),
      menuItem(
        "4. View & export data",
        tabName = "description",
        icon = icon("pencil-square-o")
      ),
      
      menuItem(
        "5. Random forest prediction",
        tabName = "predict",
        icon = icon("pencil-square-o")
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
                    selectInput('lmo_dataset_list', 'Timepoints LMO 2012', lmo2012, multiple=TRUE, selectize=FALSE, selected = c("X120314", "X120322")),
                    
                    selectInput('transect_dataset_list', 'Stations Transect 2014', transect2014, multiple=TRUE, selectize=FALSE),
                    
                    selectInput('redox_dataset_list', 'Depths Redox gradient 2014', redoxgradient2014, multiple=TRUE, selectize=FALSE)
                    
                  )
                  
                  
                  
                ),
                
                column(
                  width = 4,
                  
                  box(
                    title = "B. Select functional annotation",
                    width = NULL,
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    
                    radioButtons("annotation_data", "TPM-normalized counts of:", c("KEGG metabolic pathway modules" = "KEGG", "eggNOGs" = "eggNOG"))
                    
                    
                  ), 
                  
                  box(
                    title = "C. Filter samples by parameter range",
                    "All criteria are applied together, samples with NA values in one parameter will be kept but may be excluded based on the filtering of the other parameters.",
                    tags$br(),tags$br(),
                    width = NULL,
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    sliderInput("depth", "Depth:", min = 0, max = 439, value = c(0, 70)),
                    sliderInput("salinity", "Salinity:", min = 0, max = 35, value = c(1, 15)),
                    sliderInput("oxygen", "Oxygen:", min = 0, max = 350, value = c(200, 350))
                    
                    
                  )
                  
                  
                  
                  
                ),
                
                column(
                  width = 4,
                  
                  box(
                    
                    title = "D. Filter samples by date range",
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
                  
                  box(title = "E. Upload external sample data",
                      width = NULL,
                      solidHeader = TRUE,
                      collapsible = TRUE,
                      
                      "Tab-separated textfile. The first column must be the functional annotation category, every other column is one sample.",
                      tags$br(), tags$br(),
                      fileInput("external_data_file", "Choose data File",
                                accept = c(
                                  "text/csv",
                                  "text/comma-separated-values,text/plain",
                                  ".csv")),
                      checkboxInput("include_external", "include external sample", FALSE)
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
                    checkboxInput("cluster_modules", "Modules", TRUE)
                    
                  )
                )
                
                
              )),
      
      
      
      
      
      tabItem(tabName = "contextual",
              h2("3. Environmental data"),
              fluidRow(
                column(
                  width = 8,
                  
                  box(
                    "Space for the environmental data.",
                    width = NULL,
                    title = "Data plot",
                    solidHeader = TRUE,
                    collapsible = TRUE
                  )
                  
                ),
                
                column(
                  width = 4,
                  box(
                    title = "A. Select data to display",
                    width = NULL,
                    solidHeader = TRUE,
                    selectInput(
                      inputId = "nutrients",
                      label = "Variable:",
                      choices = c(
                        "Nitrate",
                        "Phosphate",
                        "NP",
                        "Temperature",
                        "DOC",
                        "Ammonium",
                        "Salinity",
                        "Chla",
                        "Silicate",
                        "TotalN",
                        "BacterialAbundance",
                        "BacterialProduction",
                        "BacterialProduction.1"
                      ),
                      selected = "Temperature"
                    )
                  )
                )
              )),
      
      
      
      tabItem(tabName = "description",
              h2("4. Data table"),
              
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
              )),
      
      
      tabItem(tabName = "predict",
              h2("5. Random forest prediction"),
              fluidRow(
                column(
                  width = 8,
                  
                  box(
                    plotlyOutput("pred_corr"),
                    width = NULL,
                    title = "Correlation plot measured vs. predicted",
                    solidHeader = TRUE,
                    collapsible = TRUE
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
                        "DOC"
                      ),
                      selected = "Temp"
                    )
                  )
                )
              ))
      
      
      
      
      
      
    ),
    

    div(
      class = "footer",
      HTML(
        "Baltic Sea metagenome dashboard. <a href=mailto:carlo.berg@scilifelab.se>carlo.berg@scilifelab.se</a>"
      )
    )
    
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
    select(samples)
  
  
  }) 
  
  
  
  
  selectedData <- reactive({
    
    if (input$annotation_data == "KEGG") {
      selectedData <- KEGG.tpm[c(input$heatmap_modules), intersect(filtered_samples()[,1], c(input$lmo_dataset_list, input$transect_dataset_list, input$redox_dataset_list))]
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
  
  

  
  
  # reading rf object and corresponding feature list
  # correlation plot      
  output$pred_corr <- renderPlotly({
    
    
    if (input$env_param == "Temp") { rf = readRDS("../data/predict/rf.cog.temperature.rds") }
    if (input$env_param == "NH4") { rf = readRDS("../data/predict/rf.cog.nh4.rds") }
    if (input$env_param == "Sal") { rf = readRDS("../data/predict/rf.cog.salinity.rds") }
    if (input$env_param == "DOC") { rf = readRDS("../data/predict/rf.cog.doc.rds") }
    
    features = readRDS("../data/predict/feature-list.cog.rds")
    
    # reading count data
    
    tab <- read.delim("../data/predict/Transect2014_EggNOG.tpm.annotated.tsv")
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
      pred[i] = predict(rf, t(counts[,i]), type="response")
    }
    pred
    
    # plotting correlation
    envdata <- read.delim(file = "../data/lmo2012_transect2014_redox2014.env-data.tsv")
    samples <- as.character(c("X", colnames(tab[,2:ncol(tab)])))
    
    env <- envdata[, samples] %>% 
      filter(X == input$env_param) %>% 
      select(-X)   
    
    measured <- as.numeric(t(env)[,1])
    
    pred_corr_plot <- ggplot() + geom_point(aes(x=measured, y=pred)) 
    ggplotly(pred_corr_plot)
    
  })  
  
  
  output$heatmap <- renderD3heatmap({
    d3heatmap(
      scale(selectedData()), Colv=input$cluster_samples, Rowv = input$cluster_modules
    )
  })
  
  
}

shinyApp(ui, server)