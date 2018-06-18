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

#
# setwd("~/Transfer/shinydashboard/app/")

# setting colors
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

# loading data ----
metadata.FILE <-
  read.delim(file = "../data/LMO.time.series.metadata.csv", stringsAsFactors = FALSE)
KEGG.tpm <- read.delim(file = "../data/lmo2012_transect2014_redox2014.KEGG-pathway-module.tpm.tsv")
eggNOG.tpm <- read.delim(file = "../data/lmo2012_transect2014_redox2014.eggNOG.tpm.tsv")


envdata <- read.delim(file = "../data/lmo2012_transect2014_redox2014.env-data.tsv")
# envdata.long <- envdata %>% 
#   gather(c(lmo2012, redoxgradient2014, transect2014), key="sample", value="value")

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



metadata.FILE[, 1] <- as.factor(metadata.FILE[, 1])

metadata.FILE.melt <- melt(
  metadata.FILE,
  id.vars = c("SampleID"),
  measure.vars = c(
    "Temperature",
    "DOC",
    "Ammonium",
    "Salinity",
    "Chla",
    "Nitrate",
    "Phosphate",
    "Silicate",
    "TotalN",
    "NP",
    "BacterialAbundance",
    "BacterialProduction",
    "BacterialProduction.1"
  )
)

metadata.FILE.melt[, 1] <-
  as.Date(metadata.FILE.melt[, 1], format = "%y%m%d")



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
      
      menuItem("Settings", tabName = "settings", icon = icon("bar-chart")),
      
      menuItem("Heatmap", tabName = "dashboard", icon = icon("bar-chart")),
      
      menuItem(
        "Environmental data",
        tabName = "contextual",
        icon = icon("database")
      ),
      menuItem(
        "View data table",
        tabName = "description",
        icon = icon("pencil-square-o")
      ),
      
      menuItem(
        "Random forest prediction",
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
              h2("Settings"),
              
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
                    title = "1. Select sample group",
                    width = NULL,
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    
                    # checkboxInput("lmo2012", "LMO station 2012", TRUE),
                    selectInput('lmo_dataset_list', 'Timepoints LMO 2012', lmo2012, multiple=TRUE, selectize=FALSE),
                    
                    # checkboxInput("transect2014", "Baltic Sea transect 2014", FALSE),
                    selectInput('transect_dataset_list', 'Stations Transect 2014', transect2014, multiple=TRUE, selectize=FALSE),
                    
                    # checkboxInput("gradient2014", "Baltic Sea redox gradient 2014", FALSE),
                    selectInput('redox_dataset_list', 'Depths Redox gradient 2014', redoxgradient2014, multiple=TRUE, selectize=FALSE)
                    
                  )
                  
                  
                  
                ),
                
                column(
                  width = 4,
                  
                  box(
                    title = "2. Select annotation data",
                    width = NULL,
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    
                    radioButtons("annotation_data", "Annotation data", c("eggNOG" = "eggNOK", 
                                                                         "KEGG" = "KEGG"))
                    
                    
                  ), 
                  
                  box(
                    title = "3. Select parameters to display",
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
                    
                    title = "4. Select date range",
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
                  
                  box(title = "5. Upload external sample",
                      width = NULL,
                      solidHeader = TRUE,
                      collapsible = TRUE,
                      
                      fileInput("file1", "Choose data File",
                                accept = c(
                                  "text/csv",
                                  "text/comma-separated-values,text/plain",
                                  ".csv")),
                      checkboxInput("somevalue", "include external sample", FALSE)
                  )
                  
                )
                
                
              )),
      
      
      
      
      
      
      
      tabItem(tabName = "dashboard",
              h2("Heatmap"),
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
                    title = "Settings",
                    width = NULL,
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    "Change some settings of the heatmap here.",
                    sliderInput("unifRange", "Columns", min = 1, max = ncol(KEGG.tpm), value = c(1, 15)),
                    numericInput("normCount", "Rows", 20)
                    
                  )
                )
                
                
              )),
      
      
      
      
      
      tabItem(tabName = "contextual",
              h2("Contextual data"),
              fluidRow(
                column(
                  width = 8,
                  
                  box(
                    plotlyOutput("contextual"),
                    width = NULL,
                    title = "Contextual data plot",
                    solidHeader = TRUE,
                    collapsible = TRUE
                  )
                  
                ),
                
                column(
                  width = 4,
                  box(
                    title = "Select data to display",
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
              h2("Data table"),
              
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
              h2("Random forest predictions"),
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
    
    
    fluidRow(
      infoBoxOutput("diversityBox"),
      infoBoxOutput("sampleBox"),
      infoBoxOutput("statusBox")
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
  
  
  
# subsetting data interactively
  
  selectedData <- reactive({
    KEGG.tpm[c(1:input$normCount), c(input$lmo_dataset_list, input$transect_dataset_list, input$redox_dataset_list)]
  })
  
#datatable
  
  {
    output$tbl = DT::renderDT(
      selectedData(), options = list(lengthChange = FALSE,
                                                                              dom = 'Bfrtip', 
                                                                              buttons = I(c('colvis', 'copy', 'csv', 'excel'))),
      extensions = 'Buttons'
    )
  }
  
  
  # plot on action
  
  v <- reactiveValues(doPlot = FALSE)
  
  observeEvent(input$go, {
    # 0 will be coerced to FALSE
    # 1+ will be coerced to TRUE
    v$doPlot <- input$go
  })
  
  output$plot <- renderPlot({
    if (v$doPlot == FALSE) return()
    
    isolate({
      data <- if (input$tabset == "dashboard") {
        runif(input$unifCount, input$unifRange[1], input$unifRange[2])
      } else {
        rnorm(input$normCount, input$normMean, input$normSd)
      }
      
      hist(data)
    })
  })
  
  # import csv file
  output$contents <- renderTable({
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    read.csv(inFile$datapath, header = input$header)
  })
  
  # checkboxes for data selection
  output$value <- renderText({ input$somevalue })
  output$value2 <- renderText({ input$somevalue2 })
  
  
  
  plotcontextualdata <- reactive({
    plotcontextualdata <-
      metadata.FILE.melt[which(
        metadata.FILE.melt[, "SampleID"] >= as.Date(input$dates[1], format = "%Y%m%d") &
          metadata.FILE.melt[, "SampleID"] <= as.Date(input$dates[2], format = "%Y%m%d") &
          metadata.FILE.melt[, "variable"] == input$nutrients
      ), ]
  })
  
  
  output$plot1 <- renderPlotly({
    plot1 <-  ggplot(data = plot1data()) +
      geom_bar(aes(x = date, y = value, fill = mainrole), stat = "identity") +
      theme(axis.text.x = element_text(angle = 90))
    ggplotly(plot1)
  })
  
  
  output$plot2 <- renderPlotly({
    plot2 <- ggplot(data = plot2data()) +
      geom_bar(aes(y = mean, x = mainrole, group = mainrole), stat = "identity") +
      theme(axis.text.x = element_text(angle = 90)) +
      coord_flip()
    ggplotly(plot2)
  })
  
  output$contextual <- renderPlotly({
    contextualplot <- ggplot(na.omit(plotcontextualdata())) +
      geom_line(aes(x = SampleID, y = value)) +
      geom_point(aes(x = SampleID, y = value)) +
      scale_x_date(limits = c(input$dates[1], input$dates[2])) + 
      ylab("Measured values") +
      xlab("Time")
    ggplotly(contextualplot)
  })
  

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
      scale(selectedData())
    )
  })
  
  
}

shinyApp(ui, server)