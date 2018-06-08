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
library(kableExtra)


# setting colors
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

# loading data ----
metadata.FILE <-
  read.delim(file = "../data/LMO.time.series.metadata.csv", stringsAsFactors = FALSE)
KEGG.tpm <- read.delim(file = "../data/lmo2012_transect2014_redox2014.KEGG-pathway-module.tpm.tsv")

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
      menuItem("Heatmap", tabName = "dashboard", icon = icon("bar-chart")),
      
      menuItem(
        "Contextual data",
        tabName = "contextual",
        icon = icon("database")
      ),
      menuItem(
        "View data table",
        tabName = "description",
        icon = icon("pencil-square-o")
      ),
      menuItem(
        "Individual Samples",
        tabName = "samples",
        icon = icon("list")
      ),
     
        
        
        sliderInput("unifRange", "Columns", min = 1, max = ncol(KEGG.tpm), value = c(1, 15)),
        numericInput("normCount", "Rows", 20),
        actionButton("go", "Plot"),
      
      fileInput("file1", "Choose data File",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")
      ),
      
      dateRangeInput(
        "dates",
        label = h4("Date range"),
        start = "2012-01-01",
        end = "2012-12-31",
        min = "2012-01-01",
        max = "2012-12-31"
      )
    )
  ),
  
  # body content ----
  dashboardBody(
    tabItems(
      # First tab content
      tabItem(tabName = "dashboard",
              h2("Heatmap"),
              fluidRow(
                column(
                  width = 12,
                  box(
                    title = "Information",
                    width = NULL,
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    "Select the date range in the side bar on the left side. The graph and boxes with indicator values will be updated automatically. Click and drag to select specific dates or categories."
                  ),
                  box(
                    d3heatmapOutput("heatmap"),
                    title = "Heatmap plot",
                    width = NULL,
                    solidHeader = TRUE,
                    collapsible = TRUE
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
                  ),
                  box(
                    plotlyOutput("plotderivative"),
                    width = NULL,
                    title = "1st derivative",
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
                  "This table displays the data which you selected for on the left side and which is displayed in the heatmap."
                ),
                box(
                  title = "View Data",
                  solidHeader = TRUE,
                  width = NULL,
                  DTOutput('tbl')
                )
              )),
      
      
      tabItem(
        tabName = "samples",
        h2("Features of individual samples"),
        
        column(
          width = 12,
          box(
            title = "List of samples",
            solidHeader = TRUE,
            width = NULL,
            "Indicator features calculated separately for each sample."
          )
        )
      )
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

  #datatable
  
  {
    output$tbl = renderDT(
      KEGG.tpm[c(1:input$normCount), c(1:input$unifRange[2])], options = list(lengthChange = FALSE)
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
  
 
  
  output$heatmap <- renderD3heatmap({
    d3heatmap(
      scale(KEGG.tpm[c(1:input$normCount), c(1:input$unifRange[2])])
    )
  })
  
  
}

shinyApp(ui, server)