## dashboard.app.R 
library(shiny)
library(shinydashboard)
library(ggplot2)
library(reshape)
library(d3heatmap)
require(RColorBrewer)
library(Hmisc)  # for statistics
library(igraph) # for igraph and networks
library(networkD3)
library(vegan)
library(plyr)
library(plotly)


# setting colors
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

# loading data ----
metadata.FILE <- read.delim(file = "LMO.time.series.metadata.csv", stringsAsFactors = FALSE)
lmo.mg.TIGRFAM.mainrole <- read.delim(file = "lmo.mg.TIGRFAM.mainrole.tab", stringsAsFactors = FALSE)[-c(11,19,20),] # excluding "no TIGRFAM", "Unclassified", "Unknown function"
colnames(lmo.mg.TIGRFAM.mainrole) <- gsub("X", "", colnames(lmo.mg.TIGRFAM.mainrole))
lmo.mg.TIGRFAM.subrole <- read.delim(file = "lmo.mg.TIGRFAM.subrole.tab", stringsAsFactors = FALSE)[-c(11,19,20),] # excluding "no TIGRFAM", "Unclassified", "Unknown function"
colnames(lmo.mg.TIGRFAM.subrole) <- gsub("X", "", colnames(lmo.mg.TIGRFAM.subrole))


lmo.mg.TIGRFAM.mainrole <- na.omit(lmo.mg.TIGRFAM.mainrole)

lmo.mg.TIGRFAM.mainrole.long <- melt(data = lmo.mg.TIGRFAM.mainrole, id.vars = c("mainrole"), variable_name = "date")
lmo.mg.TIGRFAM.mainrole.long[, "date"] <- as.Date(lmo.mg.TIGRFAM.mainrole.long[, "date"], format = "%Y%m%d")



lmo.mg.TIGRFAM.mainrole.sh <- na.omit(lmo.mg.TIGRFAM.mainrole)
lmo.mg.TIGRFAM.mainrole <- na.omit(lmo.mg.TIGRFAM.mainrole)
rownames(lmo.mg.TIGRFAM.mainrole.sh) <- (lmo.mg.TIGRFAM.mainrole[,1])
lmo.mg.TIGRFAM.mainrole.sh <- lmo.mg.TIGRFAM.mainrole.sh[,-1]
lmo.mg.TIGRFAM.mainrole.sh.t <- t(lmo.mg.TIGRFAM.mainrole.sh)

shannon <- vegan::diversity(lmo.mg.TIGRFAM.mainrole.sh.t, index = "shannon")
shannon <- as.data.frame(shannon)
shannon$date <- row.names(shannon)

diversity.melt <- melt(shannon, id.vars = "date")
diversity.melt$date <- as.Date(diversity.melt$date, format = "%Y%m%d")
diversity <- t(diversity.melt[,c(1,3)])


rownames(lmo.mg.TIGRFAM.mainrole) <- lmo.mg.TIGRFAM.mainrole[, 1]
lmo.mg.TIGRFAM.mainrole <- lmo.mg.TIGRFAM.mainrole[, -1]
# lmo.mg.TIGRFAM.mainrole <- lmo.mg.TIGRFAM.mainrole[-c(10,18), ]




metadata.FILE[,1] <- as.factor(metadata.FILE[,1])

metadata.FILE.melt <- melt(metadata.FILE, id.vars = c("SampleID"), 
                           measure.vars = c("Temperature", 
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
                                            "BacterialProduction.1"))

metadata.FILE.melt[,1] <- as.Date(metadata.FILE.melt[,1], format = "%y%m%d")



# BEGIN UI #############################################################################################

ui <- dashboardPage(
  skin = "green",
  
  # header content ----  
  dashboardHeader(title = "LMO metagenome indicators", 
                  titleWidth = 300),
  
  # sidebar content ----
  dashboardSidebar(
    width = 300,
    sidebarMenu(
      menuItem("Overview", tabName = "dashboard", icon = icon("bar-chart")),
      menuItem("Indicators", tabName = "indicators", icon = icon("th")),
      menuItem("TIGRFAM", tabName = "TIGRFAM", icon = icon("align-center")),
      menuItem("Contextual data", tabName = "contextual", icon = icon("database")),
      menuItem("Data description", tabName = "description", icon = icon("pencil-square-o")),
      dateRangeInput("dates", label = h4("Date range"), start = "2012-01-01", end = "2012-12-31", min = "2012-01-01", max = "2012-12-31")

            # uiOutput('resetable_input'),
      # actionButton("reset_date", label = "reset date range")
      
      
    )
  ),
  
  # body content ----  
  dashboardBody(
    tabItems(
      
      
      # First tab content
      tabItem(tabName = "dashboard",
              h2("Dashboard overview"),
              fluidRow(
                
                column(width = 6,
                       box("Click and drag to select specific dates or categories.",
                           d3heatmapOutput("heatmap"),  
                           title = "Heatmap plot", width = NULL, solidHeader = TRUE, collapsible = TRUE),
                       box(title = "Information", width = NULL, solidHeader = TRUE, collapsible = TRUE,
                           "Select the date range in the side bar on the left side. The graph and boxes with indicator values will be updated automatically.")
                ),
                
                column(width = 6,
                       box(title = "Network view", width = NULL, solidHeader = TRUE, collapsible = TRUE,
                           "Use the mouse wheel to zoom.",
                           # sliderInput("opacity", label = h4("Network Opacity"), "Decimal:", min = 0, max = 1, value = 0.8, step= 0.1, dragRange = FALSE),
                           forceNetworkOutput("force"))
                       
                )                
              )
              
              
              
              
      ),
      
      
      
      
      
      
      tabItem(tabName = "TIGRFAM", 
              h2("TIGRFAM categories"),
              
              fluidRow(
                column(width = 8,
                       box(plotlyOutput("plot1"), width = NULL,  title = "TIGRFAM bar plot", solidHeader = TRUE, collapsible = TRUE),
                       box(plotlyOutput("plot2"), width = NULL,  title = "TIGRFAM bar plot s2", solidHeader = TRUE, collapsible = TRUE)
                ),
                column(width = 4,
                       box(title = "Protein families", width = NULL, solidHeader = TRUE,
                           "Shows the TIGRFAM protein families.")
                )
                
              )
      ),
      
      
      
      tabItem(tabName = "indicators",
              h2("Indicators"),
              
              column(width = 8,
                     box(title ="Overview", solidHeader = TRUE, width = NULL, collapsible = TRUE,
                         
                         h3("Main objectives"),
                         "The main objective is to provide a description of the current ecosystem status that integrates metagenomic data of microbial populations.",
                         
                         h3("Core indicators"),
                         
                         "Indicators breaking down characteristics of the dataset. Core indicators with relevance to the Baltic Sea are: ", 
                         HTML("
                              <br><br>
                              <ul>
                              <li>Microbial abundance</li>
                              <li>Taxonomic diversity<b>*</b></li>
                              <li>Cyanobacteria abundance<b>*</b></li>
                              </ul><br>
                              <b>*</b> = based on metagenomic data
                              "),
                         
                         h3("Additional indicators"),
                         
                         HTML("
                              <ul>
                              <li>Temperature onset of nodularin gene-expression<b>*</b></li>
                              
                              </ul>
                              ")
                         
                         ))
              
                     ),
      
      
      
      tabItem(tabName = "contextual", 
              h2("Contextual data"),
              fluidRow(
                column(width = 8,
                       
                       box(plotlyOutput("contextual"), width = NULL, title = "Contextual data plot", solidHeader = TRUE, collapsible = TRUE)
                       
                ),
                
                column(width = 4,
                       box(title = "Select data to display", width = NULL, solidHeader = TRUE, selectInput(inputId = "nutrients",
                                                                                                           label = "Variable:",
                                                                                                           choices = c("Nitrate", "Phosphate", "NP", "Temperature", "DOC", "Ammonium", "Salinity", "Chla", "Silicate", "TotalN", "BacterialAbundance", "BacterialProduction", "BacterialProduction.1"),
                                                                                                           selected = "Temperature"))
                )
              )
              
              
      ),
      
      
      
      tabItem(tabName = "description",
              h2("Description of the data"),
              
              column(width = 12, box(title = "Description", solidHeader = TRUE, width = NULL,
                                     "Metagenomic data from 2012 retrieved from sampling at the Linnaeus Microbial Observatory (LMO) station. (...)"))
              
              
      )
      
      
                         ),
    
    
    fluidRow(
      infoBoxOutput("diversityBox"),
      infoBoxOutput("sampleBox"),
      infoBoxOutput("statusBox")
    ),
    
    div(class = "footer", p("LMO indicator dashboard. Version as of 2016-12-16, carlo.berg@scilifelab.se"))

                     )

                     )


# BEGIN SERVER ###################################################
server <- function(input, output) { 

    # Network graph, in progres....
    
  output$force <- renderForceNetwork({
    
    lmo.mg.TIGRFAM.mainrole.sh.t.x <- lmo.mg.TIGRFAM.mainrole.sh.t[, -c(10,18)]
    
    lmo.corrtable <- rcorr(x = as.matrix(lmo.mg.TIGRFAM.mainrole.sh.t.x[,c(2:ncol(lmo.mg.TIGRFAM.mainrole.sh.t.x))]), 
                           type = "spearman") # use "spearman" or "pearson"
    
    # extract correlation tables with r and P values ----
    lmo.corrtable.r <- lmo.corrtable$r
    lmo.corrtable.P <- lmo.corrtable$P
    
    #making edgelists
    el.r <- as.data.frame.table(lmo.corrtable.r)
    el.P <- as.data.frame.table(lmo.corrtable.P)
    el <- cbind(el.r, el.P[, 3])
    colnames(el) <- c("Var1", "Var2", "r", "P")
    
    
    # remove self-loops from edgelist ----
    el <- el[which(el[,1] != el[,2]), ]
    
    # setting cutoffs for r and P
    el <- el[which((el["r"] > 0.3 | el["r"] < -0.2) & el["P"] < 0.001), ]
    
    # assigning red color to negative correlations ----
    el["color"] <- "darkgray"
    el[which((el["r"] < -0.4) & el["P"] < 0.0001), "color"] <- "red"
    el.color <- el["color"]
    
    
    nodes <- as.data.frame(rownames(lmo.corrtable.r))
    colnames(nodes) <- "nodes"
    nodes["color"] <- c("green", "green", "green", "green", "green", "black", "green", "black", "green", "black", "green", "green", "green", "green", "green")
    
    el[, 3] <- as.integer(rep(6, nrow(el)))
    el <- el[, c(1:3)]
    
    g <- graph.edgelist(as.matrix(el[,1:2]), directed = FALSE) 
    
    g_d3 <- igraph_to_networkD3(g, group = membership(cluster_walktrap(g)))
    
    forceNetwork(Links = g_d3$links, Nodes = g_d3$nodes, 
                 Source = 'source', Target = 'target', 
                 NodeID = 'name', Group = 'group', zoom =T)    

  })  
  
  
  
  
  
  plot1data <- reactive({
    plot1data <- lmo.mg.TIGRFAM.mainrole.long[ which(lmo.mg.TIGRFAM.mainrole.long[, "date"] >= as.Date(input$dates[1], format = "%Y%m%d") & 
                                                    lmo.mg.TIGRFAM.mainrole.long[, "date"] <= as.Date(input$dates[2], format = "%Y%m%d")), ] 
    })
  
  plotcontextualdata <- reactive({ 
    plotcontextualdata <- metadata.FILE.melt[ which(metadata.FILE.melt[, "SampleID"] >= as.Date(input$dates[1], format = "%Y%m%d") & 
                                                      metadata.FILE.melt[, "SampleID"] <= as.Date(input$dates[2], format = "%Y%m%d") &
                                                      metadata.FILE.melt[, "variable"] == input$nutrients), ] 
  })
  
diversity.index <- reactive({ 
  # calculate average shannon index
  diversity.index <- mean(diversity.melt[which(diversity.melt[, "date"] >= as.Date(input$dates[1], format = "%Y%m%d") & 
                                                 diversity.melt[, "date"] <= as.Date(input$dates[2], format = "%Y%m%d")), "value"])
  })
  
num.samples <- reactive({ 
  num.samples <- length(diversity.melt[which(diversity.melt[, "date"] >= as.Date(input$dates[1], format = "%Y%m%d") & 
                                               diversity.melt[, "date"] <= as.Date(input$dates[2], format = "%Y%m%d")), "value"])
  })
  
output$plot1 <- renderPlotly({
   plot1 <-  ggplot(data = plot1data()) + 
     geom_bar(aes(x = date, y = value, fill = mainrole), stat = "identity") + 
     theme(axis.text.x = element_text(angle = 90))
   ggplotly(plot1)
  })
  
plot2data <- reactive({
  plot2data <- ddply(lmo.mg.TIGRFAM.mainrole.long, ~ mainrole, summarise, mean = mean(value), sd = sd(value))
  })
  
output$plot2 <- renderPlotly({
        plot2 <- ggplot(data = plot2data()) + 
          geom_bar(aes(y = mean, x = mainrole, group = mainrole), stat = "identity") + 
          theme(axis.text.x = element_text(angle = 90)) + 
          coord_flip()
        ggplotly(plot2)
        })

output$contextual <- renderPlotly({
  contextualplot <- ggplot(plotcontextualdata()) + 
        geom_line(aes(x = SampleID, y = value)) + 
        geom_point(aes(x = SampleID, y = value)) +
        ylab("Measured values") +
        xlab("Time")
      ggplotly(contextualplot)
    })

  plotheatmapdata <- reactive({ 
    plotheatmapdata <- lmo.mg.TIGRFAM.mainrole[, c(1, which(as.Date(colnames(lmo.mg.TIGRFAM.mainrole), format = "%Y%m%d") >= as.Date(input$dates[1], format = "%Y%m%d") & 
                                                              as.Date(colnames(lmo.mg.TIGRFAM.mainrole), format = "%Y%m%d") <= as.Date(input$dates[2], format = "%Y%m%d")
    )) ] 

  })  

  
  output$heatmap <- renderD3heatmap({
    d3heatmap(plotheatmapdata(), colors = cols, scale="column")
  })
  

# InfoBoxes
  output$diversityBox <- renderInfoBox({
    infoBox(
      "Average Diversity Index", round(diversity.index(), digits = 3), icon = icon("pie-chart"),
      color = "navy", fill = TRUE
    )
  })
  
  output$sampleBox <- renderInfoBox({
    infoBox(
      "Number of samples", num.samples(), icon = icon("list"),
      color = "maroon", fill = TRUE
    )
  })
  output$statusBox <- renderInfoBox({
    infoBox(
      "Summary", "80%", icon = icon("thumbs-up", lib = "glyphicon"),
      color = "green", fill = TRUE
    )
  })
  
}

shinyApp(ui, server)