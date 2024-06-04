list.of.packages <- c("shiny", "shinydashboard", "Seurat", "ggplot2", "plotly")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(shiny)
library(shinydashboard)
library(Seurat)
library(ggplot2)
library(plotly)
data <- readRDS("new_44_rep1.rds")
data@active.assay = "RNA"
data = NormalizeData(data)

ui = fluidPage(
  tabsetPanel(type = "pills",
    tabPanel("Introduction", fluid = TRUE,
        
               mainPanel(fluidPage(
                 uiOutput("img"),
                 fluidPage(tags$div(tags$h2("Beta Version - Description Still Needs to Be Determined")))
                 
               ))
             
    ),
    tabPanel("Metdata", fluid = TRUE,  sidebarLayout(
      sidebarPanel(downloadButton('mety',"Download Data")),
      fluidRow(column(7,shinycssloaders::withSpinner(dataTableOutput('meta'), size = 2, type = 8)), style = "max-height: 80vh; overflow-y: auto;" ))),
    tabPanel("Expression", fluid = TRUE, sidebarLayout(
      sidebarPanel(downloadButton('dept',"Download Data")),
      fluidPage(fluidRow(column(7,shinycssloaders::withSpinner(dataTableOutput('tab'), size = 2, type = 8 )), style = "max-height: 80vh; overflow-y: auto;" )))),
    
    tabPanel("Differential Expression", fluid = TRUE, sidebarLayout(
      sidebarPanel(selectInput('celltype', 'Cell Type 1',choices = c("EPI-like","Mesoderm-1", "Mesoderm-2" ,"ExE-like" , "Ectoderm", "Endoderm","hPGCLC"), selected = "EPI-like"),
                   numericInput('logfc', 'Log Fold Change',0.5,min = 0, max = 50), numericInput('pct', 'Minimum Expression',0.1,min = 0, max = 50),         
                   downloadButton('down23',"Download Data")),
      fluidPage(fluidRow(column(7,shinycssloaders::withSpinner(dataTableOutput('diffyy'), size = 2, type = 8)), style = "max-height: 80vh; overflow-y: auto;")))),
    
    tabPanel("Principal Component Analysis", fluid = TRUE,sidebarLayout(
      sidebarPanel(h4("Settings"),width = 3 ,
                   numericInput('Dim_X', 'Dimension X',1,min = 0, max = 50), 
                   numericInput('Dim_Y', 'Dimension Y',2,min = 0, max = 50),
                   selectInput('label.color', 'Label Color',choices = c("red", "black", "blue")),
                   numericInput("label.size", "Label Size", 4, min = 0, max = 10),
                   numericInput("pt.size", "Point Size", 0.8, min = 0, max = 10),
                   actionButton("RUN", "Make Changes"),
                   downloadButton('down',"Download")
      ),
      mainPanel(tabsetPanel(id = "tabsetPanelID",type = "tabs",tabPanel("2D Plot", plotOutput('pcaplot')),tabPanel("Interactive", plotlyOutput('pcp')))))),
    
    
    tabPanel("UMAP", fluid = TRUE,sidebarLayout(
      sidebarPanel(h4("Settings"),width = 3 ,
                   selectInput('label.coloro', 'Label Color',choices = c("red", "black", "blue")),
                   numericInput("label.sizeo", "Label Size", 4, min = 0, max = 10),
                   numericInput("pt.sizeo", "Point Size", 0.8, min = 0, max = 10),
                   actionButton("RUNo", "Make Changes"),
                   downloadButton('deptyuo',"Download")
      ),
      mainPanel(tabsetPanel(id = "umaptabid",type = "tabs",tabPanel("2D Plot", plotOutput('umapploy')),tabPanel("Interactive", plotlyOutput('ump')))))),
    
    tabPanel("Diffusion Map", fluid = TRUE,sidebarLayout(
      sidebarPanel(h4("Settings"),width = 3 ,numericInput('Dim_Xd', 'Dimension X',1,min = 0, max = 50), 
                   numericInput('Dim_Yd', 'Dimension Y',2,min = 0, max = 50),
                   selectInput('label.colorod', 'Label Color',choices = c("red", "black", "blue")),
                   numericInput("label.sizeod", "Label Size", 4, min = 0, max = 10),
                   numericInput("pt.sizeod", "Point Size", 0.8, min = 0, max = 10),
                   actionButton("RUNod", "Make Changes"),
                   downloadButton('deptyuos',"Download")
      ),
      mainPanel(tabsetPanel(id = "diffusi",type = "tabs",tabPanel("2D Plot", plotOutput('diff')),tabPanel("Interactive", plotlyOutput('dify')))))),
    
    tabPanel("Feature Plotting", fluid = TRUE,sidebarLayout(
      sidebarPanel(selectizeInput("foo1", label = "Gene 2",choices = NULL ), selectizeInput("foo2", label = "Gene 2",choices = NULL ),
                   numericInput('Dim_Xa', 'Dimension X',1,min = 0, max = 50), 
                   numericInput('Dim_Ya', 'Dimension Y',2,min = 0, max = 50),
                   selectizeInput('col1', 'Color 1',choices = c("grey", "red", "blue", "yellow")),
                   selectizeInput('col2', 'Color 2',choices = c("red", "blue", "yellow","grey")),
                   selectInput('label.colora', 'Label Color',choices = c("red", "black", "blue")),
                   numericInput("pt.sizea", "Point Size", 0.8, min = 0, max = 10),
                   actionButton("rfeature", "Make Changes"),
                   downloadButton('down_10',"Download")
      ),
      
      mainPanel(plotOutput("featuooreplot"),plotOutput("featuooreplot1"))))))



server <- function(input, output, session) {
  
  output$img <- renderUI({
    tags$img(src = "https://medicine.wustl.edu/wp-content/uploads/medicine-logo-sample.gif", height="80%", width="90%", align="down-left")
  })
  
  
  
  
  
  
  ## Metadata
  
  output$meta = renderDataTable({met})
  met = data@meta.data
  
  output$mety <- downloadHandler(
    filename = function(){"thename.csv"}, 
    content = function(fname){
      write.csv(met, fname)
    })
  
  ## Average Expression
  
  output$tab = renderDataTable({
    fg = AverageExpression(data, assay = "RNA", group.by = "gast.id")
    fg = fg$RNA
    fg = cbind(rownames(fg), fg)
    fg
  })
  
  output$dept <- downloadHandler(
    filename = function(){"thename.csv"}, 
    content = function(fname){
      write.csv(fg, fname)
    }
  )
  
  
  
  # Differential Expression
  
  cds.markers = reactive({
    cds.markers = FindMarkers(data, ident.1 = input$celltype,assay = "RNA", group.by = "gast.id", logfc.threshold = input$logfc, min.pct = input$pct)
    cds.markers = cbind(rownames(cds.markers), cds.markers)
    cds.markers
  })
  output$diffyy = renderDataTable(cds.markers())
  
  output$down23 <- downloadHandler(
    filename = function(){"thename.csv"}, 
    content = function(fname){
      write.csv(cds.markers(), fname)
    }
  )
  
  
  
  
  # PCA
  
  Dim_X = isolate(reactive({input$Dim_Y}))
  Dim_Y = isolate(reactive({input$Dim_X}))
  
  output$pcaplot<-renderPlot({
    input$RUN
    isolate(DimPlot(data, reduction = "pca",pt.size=input$pt.size, dims = c(Dim_X(),Dim_Y()), label = TRUE, label.size = input$label.size, label.color = input$label.color, repel = T, group.by = "gast.id"))
  })
  
  p <- reactive({
    input$RUN
    isolate(DimPlot(data, reduction = "pca",pt.size=input$pt.size, dims = c(Dim_X(),Dim_Y()), label = TRUE, label.size = input$label.size, label.color = input$label.color, repel = T, group.by = "gast.id"))
  })
  
  dt = reactive(DimPlot(data, reduction = "pca",pt.size=input$pt.size, dims = c(Dim_X(),Dim_Y()), label = TRUE, label.size = input$label.size, label.color = input$label.color, repel = T,  group.by = "gast.id"))
  dt = reactive(dt$data)
  
  output$pcp <- renderPlotly({
    dt = DimPlot(data, reduction = "pca",pt.size=input$pt.size, dims = c(Dim_X(),Dim_Y()), label = TRUE, label.size = input$label.size, label.color = input$label.color, repel = T, group.by = "gast.id")
    dt = dt$data
    plot_ly(dt, x = dt[,1], y = dt[,2], color = ~gast.id)})
  
  output$down <- downloadHandler(
    file = "save.pdf" , # variable with filename
    content = function(file) {
      ggsave(p(), filename = file)
    })
  
  
  # UMAP
  
  output$umapploy<-renderPlot({
    input$RUNo
    isolate(DimPlot(data, reduction = "umap",pt.size=input$pt.sizeo, cells = NULL, cols = NULL, label = TRUE, label.size = input$label.sizeo, label.color = input$label.coloro, repel = T, group.by = "gast.id"))
  })
  
  d <- reactive({
    input$RUNo
    isolate(
      DimPlot(data, reduction = "umap",pt.size=input$pt.sizeo, cells = NULL, cols = NULL, label = TRUE, label.size = input$label.sizeo, label.color = input$label.coloro, repel = T, group.by = "gast.id"))
  })
  
  
  output$ump <- renderPlotly({
    dt = DimPlot(data, reduction = "umap",pt.size=input$pt.size, label = TRUE, label.size = input$label.size, label.color = input$label.color, repel = T, group.by = "gast.id")
    dt = dt$data
    plot_ly(dt, x = dt[,1], y = dt[,2], color = ~gast.id)})
  
  output$deptyuo <- downloadHandler(
    file = "save.pdf" , # variable with filename
    content = function(file) {
      ggsave(d(), filename = file)
    })
  
  # Diffusion
  
  Dim_Xd = isolate(reactive({input$Dim_Yd}))
  Dim_Yd = isolate(reactive({input$Dim_Xd}))
  
  output$diff<-renderPlot({
    input$RUNod
    isolate(DimPlot(data, reduction = "diff",dims = c(Dim_Xd(),Dim_Yd()),pt.size=input$pt.sizeod, cells = NULL, cols = NULL, label = TRUE, label.size = input$label.sizeod, label.color = input$label.colorod, repel = T, group.by = "gast.id"))
  })
  
  dp <- reactive({
    input$RUNod
    isolate(
      DimPlot(data, reduction = "diff",dims = c(Dim_Xd(),Dim_Yd()),pt.size=input$pt.sizeod, cells = NULL, cols = NULL, label = TRUE, label.size = input$label.sizeod, label.color = input$label.colorod, repel = T, group.by = "gast.id"))
  })
  
  
  output$dify <- renderPlotly({
    dt.2 = DimPlot(data, reduction = "diff",dims = c(Dim_Xd(),Dim_Yd()),pt.size=input$pt.sizeod, label = TRUE, label.size = input$label.sizeod, label.color = input$label.colorod, repel = T, group.by = "gast.id")
    dt.2 = dt.2$data
    plot_ly(dt.2, x = dt.2[,1], y = dt.2[,2], color = ~gast.id)})
  
  output$deptyuos <- downloadHandler(
    file = "save.pdf" , # variable with filename
    content = function(file) {
      ggsave(dp(), filename = file)
    })
  
  # Feature Plot
  
  
  observeEvent(data, {
    updateSelectizeInput(session, 'foo1', choices = rownames(data), server = TRUE)
  })
  
  observeEvent(data, {
    updateSelectizeInput(session, 'foo2', choices = rownames(data), server = TRUE)
  })
  
  D_X1 = isolate(reactive({input$Dim_Xa}))
  D_Y1 = isolate(reactive({input$Dim_Ya}))
  
  
  output$featuooreplot<-renderPlot({
    input$rfeature
    isolate(
      
      FeaturePlot(object = data,  dims = c(D_X1(),D_Y1()), cols = c(input$col1, input$col2), features = input$foo1,
                  split.by = NULL, reduction = "umap",label.color = input$label.colora, pt.size = input$pt.sizea)
      
    )
  })
  
  
  output$featuooreplot1<-renderPlot({
    input$rfeature
    isolate(
      
      FeaturePlot(object = data,  dims = c(D_X1(),D_Y1()), cols = c(input$col1, input$col2), features = input$foo2,
                  split.by = NULL, reduction = "umap",label.color = input$label.colora, pt.size = input$pt.sizea) 
      
      
    )
  })
  
  
  ddp <- reactive({
    input$rfeature
    isolate(FeaturePlot(object = data,dims = c(D_X1(),D_Y1()), cols = c(input$col1, input$col2), features = input$foo2,
                        split.by = NULL, reduction = "umap",label.color = input$label.colora, pt.size = input$pt.sizea)) +
      
      FeaturePlot(object = data,  dims = c(D_X1(),D_Y1()), cols = c(input$col1, input$col2), features = input$foo1,
                  split.by = NULL, reduction = "umap",label.color = input$label.colora, pt.size = input$pt.sizea)
    
  })
  
  output$down_10 <- downloadHandler(
    file = "save.pdf" , # variable with filename
    content = function(file) {
      ggsave(ddp(), filename = file)
    })
}


shinyApp(ui = ui, server = server)


