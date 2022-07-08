# NOTE: Crashes every other launch because of retriculate
# Libraries
library(ggplot2)
library(htmlwidgets)
library(shiny)
library(shinyjs)
# GEOquery
# install.packages("BiocManager")
# BiocManager::install("GEOquery")
library(GEOquery)

# Set up reticulate (py3)
PYTHON_DEPENDENCIES = c('pip', 'matplotlib', 'pandas', 'scipy', 'sklearn')
virtualenv_dir = Sys.getenv('VIRTUALENV_NAME')
python_path = Sys.getenv('PYTHON_PATH')

# Create virtual env and install dependencies
reticulate::virtualenv_create(envname=virtualenv_dir, python=python_path)
reticulate::virtualenv_install(virtualenv_dir, packages=PYTHON_DEPENDENCIES, ignore_installed=FALSE)
reticulate::use_virtualenv(virtualenv_dir, required=T)

# Import python functions
reticulate::source_python('./maniNetCluster/pyManifold.py')

# # Testing
# gse = getGEO("GSE21653", GSEMatrix = TRUE)[[1]]
# samp = phenoData(gse)$geo_accession[c(1:10)]
# data = NULL
# for (s in samp) {
#   s_gse = getGEO(s, GSEMatrix = TRUE)
#   if (is.null(data)) {
#     data = as.data.frame(s_gse@dataTable@table[,2])
#     rownames(data) = s_gse@dataTable@table[,1]
#   }
#   else
#     data = cbind(data, s_gse@dataTable@table[,2])
# }
# colnames(data) = samp
# t_gse1 = gse
# t_gse1_data = data[c(1:100), ]
# 
# gse = getGEO("GSE21654", GSEMatrix = TRUE)[[1]]
# samp = phenoData(gse)$geo_accession[c(1:10)]
# data = NULL
# for (s in samp) {
#   s_gse = getGEO(s, GSEMatrix = TRUE)
#   if (is.null(data)) {
#     data = as.data.frame(s_gse@dataTable@table[,2])
#     rownames(data) = s_gse@dataTable@table[,1]
#   }
#   else
#     data = cbind(data, s_gse@dataTable@table[,2])
# }
# colnames(data) = samp
# t_gse2 = gse
# t_gse2_data = data[c(1:100), ]

# # Sample
# mat1 = t_gse1_data
# mat2 = t_gse2_data
#
# aligned = ManiNetCluster(
#   t(mat1), t(mat2),
#   nameX='s1',nameY='s2',
#   corr=Correspondence(matrix=cor(mat1, mat2)),
#   d=as.integer(2),
#   method='nonlinear manifold aln',
#   k_NN=as.integer(3),
#   k_medoids=as.integer(6)
# )

# UI
ui <- shinyUI(
  navbarPage(
    "ML Webapp",
    tabPanel(
      "Flowcharts",
      img(src="fig1.png", align="left")
    ),
    tabPanel(
      "Multi-Modal Integration",
      sidebarLayout(
        sidebarPanel(
          selectizeInput(
            "dataset1", label="Dataset 1", choices=list("GSE21653"),
            multiple=F, options=list(create=TRUE, placeholder="Select Datasets")
          ),
          selectizeInput(
            "dataset2", label="Dataset 2", choices=list("GSE21654"),
            multiple=F, options=list(create=TRUE, placeholder="Select Datasets")
          ),
          selectizeInput(
            "color1", label="Color 1", choices=list(), multiple=F
          ),
          selectizeInput(
            "color2", label="Color 2", choices=list(), multiple=F
          ),
          selectizeInput("cluster", label="Cluster Method",
                         choices=list("PAM", "K-Means", "KNN")
          ),
          selectizeInput("method", label="Integration Method",
                         choices=list("BOMA", "JAMIE", "Unioncom",
                                      "Non-Linear Manifold Alignment",
                                      "Linear Manifold Alignment",
                                      "Non-Linear Manifold Warping",
                                      "Manifold Warping", "CCA")
          ),
          imageOutput("butterfly"),
        ),
        mainPanel(
          plotOutput("integrated1"),
          plotOutput("integrated2"),
          strong("Dataset 1:"),
          textOutput("summary1"),
          strong("Dataset 2:"),
          textOutput("summary2"),
        ),
      )
    ),
    tabPanel(
      "Cross-Modal Imputation",
      sidebarLayout(
        sidebarPanel(
          selectizeInput(
            "dataset1", label="Dataset 1", choices=list("GSE21653"),
            multiple=F, options=list(create=TRUE, placeholder="Select Datasets")
          ),
          selectizeInput(
            "dataset2", label="Dataset 2", choices=list("GSE21653"),
            multiple=F, options=list(create=TRUE, placeholder="Select Datasets")
          ),
          selectizeInput(
            "color", label="Color", choices=list("cell_type", "time"), multiple=F
          ),
          selectizeInput("cluster", label="Cluster Method",
                         choices=list("PAM", "K-Means", "KNN")
          ),
          selectizeInput("method", label="Imputation Method",
                         choices=list("MOFA+", "JAMIE", "Polarbear")
          )
        ),
        mainPanel( plotOutput("sampleGraph2") )
      )
    ),
    selected="Multi-Modal Integration"
  )
)

# Server
server <- function(input, output, session) {
  # DATA GENERATION
  r_data <- reactive({
    set.seed(42)
    data.frame(x=rnorm(1000), y=rnorm(1000))
  })
  
  gse1 <- reactive({ # Crashes on invalid entry
    return(t_gse1)
    getGEO(input$dataset1, GSEMatrix = TRUE)[[1]]
  })
  
  gse2 <- reactive({
    return(t_gse2)
    getGEO(input$dataset2, GSEMatrix = TRUE)[[1]]
  })
  
  gse_to_data <- function(gse) {
    samp = phenoData(gse)$geo_accession[c(1:10)] # TODO: Remove limit
    data = NULL
    for (s in samp) {
      s_gse = getGEO(s, GSEMatrix = TRUE)
      if (is.null(data)) {
        data = as.data.frame(s_gse@dataTable@table[,2])
        rownames(data) = s_gse@dataTable@table[,1]
      }
      else
        data = cbind(data, s_gse@dataTable@table[,2])
    }
    colnames(data) = samp
    data
  }
  
  gse1_data <- reactive({
    return(t_gse1_data) # testing
    gse_to_data(gse1())
  })
  
  gse2_data <- reactive({
    return(t_gse2_data) # testing
    gse_to_data(gse2())
  })
  
  # DATA PROCESSING
  integrated_data <- reactive({
    mat1 = gse1_data()
    mat2 = gse2_data()
    aligned = ManiNetCluster(
      t(mat1), t(mat2),
      nameX='s1',nameY='s2',
      corr=Correspondence(matrix=cor(mat1, mat2)),
      d=as.integer(2),
      method='nonlinear manifold aln',
      k_NN=as.integer(3),
      k_medoids=as.integer(6)
    )
  })
  
  # OUTPUTS
  output$summary1 <- renderText({
    experimentData(gse1())@other$summary
  })
  
  output$summary2 <- renderText({
    experimentData(gse2())@other$summary
  })
  
  output$butterfly <- renderImage({
    # TODO - generate through python
    list(src="www/ex_butterfly.png", width="100%")
  }, deleteFile=FALSE)
  
  output$integrated1 <- renderPlot({
    gse = gse1()
    data = integrated_data()
    data = data[data$data == 's1',]
    cl = phenoData(gse)@data[,input$color1][c(1:10)] # TODO: Remove limit
    data = cbind(as.data.frame(data), Group=cl)
    pl = ggplot(data=data, mapping=aes(x=Val0, y=Val1, color=Group))
    pl = pl + geom_point()
    pl
  })
  
  output$integrated2 <- renderPlot({
    gse = gse2()
    data = integrated_data()
    data = data[data$data == 's2',]
    cl = phenoData(gse)@data[,input$color2][c(1:10)] # TODO: Remove limit
    data = cbind(as.data.frame(data), Group=cl)
    pl = ggplot(data=data, mapping=aes(x=Val0, y=Val1, color=Group))
    pl = pl + geom_point()
    pl
  })
  
  output$sampleGraph <- renderPlot({
    pl = ggplot(data=r_data(), mapping=aes(x=x, y=y))
    pl = pl + geom_point(colour='red')
    pl
  })
  
  output$sampleGraph2 <- renderPlot({
    pl = ggplot(data=r_data(), mapping=aes(x=x, y=y))
    pl = pl + geom_point(colour='red')
    pl
  })
  
  observeEvent(gse1(), updateSelectizeInput(session, "color1", choices=colnames(phenoData(gse1())@data)))
  observeEvent(gse2(), updateSelectizeInput(session, "color2", choices=colnames(phenoData(gse2())@data)))
}

# Run
shinyApp(ui=ui, server=server)
