# NOTE: Crashes every other launch because of retriculate
# Libraries
library(ggplot2)
library(htmlwidgets)
library(shiny)
library(shinyjs)

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

# UI
ui <- shinyUI(
  navbarPage(
    "ML Webapp",
    tabPanel(
      "Integration",
      sidebarLayout(
        sidebarPanel(
          selectizeInput(
            "dataset", label="Data", choices=list("D1", "D2", "D3", "D4", "D5"),
            multiple=T, options=list(create=TRUE, placeholder="Select Datasets")
          ),
          selectizeInput(
            "color", label="Color", choices=list("cell_type", "time"), multiple=F
          ),
          radioButtons("cluster", label="Cluster Method",
                       choices=list("PAM", "K-Means", "KNN")
          ),
          radioButtons("method", label="Integration Method",
                       choices=list("BOMA", "JAMIE", "Unioncom",
                                    "Non-Linear Manifold Alignment",
                                    "Linear Manifold Alignment",
                                    "Non-Linear Manifold Warping",
                                    "Manifold Warping", "CCA")
          )
        ),
        mainPanel( plotOutput("sampleGraph") )
      )
    ),
    tabPanel(
      "Imputation",
      sidebarLayout(
        sidebarPanel(
          selectizeInput(
            "Dataset", label="Data", choices=list("D1", "D2", "D3", "D4", "D5"),
            multiple=T, options=list(create=TRUE, placeholder="Select Datasets")
          ),
          selectizeInput(
            "color", label="Color", choices=list("cell_type", "time"), multiple=F
          ),
          radioButtons("cluster", label="Cluster Method",
                       choices=list("PAM", "K-Means", "KNN")
          ),
          radioButtons("method", label="Imputation Method",
                       choices=list("MOFA+", "JAMIE", "Polarbear")
          )
        ),
        mainPanel( plotOutput("sampleGraph2") )
      )
    ),
  )
)

# Server
server <- function(input, output, session) {
  r_data <- reactive({
    set.seed(42)
    data.frame(x=rnorm(1000), y=rnorm(1000))
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
}

# Run
shinyApp(ui=ui, server=server)

