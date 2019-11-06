ui <- fluidPage(
  theme = shinytheme("united"),
  shinyjs::useShinyjs(),
  introjsUI(),
  
  
  # App title ----
  titlePanel(HTML("<div><b>Complex+ <b> <small>: Complex Aided Decision-Making for the Study of Protein Complexes</small></div>")),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      # must include in UI
      strong(actionButton("help", "Press for introduction", class = "btn-primary")),
      br(),br(),
      # Input: Select the random distribution type ----
      introBox(
        h5(radioButtons("species", "1. Dataset:",c("Homo Sapiens" = "hm"
                                                   #, "Saccharomyces cerevisiae" = "sc"
                                                   ))),
        data.step = 1,
        data.intro = "The first step is to choose your dataset"
      ),
      br(),
      # Help text 
      h5(strong("2. Choose protein complex:")),
      introBox(
        helpText(textOutput("complex")),
        data.step = 3,
        data.intro = "The selected complex and its member will show up here."
      ),
      br(),
      # Add Options
      h5(strong("3. Options:")),
      fluidPage(
        introBox(
          h5(strong(sliderInput("num", "Number of close interactors:", 30, min = 5, max = 100))),
          h5(strong(uiOutput("CalibrateOptions"))),
          data.step = 4,
          data.intro = "1. Select the number of potential interactors you are interested to see.\n 2. Add other proteins to the complex."
        )),
      introBox(
        strong(actionButton("go", "Calculate scores!")),
        data.step = 5,
        data.intro = "Calculate discriminant score of interactors."
      )
      #,
      #strong(actionButton("calibrate", "calibrate the result!"))
      
      , width = 3),
    
    # Main panel for displaying outputs ----
    mainPanel(
      tabsetPanel(type = "tabs",
                  
                  tabPanel("Protein Complex",
                           br(),
                           fluidPage(
                             introBox(
                               DT::dataTableOutput("table"),
                               data.step = 2,
                               data.intro = "Type your complex of interest in the search box on top right side of the table and click on the related record in table"
                             ),
                             br(),br(),
                             visNetworkOutput("complexMembers"),
                             br(),br(),
                             
                             htmlOutput("Description"),br(),br(),
                             fluidRow(
                               plotOutput(outputId="Schematics", height = "auto")
                             ), 
                             br(),br(),
                             fluidRow(plotOutput("DiseasePlot"))
                           )
                  ),
                  tabPanel("Interactors Network", 
                           br(),
                           conditionalPanel(
                             condition = "input.go",
                             downloadButton("downloadData", "Download result"),
                             fluidRow(
                               column(7, offset = 1,
                                      br(),br(),
                                      introBox(
                                        visNetworkOutput("resultNetwork", width = 1100, height = 700),
                                        data.step = 6,
                                        data.intro = "Check the protein interactors network")),
                               column(2, offset = 1,
                                      br(),br(),
                                      htmltools::div(style = "display:inline-block",
                                                     br()
                                                     ,plotlyOutput("plot", width = 200, height = 700)))),
                             br(),br(),
                             introBox(
                               br(),br(),br(),
                               div(DT::dataTableOutput("table3"), style = "font-size:90%"),
                               data.step = 7,
                               data.intro = "Search physical interaction of protein of interest and view the paper by clicking on related record")
                           ),
                           conditionalPanel(
                             condition = "input.table3_rows_selected>0",
                             htmlOutput("frame"))
                  )
      )
      
      , width = 9)
    , fluid=TRUE)
)
