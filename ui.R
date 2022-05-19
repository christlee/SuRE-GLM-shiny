library(plotly)

shinyUI(
  navbarPage("SuRE GLM",
    tags$head(tags$script(src="https://cdn.plot.ly/plotly-latest.min.js")),
    tabPanel("How-To",
        sidebarLayout(
            sidebarPanel(
                HTML(paste('<b><h1>How to:</h1>',
                           '<h3>Inspect a promoter:</h3>',
                           '<ol>',
                           '  <li>Go to the "inspect promoter" tab</b>',
                           '    <dd>- The inspect promoter tab allows you to inspect promoter activity around a position of interest',
                           '    <dd>- The predict fragments tab allows you to predict expression for a list of fragments (by genomic coordinate)',
                           '    </br></br>',
                           '  <b><li>Enter a specific locus in the "position of interest" field using:</li></b>',
                           '    <dd>- gene symbol (e.g. NUP214)',
                           '    <dd>- ensembl/gencode gene id (e.g. ENSG00000126883[.16])',
                           '    <dd>- ensemble/gencode transcript id (e.g. ENST00000359428[.5])',
                           '    <dd>- position on the genome (hg19) in the form of [chromosome]:[position]:[strand] (used as center)</br></br>',
                           '  <b><li>Press Submit</li></br>',
                           '  <li>Select a promoter fragment to predict its expression by:</li></b>',
                           '    <dd>- Clicking on a point within the triangle plot</dd>',
                           '    <dd>- Selecting a region in the coefficient plot <i>(click + drag)</i></dd>',
                           '    <dd>- Filling in promoter fragment coordinates <i>(and pressing update)</i></dd></br>')),
            ),
            mainPanel(
                # HTML(paste('<img src="https://github.com/christlee/SuRE-GLM-shiny/raw/main/Tutorial1.png">',
                #            "<b><h3>Extra information</h3></b>",
                #            '<img src="https://github.com/christlee/SuRE-GLM-shiny/raw/main/Tutorial2.png">'))
                HTML('<img src="https://github.com/christlee/SuRE-GLM-shiny/raw/main/Tutorial2.png">')
            )

        )),
    tabPanel("Inspect promoter",
      # tags$style("#hg19_selection {font-size:12px;}"),
      # tags$style("#hg38_selection {font-size:12px;}"),
      # App title ----
      # Sidebar layout with input and output definitions ----
      fluidRow(
        # Sidebar panel for inputs ----
        column(2,
            radioButtons('lib', 'SuRE Library',
                         choices=c('K562', 'HEPG2')),
            checkboxInput("advanced", strong("show advanced options"),
                          value=FALSE,width="100%"),
            conditionalPanel(condition = "input.advanced == 1",
              # Input: Slider for the number of bins ----
              sliderInput("cutoff42",
                    label = "Promoter sizes to consider:",
                    min = 1, max = 600, value = 400),
              # Input: Slider for the number of bins ----
              sliderInput(inputId = "window",
              label = "region up and downstream from POI",
              min = -2000,
              max = 2000,
              value = c(-1000,1000)
              ),
              checkboxInput("show_minus", strong("show anti-sense orientation"),
                            value=FALSE,width="100%"),
          ),
          helpText("Use either gene symbol, ensembl ID (using gencode v27) or",
                   br(),
                   "[chromosome]:[position]:[strand] (hg19)."),
          textInput("ROI", h4("Position of interest"),
                    value = "e.g. NUP214, ENSG00000126883[.16], chr9:134000948:+"),
          actionButton("go", "Submit")
        ),
        column(1, plotOutput(outputId = "legend", height=150)),

        # Main panel for displaying outputs ----q
        column(6,
              uiOutput("text_warning"),
              fluidRow(column(6,uiOutput("text_plus")),
                       column(4,uiOutput("ucsc"))),
              # checkboxInput("show_plus", strong("sense orientation"),
              #               value=TRUE, width="100%"),
              plotOutput(outputId = "trianglePlot", click="plot_click", width=800, height=300),
              # plotOutput(outputId = "peakPlot", click="plot_click_peak", width=800, height=100),
              plotOutput(outputId = "flatPlot", brush="plot_brush", width=800, height=200),

              plotlyOutput(outputId="frame", width=810, height='100'),

              conditionalPanel(condition = "input.show_minus == 1",
                  uiOutput("text_minus"),
                  plotOutput(outputId = "trianglePlot_rev", click="plot_click" , width=800, height=300),
                  # plotOutput(outputId = "peakPlot_rev", click="plot_click_peak", width=800, height=100),
                  plotOutput(outputId = "flatPlot_rev", brush="plot_brush", width=800, height=200),
              ),

              fluidRow(column(6, uiOutput("hg19_sel")),
                       column(4, div(tableOutput('selection'), style="font-size:150%"),
                                 uiOutput("help_selection"))),



          ),
    column(3,
    ))),
    tabPanel("Predict fragments",
      sidebarLayout(
       # Sidebar panel for inputs ----
       sidebarPanel(
         radioButtons('hg_version', 'Human Genome version:',
                      choices=c('hg19', 'hg38')),
         # Input: Select a file ----
         fileInput("fragment_file", "Choose BED File",
             multiple = TRUE,
             accept = c("text/csv",
                        "text/comma-separated-values,text/plain",
                        ".bed")
          ),
          textAreaInput("text_fragments", "Fragment locations (bed format)",
                        placeholder=paste0("chr13\t32889034\t32889744\tProm_delta_1\t.\t+\n",
                                           "chr13\t32889421\t32889744\tProm_delta_2\t.\t+")),
          actionButton("go_fragments", "Submit", )

        ),
        mainPanel("",
          DT::dataTableOutput(outputId = "fragment_result"),
          uiOutput("download")
        )
      )
    )
  )
)
