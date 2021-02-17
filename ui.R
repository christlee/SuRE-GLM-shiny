

shinyUI(
    fluidPage(
      tags$style("#hg19_selection {font-size:9px;}"),
      tags$style("#hg38_selection {font-size:9px;}"),
      # App title ----
      titlePanel("Hello Shiny!"),
      # Sidebar layout with input and output definitions ----
      fluidRow(
        # Sidebar panel for inputs ----
        column(2,
          radioButtons('lib', 'SuRE Library',
                       choices=c('K562', 'HEPG2')),
          # Input: Slider for the number of bins ----
          conditionalPanel(
            condition = "input.lib == 'HT1080'",
            sliderInput("cutoff23",
                        label = "Promoter sizes to consider:",
                        min = 1, max = 2000, value = 1500)
          ),
          conditionalPanel(
            condition = "input.lib %in% c('K562', 'HEPG2')",
            sliderInput("cutoff42",
                        label = "Promoter sizes to consider:",
                        min = 1, max = 600, value = 400)
          ),
          # Input: Slider for the number of bins ----
          sliderInput(inputId = "binsize",
            label = "binsize used to display output",
            min = 1,
            max = 100,
            value = 10
          ),
          # Input: Slider for the number of bins ----
          sliderInput(inputId = "window",
            label = "region up and downstream from ROI",
            min = -2000,
            max = 2000,
            value = c(-1000,1000)
          ),

          helpText(paste0("Use either gene symbol, ensembl ID or chromosome ",
                          "position (using gencode v27).")),
          textInput("ROI", h4("Region of interest"),
                    value = "e.g. NUP214, ENSG00000126883[.1], chr9:134000948:+"),
          actionButton("go", "Submit")
        ),
        # Main panel for displaying outputs ----q
        column(8,
          #
          # Output: Histogram ----
          fluidRow(
              column(3,
                  fluidRow(textInput("hg19_selection",
                                     "ROI hg19:"),
                           textInput("hg38_selection",
                                     "ROI hg38:")),
                           actionButton("update", "Update")),
              column(8,
                  tableOutput('selection'))),

          h3("sense orientation"),
          plotOutput(outputId = "trianglePlot", click="plot_click", width=800, height=300),
          plotOutput(outputId = "peakPlot", click="plot_click_peak", width=800, height=100),
          plotOutput(outputId = "flatPlot", brush="plot_brush", width=800, height=200),
          htmlOutput(outputId = "jbrowse", width=800, height=200),
          h3("anti-sense orientation"),
          plotOutput(outputId = "trianglePlot_rev", click="plot_click", width=800, height=300),
          plotOutput(outputId = "peakPlot_rev", click="plot_click_peak", width=800, height=100),
          plotOutput(outputId = "flatPlot_rev", brush="plot_brush", width=800, height=200),
          # plotOutput(outputId = "trianglePlot2", click="plot_click", width=800, height=800),
        ),
        column(1,
          plotOutput(outputId = "legend", width=50, height=100)
        )
      )
  )
)
