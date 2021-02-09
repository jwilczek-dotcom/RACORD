#..............................................................................................
# R A C O R D  -  UI
#..............................................................................................
# Last update: 2021/02/04

library(shiny)
library(rgl)

width.adj <- c(2, 10)                       # width panels adjustement (side panel, main panel)
mainPanel.adj <- c("800px", "520px")        # main panel adjustement (width, height)


shinyUI(
  navbarPage("RACORD",
    tabPanel("Best-match Finder",
       fluidPage(
         sidebarLayout(
           sidebarPanel(
             strong('1. Load sources and targets'),
             fileInput('fileSources', label=NULL, accept=c('.txt')),
             fileInput('fileTargets', label=NULL, accept=c('.txt')),
             checkboxInput('cbxBreakProfile', 'Virtual fragmentation'),
               conditionalPanel(
                 condition = "input.cbxBreakProfile",
                 conditionalPanel(
                   condition = "input.cbxBreakProfile",
                   sliderInput('sldBreaks', 'Number of fragments', min=2, max=10, step=1, value=4),
                   sliderInput('sldMinPerc', 'Minimum % size of one fragment', min=5, max=50, step=1, value=25),
                   numericInput('numPtsSetSeed', 'Set seed', min = 1, max = 100000000, value = 1234),
                   checkboxInput('cbxOnlyRim', 'Extract only rim', value=F),
                   actionButton('btnBreakPottery', label='Break profile')
                 )
             ),
             br(),
             br(),
             strong('2. Data preparation'),
             selectInput('selNormalisation', label=NULL,
                         c("Select normalisation" = "",
                           "Real radius" = "Normalisation_radius",
                           "Rim units" = "Normalisation_rim"),
                         selected="Normalisation_radius"
             ),
             selectInput('selSampling', label=NULL,
                         c("Select sampling" = "",
                           "Segments" = "Sampling_segments",
                           "Equidistant" = "Sampling_equidistant"),
                         selected="Sampling_equidistant"
             ),
             checkboxInput('cbxAdvancedSampling', 'Advanced sampling'),
             conditionalPanel(
               condition = "input.cbxAdvancedSampling",
               conditionalPanel(
                 condition = "input.selSampling == 'Sampling_segments'",
                 sliderInput('sldNbOfPoints', label='Number of points', min=20, max=600, step=20, value = 20)
               ),
               conditionalPanel(
                 condition = "input.selSampling == 'Sampling_equidistant'",
                 sliderInput('sldPtsDistance', label='Distance between points', min=0.01, max=1, step=0.01, value = 0.01)
               )
             ),
             br(),
             br(),
             strong('3. Method selection'),
             selectInput('selMethod', label=NULL,
                         c("Select Matching method" = "",
                           "DCT (rim)" = "selDCT",
                           "RDP (rim)" = "selRDP",
                           "RTC (rim)" = "selRTC",
                           "Bezier (rim)" = "selBEZ",
                           "ICP (rim)" = "selICP_Rim",
                           "ICPz (all)" = "selICP_Rigid",
                           "ICP (all)" = "selICP",
                           "All Methods (rim)" = "selALL"
                           ),
                         selected="selRDP"
             ),
             conditionalPanel(
               condition = "input.selMethod == 'selDCT'",
               checkboxInput('cbxAdvancedDCT', 'Advanced options'),
               conditionalPanel(
                 condition = "input.cbxAdvancedDCT",
                 sliderInput('sldDCT_nbHarm', 'Number of harmonics', min = 2, max = 100, step = 1, value=20)
               )
             ),
             conditionalPanel(
               condition = "input.selMethod == 'selRDP'",
               checkboxInput('cbxAdvancedRDP', 'Advanced options'),
               conditionalPanel(
                 condition = "input.cbxAdvancedRDP",
                 sliderInput('sldRDP_nbPts', 'Number of points', min = 3, max = 100, step = 1, value=20)
               )
             ),
            conditionalPanel(
              condition = "input.selMethod == 'selRTC'",
              checkboxInput('cbxAdvancedRTC', 'Advanced options'),
              conditionalPanel(
                condition = "input.cbxAdvancedRTC",
                  sliderInput('sldRTC_W', 'Function weights (Rs,Ts,Ks)', min = 0, max = 1, step = 0.01, dragRange = T, value=c(0.33, 0.66))
              )
            ),
            conditionalPanel(
              condition = "input.selMethod == 'selICP_Rigid'",
              checkboxInput('cbxAdvancedICP_Rigid', 'Advanced options'),
              conditionalPanel(
                condition = "input.cbxAdvancedICP_Rigid",
                sliderInput('sldICP_rigid_z', 'Translation parameter (z)', min = 0.01, max = 1, step = 0.01, value=0.1)
              )
            ),
            conditionalPanel(
              condition = "input.selMethod == 'selICP'",
              checkboxInput('cbxAdvancedICP', 'Advanced options'),
              conditionalPanel(
                condition = "input.cbxAdvancedICP",
                               strong("Optimisation parameters (r, theta, size)"),
                               sliderInput('sldICP_r', label=NULL, value=c(-1,1), min=-10, max=10, step=0.1),
                               sliderInput('sldICP_theta', label=NULL, value=c(-5,5), min=-20, max=20, step=1),
                               sliderInput('sldICP_size', label=NULL, value=c(0.95,1.05), min=0.5, max=1.5, step=0.05),
                               
                               strong("SAN parameters (MaxIter, InitTemp, MaxTemp)"),
                               sliderInput('sldICP_MaxIter', label=NULL, value=2000, min=100, max=10000, step=100),
                               sliderInput('sldICP_Temp', label=NULL, value=100, min=100, max=10000, step=10),
                               sliderInput('sldICP_Tmax', label=NULL, value=40, min=10, max=1000, step=10)
              )
            ),
            actionButton('btnRun', label='RUN'),
            br(),
            br(),
            selectInput('selShow', label=NULL,
                       c("Show Results" = "",
                         "Analyse" = "selShowAnalyse",
                         "Results" = "selShowResults"),
                       selected="selShowAnalyse"),
            
            br(),
            br(),
            strong('4. Attribution'),
            br(),
            checkboxInput('cbxLabelling', 'Advanced options'),
            conditionalPanel(
              condition = "input.cbxLabelling",
              sliderInput('sldK', 'Number of best matches', min=1, max=10, step=1, value=1),
              selectInput('selLabellingMethod', label=NULL,
                          c("All the same class" = "selLabelSAME",
                            "Majority voting" = "selLabelMAJOR"
                          ),
                          selected="selLabelMAJOR"
              )
            ),
            actionButton('btnRefreshLabelling', label='Refresh'),
            
            width=width.adj[1]),

           mainPanel(
             splitLayout(
               DT::dataTableOutput('dataTableSources'),
               DT::dataTableOutput('dataTableTargets'),
               plotOutput('plotTemp1', height = 500)
             ),
             splitLayout(
                 plotOutput('plotTemp2'),
                 plotOutput('plotTemp3'),
                 plotOutput('plotTemp4')
             ),
             width=width.adj[2]
           )
         )
       )
    )
  )
)
