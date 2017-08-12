#Assay_Calculator 130517
#version 1.4
library(shiny)
library(shinydashboard)
options(java.parameters = "-Xss2560k")

dashboardPage(skin = "green",
  dashboardHeader(
    title = "Assay Calculator"  
    
  ),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Instructions", tabName = "instructions", icon = icon("commenting-o"), selected = FALSE),
      menuItem("Data Analysis", tabName = "analysis", icon = icon("calculator"), selected = TRUE,
               menuSubItem("Standart Analysis", tabName = "rawdata", icon = icon("file-o")),
               menuSubItem("Multiple Measurements", tabName = "multi", icon = icon("clone"))),
      menuItemOutput("menu1"),
      menuItem("Mapping Tools", tabName = "mapping", icon = icon("map-o"), selected = FALSE,
               menuSubItem("Combine Files", tabName = "combine", icon = icon("clone")),
               menuSubItem("Mapping Wellcurves", tabName = "mapping_wellcurves", icon = icon("calculator"))),
      menuItemOutput("menu2")

    ),
    fileInput("data_file", label = h4("Data file input"), accept = c(".xls", ".xlsx", ".csv")),
    fileInput("layout_file", label = h4("Layout file input (optional)"), accept = c(".xls", ".xlsx"))
  ),
  
  dashboardBody(
    tags$head(tags$style(HTML('
#      .skin-green .main-sidebar {
#        background-color: #666666;
#      }
#      .skin-green .sidebar-menu>li.active>a, .skin-green .sidebar-menu>li:hover>a {
#        background-color: #444444;
#      }
      .content-wrapper,
      .right-side {
        background-color: #ffffff;
      }
    '))),
    tabItems(
      tabItem(tabName = "instructions",
          HTML("<center><h1><strong>Welcome to Assay Calculator</strong></h1><h6>Version 2.1 (21.05.2017)</h6></center>"),
          br(),
          div(class="body", style="font-size:120%",
              HTML("<p>This web application is designed to help you to analyse your data.
                   <br> 
                   For now it can process data you obtained from aequorin luminescence measurements (<strong>'Calcium Assay'</strong>) or ROS measurements (<strong>'ROS Assay'</strong>)
                   <br>
                   <br>
                   To view well plots and to start the calculations open the 'Data Analysis' tab.
                   <strong><h3>Here are some things you have to keep in mind:</h3></strong>
                   <li>Input files have to be .xls or .xlsx files</li>
                   <li>Data has to be formatted accordingly</li>
                   <li>Plate layout is optional but required for ROS Assay normalization and mean/maxima calculation</li>
                   <li><strong>ROS layout: Put 'control' in elicitor column for measurements used for blanking</strong></li>
                   <li>Calculations take their time! Be patient!</li></p>"),
              HTML("<h3><strong>Download example data or empty templates:</strong></h3>"),
              tags$li(div(style="display:inline-block",
                    p("Get")),
                downloadLink('downloadCaExample', 'Calcium Assay data'),
                div(style="display:inline-block",
                    p("and")),
                downloadLink('downloadCaExample_layout', 'Calcium Assay plate layout')),
              tags$li(div(style="display:inline-block",
                    p("Get")),
                downloadLink('downloadROSExample', 'ROS Assay data'),
                div(style="display:inline-block",
                    p("and")),
                downloadLink('downloadROSExample_layout', 'ROS Assay plate layout')),
              tags$li(div(style="display:inline-block",
                    p("Get empty")),
                downloadLink('downloadDatatemplate', 'data template'),
                div(style="display:inline-block",
                    p("and")),
                downloadLink('downloadLayouttemplate', 'plate layout template')  
              )),
          hr(),
          div(style="display:inline-block",
            p("This web application was created by Alexander Kutschera (TU Munich) as a OpenPlantScience community project.")),
          div(style="display:inline-block",
              p("Please feel free and send me your ")),
          div(style="display:inline-block",
              a(href="mailto:alexander.kutschera@tum.de", "feedback!")
          ),
          div(style="display:inline-block",
              p("Assay Calculator is licenced with the ")),
          downloadLink('downloadLicense', 'GNU General Public License v3.0')



      ),
      tabItem(tabName = "rawdata",
        fluidRow(
          box(title = "Well Curves",
              solidHeader = TRUE,
              status = "success",
#              background = "black",
              width = 9,
              collapsible = FALSE,
              plotOutput("plot", height = 600)),
          box(title = "Well Curves Settings",
              solidHeader = TRUE,
              status = "success",
              width = 3,
#              collapsible = TRUE,
              radioButtons(inputId="assay_type", label="Assay Type:", 
                          choices = list("Calcium Assay" = 1, "ROS Assay" = 2),
                          selected = 1),
              uiOutput("ui.settings1"),
              hr(),
              h4("Options"),
              checkboxInput("file_name", label = "Write filename as header", value = TRUE),
              uiOutput("ui.settings5"),
              uiOutput("ui.settings2"),
              uiOutput("ui.settings3"),
              uiOutput("ui.settings4")
          )
        ),
        fluidRow(
          uiOutput("ui.plate_layout")
        )
      ),
      
      tabItem(tabName = "norm_data1",
        fluidRow(
          box(title = "Maxima of the Mean with SD",
              width = 9,
#              collapsible = TRUE,
              solidHeader = TRUE,
              status = "success",
              plotOutput("bar_max", height = 600)
              ),
          box(title = "Plot Settings",
              width = 3,
              solidHeader = TRUE,
              status = "success",
              radioButtons(inputId="bar_rotation", label=h4("Bar Rotation"), 
                           choices = list("Vertical Bars" = 1, "Horizontal Bars" = 2),
                           selected = 1),
              radioButtons(inputId="bar_columns", label=h4("Plot Arrangement"), 
                           choices = list("Vertical Alignment" = 1, "Horizontal Alignment" = 2),
                           selected = 1)
            )
          )
        ),
      tabItem(tabName = "norm_data2",
          fluidRow(
            box(title = "Mean Kinetics",
                width = 9,
#                collapsible = TRUE,
                solidHeader = TRUE,
                status = "success",
                plotOutput("graph_mean", height = 600)
                ),
            box(title = "Plot Settings",
                width = 3,
                solidHeader = TRUE,
                status = "success",
                checkboxGroupInput("settings_mean", label = h5("General Setting"), 
                                   choices = list("Show SD (takes time!)" = 1)),
                radioButtons(inputId="graph_sorting", label=h4("Sort graphs by"), 
                             choices = list("Elicitor" = 1, "Genotype" = 2),
                             selected = 1)
                )
          )
        ),


      tabItem(tabName = "combine",
          fluidRow(  
            box(title = "Combine raw files",
                solidHeader = TRUE,
                status = "success",
                width = 12,
                collapsible = TRUE,
                div(style="display: inline-block;vertical-align:top; width: 300px;",fileInput("Q1E_file", label = h4("Quadrant 1 measurement data input"), accept = c(".txt"))),
                div(style="display: inline-block;vertical-align:top; width: 300px;",fileInput("Q2E_file", label = h4("Quadrant 2 measurement data input"), accept = c(".txt"))),
                div(style="display: inline-block;vertical-align:top; width: 300px;",fileInput("Q3E_file", label = h4("Quadrant 3 measurement data input"), accept = c(".txt"))),
                div(style="display: inline-block;vertical-align:top; width: 300px;",fileInput("Q4E_file", label = h4("Quadrant 4 measurement data input"), accept = c(".txt"))),
                br(),
                div(style="display: inline-block;vertical-align:top; width: 300px;",fileInput("Q1D_file", label = h4("Quadrant 1 discharge data input"), accept = c(".txt"))),
                div(style="display: inline-block;vertical-align:top; width: 300px;",fileInput("Q2D_file", label = h4("Quadrant 2 discharge data input"), accept = c(".txt"))),
                div(style="display: inline-block;vertical-align:top; width: 300px;",fileInput("Q3D_file", label = h4("Quadrant 3 discharge data input"), accept = c(".txt"))),
                div(style="display: inline-block;vertical-align:top; width: 300px;",fileInput("Q4D_file", label = h4("Quadrant 4 discharge data input"), accept = c(".txt"))),
                br(),
                downloadButton('downloadCombined', 'Download File')

          )
        )
),
      tabItem(tabName = "mapping_wellcurves",
            fluidRow(
              box(title = "Well Curves Settings",
                  solidHeader = TRUE,
                  status = "success",
                  width = 12,
                  collapsible = TRUE,
                  h4("Upload assay data with wildtype plants for comparison"),
                  div(style="display: inline-block;vertical-align:top; width: 300px;",fileInput("WT_file", label = h6("Wildtype data input (optional)"), accept = c(".csv"))),
                  br(),
                  div(style="display: inline-block;vertical-align:mid; width: 150px;", checkboxInput("mean_overlay_WT", label = "plot WT mean", value = FALSE)),
                  div(style="display: inline-block;vertical-align:mid; width: 150px;",checkboxInput("mean_overlay_plate", label = "plot plate mean", value = FALSE))
              )
            ),          
              
          fluidRow(
            box(title = "Mapping Well Curves",
                solidHeader = TRUE,
                status = "success",
                #              background = "black",
                width = 12,
                collapsible = FALSE,
                plotOutput("mappingplot", height = 900))
          )

),

tabItem(tabName = "multi",
        fluidRow(
          box(title = "Well Curves Settings",
              solidHeader = TRUE,
              status = "success",
              width = 12,
              collapsible = TRUE,
              h4("Upload additional datasets for compairson"),
              div(style="display: inline-block;vertical-align:top; width: 300px;",fileInput("data2", label = h4("Data file input 2"), accept = c(".xls", ".xlsx"))),
              div(style="display: inline-block;vertical-align:top; width: 300px;",fileInput("data3", label = h4("Data file input 3"), accept = c(".xls", ".xlsx"))),
              div(style="display: inline-block;vertical-align:top; width: 300px;",fileInput("data4", label = h4("Data file input 4"), accept = c(".xls", ".xlsx"))),
              br(),
              div(style="display: inline-block;vertical-align:top; width: 300px;",fileInput("layout2", label = h4("Layout file input 2"), accept = c(".xls", ".xlsx"))),
              div(style="display: inline-block;vertical-align:top; width: 300px;",fileInput("layout3", label = h4("Layout file input 3"), accept = c(".xls", ".xlsx"))),
              div(style="display: inline-block;vertical-align:top; width: 300px;",fileInput("layout4", label = h4("Layout file input 4"), accept = c(".xls", ".xlsx")))
          )
        ),
          
          box(title = "Multiple Well Curves",
              solidHeader = TRUE,
              status = "success",
              #              background = "black",
              width = 12,
              collapsible = FALSE,
              plotOutput("multiple_wellplots")
              )
),


      tabItem(tabName = "download", 
             h3("Download Everything!"),
             box(title = "Download Data",
                 width = 4,
                 solidHeader = TRUE,
                 status = "success",
                 uiOutput("ui.downloaddata"),
                 checkboxGroupInput("settings_data_download1", label = h5("General Setting"),
                                    selected = "1",
                                    choices = list("Add graphs to Excel file" = 1))
                                                   #"Don't add mean data" = 2))
             ),
             box(title = "Download Plots",
                 width = 4,
                 solidHeader = TRUE,
                 status = "success",
                 uiOutput("ui.downloadplot")
                 #checkboxGroupInput("settings_plot_download1", label = h5("General Setting"), 
                 #                  choices = list("Don't add mean data" = 1))
             )
             
        )
      )
    )
  )
