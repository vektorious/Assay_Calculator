#Assay_Calculator 130517
#version 1.4
library(shiny)
library(shinydashboard)

dashboardPage(skin = "green",
  dashboardHeader(
    title = "Assay Calculator"  
    
  ),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Raw Data", tabName = "rawdata", icon = icon("calculator"), selected = TRUE),
      menuItemOutput("menu")
#      menuItem("Download", tabName = "download", icon = icon("download"))
    ),
    fileInput("data_file", label = h4("Data file input"), accept = c(".xls", ".xlsx")),
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
              radioButtons(inputId="assay_type", label=h4("Assay Type"), 
                          choices = list("Calcium Assay" = 1, "ROS Assay" = 2),
                          selected = 1),   
              uiOutput("ui.settings1"),
              hr(),
              checkboxInput("file_name", label = "write filename as header", value = TRUE),
              uiOutput("ui.settings5"),
              hr(),
              uiOutput("ui.settings2"),
              uiOutput("ui.settings3"),
              uiOutput("ui.settings4"),
              uiOutput("ui.downloaddata"),
              uiOutput("ui.downloadplot")
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
                           choices = list("vertical bars" = 1, "horizontal bars" = 2),
                           selected = 1),
              radioButtons(inputId="bar_columns", label=h4("Plot Arrangement"), 
                           choices = list("vertical alignment" = 1, "horizontal alignment" = 2),
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
                                   choices = list("Show SD (takes time!)" = 1))
                )
          )
        )
      )
    )
  )
