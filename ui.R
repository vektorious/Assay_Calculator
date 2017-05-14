#Assay_Calculator 130517
#version 1.4
library(shiny)
library(shinydashboard)

dashboardPage(
  dashboardHeader(
    title = "Assay Calculator"  
    
  ),
  
  dashboardSidebar(
    sidebarMenuOutput("menu"),
    
    fileInput("data_file", label = h4("Data file input"), accept = c(".xls", ".xlsx")),
    fileInput("layout_file", label = h4("Layout file input (optional)"), accept = c(".xls", ".xlsx"))
  ),
  
  dashboardBody(
    tabItems(
      tabItem(tabName = "rawdata",
        fluidRow(
          box(title = "Well Curves",
              width = 9,
              collapsible = TRUE,
              plotOutput("plot", height = 600)),
          box(title = "Well Curves Settings",
              width = 3,
              collapsible = TRUE,
              radioButtons(inputId="assay_type", label=h4("Assay Type"), 
                          choices = list("Calcium Assay" = 1, "ROS Assay" = 2),
                          selected = 1),   
              uiOutput("ui.settings1"),
              uiOutput("ui.settings2"),
              uiOutput("ui.downloaddata"),
              uiOutput("ui.downloadplot")
          )
        ),
        fluidRow(
          box(title = "Plate Layout",
              width = 12,
              collapsible = TRUE,
              plotOutput("plate_layout", height = 300))
        )
      ),
      
      tabItem(tabName = "norm_data1",
        fluidRow(
          box(title = "Maxima of the Mean with SD",
              width = 9,
              collapsible = TRUE,
              plotOutput("bar_max", height = 600)
              ),
          box(title = "Plot Settings",
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
                collapsible = TRUE,
                plotOutput("graph_mean", height = 600)
                ),
            box(title = "Plot Settings")
          )
        )
      )
    )
  )
