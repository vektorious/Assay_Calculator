#test_app 180317
library(shiny)
library(ggplot2)
library(XLConnect)
library(plyr)
library(reshape2)
library(grid)
library(shinydashboard)

shinyServer(function(input, output) {

  calculate <- reactive({
    
    ##################### functions #####################
    
    #normalising one well
    norm.well <- function(x){
      cs <- rev(cumsum(rev(x))); #cumulative sum
      xn <- x/(10*(cs))*1000;
      return(xn)
    }
    
    ####################################################
    
    if (input$assay_type == 2){
      dc <- 0
    }
    else { 
      dc <- input$dc
    }
    
    inputfile <- input$data_file
    
    if (is.null(inputfile))
      return(NULL)
    
    rawdata <- readWorksheetFromFile(inputfile$datapath, sheet=1)
    colnames(rawdata) <- gsub("M.", "", colnames(rawdata))
    rawdata[is.na(rawdata)]=0
  
    
    #normdata <- apply(rawdata, 2, norm.well) #w/o dependencies 
    normdata <- colwise(norm.well)(rawdata) # colwise needs library(plyr)
    
    ms_rawdata <- rawdata[1:(nrow(rawdata)-dc),] # isolates MS values of raw data
    ms_normdata <- normdata[1:(nrow(normdata)-dc),] # isolates MS values of normalized data
    
    ms_rawdata$time <- seq(0, (nrow(ms_rawdata)-1)*10, 10) # adds time values in seconds
    ms_normdata$time <- seq(0, (nrow(ms_normdata)-1)*10, 10) # adds time values in seconds
    
    ms_rd_melted <- melt(ms_rawdata, id.vars = "time", variable.name = "well") # formatting for ggplot
    ms_nd_melted <- melt(ms_normdata, id.vars = "time", variable.name = "well") # formatting for ggplot
    
    calc_output <- list("ms_rawdata" = ms_rawdata, 
                        "ms_normdata" = ms_normdata, 
                        "ms_rd_melted" = ms_rd_melted, 
                        "ms_nd_melted" = ms_nd_melted,
                        "yaxis_max" = max(ms_nd_melted$value))
    return(calc_output)
    
  })

  get_layout <- reactive({
    inputlayout <- input$layout_file
    
    if (is.null(inputlayout))
      return(NULL)
    
    raw_plate_layout <- readWorksheetFromFile(inputlayout$datapath, sheet=1)
    
    plate_layout <- mutate(raw_plate_layout,
                           Row=as.numeric(match(toupper(substr(well, 1, 1)), LETTERS)),
                           Column=as.numeric(substr(well, 2, 5)))
    
    return(plate_layout)
    
  })
  
  plate_layout <- reactive({
    
    plate_layout <- get_layout()
    
    if (is.null(plate_layout))
      return(NULL)
    
    glayout <- ggplot(data=plate_layout, aes(x=Column, y=Row)) +
      geom_point(size=9) +
      geom_point(size=7, aes(colour = elicitor)) +
      geom_point(size=9, aes(alpha = genotype)) +
      scale_alpha_discrete(range = c(0.4, 0.05)) +
      labs(title="Plate Layout")
    glayout <- glayout +
      coord_fixed(ratio=(13/12)/(9/8), xlim=c(0.8, 12.2), ylim=c(0.6, 8.4)) +
      scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) +
      scale_x_continuous(breaks=seq(1, 12))
    glayout <- glayout + theme_bw() + guides(colour = guide_legend(title = "Elicitors", ncol = 2, byrow = TRUE), alpha = guide_legend(title = "Genotypes", ncol = 2, byrow = TRUE)) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.key=element_blank()
      )
    return(glayout)
    
  })
  
  mean_max_calculate <- reactive({
    
    plate_layout <- get_layout()
    
    if (is.null(plate_layout))
      return(NULL)
    
    raw_plate_layout <- plate_layout[,1:4]
    raw_plate_layout <- raw_plate_layout[complete.cases(raw_plate_layout),] #delete rows with NAs
    
    elicitors <- unique(raw_plate_layout$elicitor)
    genotypes <- unique(raw_plate_layout$genotype)
    
    data <- calculate()
    
    if (is.null(data))
      return(NULL)
    
    ms_normdata <- data$ms_normdata
    
    all_means <- matrix(nrow = length(elicitors)*length(genotypes), ncol = 4, dimnames = NULL)
    all_means_graph <- matrix(nrow = length(elicitors)*length(genotypes)*length(ms_normdata[,1]), ncol = 5, dimnames = NULL)
    i = 1
    i2 = 1
    
    for(eli in elicitors){
      for(gen in genotypes){
        wells <- raw_plate_layout$well[raw_plate_layout$elicitor == eli & raw_plate_layout$genotype == gen]
        current_set <- ms_normdata[colnames(ms_normdata) %in% wells]
        
        current_rowmeans <- rowMeans(current_set)
        current_rowsds <- apply(current_set, 1, sd)
        current_max <- max(current_rowmeans[16:length(current_rowmeans)])
        current_sd <- current_rowsds[current_rowmeans==current_max]
        #print(length(c(eli, gen, current_max, current_sd)))
        #print(c(paste(gen, eli, sep = " "), current_max, current_sd))
        all_means[i,] <- c(eli, gen, current_max, current_sd)
        all_means_graph[i2:(i2+length(current_rowmeans)-1),1] <- rep(eli, length(current_rowmeans))
        all_means_graph[i2:(i2+length(current_rowmeans)-1),2] <- rep(gen, length(current_rowmeans))
        all_means_graph[i2:(i2+length(current_rowmeans)-1),3] <- current_rowmeans
        all_means_graph[i2:(i2+length(current_rowmeans)-1),4] <- current_rowsds
        all_means_graph[i2:(i2+length(current_rowmeans)-1),5] <- seq(0, (nrow(ms_normdata)-1)*10, 10)
        
        i = i+1
        i2 = i2 + length(current_rowmeans)
      }
    }
    
    total_output <- list("all_means" = all_means,
                         "all_means_graph" = all_means_graph)
    
    return(total_output)
    
  })
  
  well_plots <- reactive({
    
    data <- calculate()
    
    
    if (is.null(data))
      return(NULL)
    
    if (input$data_type == 1){
      chosen_set <- data$ms_rd_melted
    }
    
    if (input$assay_type == 2){
      chosen_set <- data$ms_rd_melted
    }
    
    else {
      chosen_set <- data$ms_nd_melted
    }
    
    g <- ggplot(chosen_set, aes(time, value)) + geom_line() + facet_wrap(~well, ncol = 12)
    
    g <- g + theme_bw() + theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      strip.text.x = element_blank(),
      strip.background = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(colour = "black")
    )
    
    #g <- g + geom_line(data = data$ms_normdata, aes(time, A1), color = "red")
    g <- g + annotate("text",
                      x = (max(chosen_set$time)*0.85),
                      y = (max(chosen_set$value, na.rm = TRUE)*0.9),
                      label = colnames(data$ms_normdata)[1:(length(data$ms_normdata)-1)],
                      size = 2.5)
    
    
    #if (data$yaxis_max != input$ylim){
    #  g <- g + ylim(0, input$ylim)
    #}
    
    return(g)

  })
  
  bar_plots_max <- reactive({
    
    plate_layout <- get_layout()
    
    if (is.null(plate_layout))
      return(NULL)
    
    raw_plate_layout <- plate_layout[,1:4]
    raw_plate_layout <- raw_plate_layout[complete.cases(raw_plate_layout),] #delete rows with NAs
    
    elicitors <- unique(raw_plate_layout$elicitor)
    genotypes <- unique(raw_plate_layout$genotype)
    
    complete_data <- mean_max_calculate()
    
    all_means <- complete_data$all_means
    
    num_all_means <- all_means[,3:4]
    class(num_all_means) <- "numeric" # sonst werden max und se nicht alszahlen erkannt...
    max_barplot <- data.frame(
      eli = all_means[,1],
      gen = factor(all_means[,2], levels = genotypes),
      maxima = num_all_means[,1],
      se = num_all_means[,2]
    )
    
    if (input$bar_columns == 2){
      bpmax <- ggplot(max_barplot, aes(fill=eli, y=maxima, x=eli)) + facet_wrap(~gen, ncol = 2) + ggtitle("Maxima with SD")
    } else {
      bpmax <- ggplot(max_barplot, aes(fill=eli, y=maxima, x=eli)) + facet_wrap(~gen, ncol = 1) + ggtitle("Maxima with SD")
      }
    
    if(input$bar_rotation == 2){bpmax <- bpmax + coord_flip() + theme(legend.position = "none")}
    
    bpmax <- bpmax + geom_bar(position="dodge", stat="identity")
    bpmax <- bpmax + geom_errorbar(aes(ymin=maxima-se, ymax=maxima+se),
                                   width=.2,                    # Width of the error bars
                                   position=position_dodge(.9))
    
    
    bpmax <- bpmax + labs(x="", y="L/Lmax") + theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5),
                                                    legend.title=element_blank(),
                                                    panel.background = element_rect(fill = "white", colour = "grey90"),
                                                    panel.grid.major = element_line(color = "grey90"),
                                                    strip.background = element_rect(fill = "grey90", colour = NA),
                                                    panel.grid.minor = element_line(colour = "grey90", size = 0.25))
    
    return(bpmax)
    
  })
  
  
  mean_graphs <- reactive({
    
    plate_layout <- get_layout()
    
    if (is.null(plate_layout))
      return(NULL)
    
    raw_plate_layout <- plate_layout[,1:4]
    raw_plate_layout <- raw_plate_layout[complete.cases(raw_plate_layout),] #delete rows with NAs
    
    elicitors <- unique(raw_plate_layout$elicitor)
    genotypes <- unique(raw_plate_layout$genotype)
    
    complete_data <- mean_max_calculate()
    
    all_means_graph <- complete_data$all_means_graph
    
    num_all_means_graph <- all_means_graph[,3:5]
    class(num_all_means_graph) <- "numeric" # sonst werden max und se nicht alszahlen erkannt...
    max_lineplot <- data.frame(
      eli = all_means_graph[,1],
      gen = factor(all_means_graph[,2], levels = genotypes),
      mean_value = num_all_means_graph[,1],
      se = num_all_means_graph[,2],
      time = num_all_means_graph[,3]
    )
    
    lpmax2 <- ggplot(max_lineplot, aes(color = gen,x=time, y=mean_value)) + geom_line() + facet_wrap(~eli) + 
      scale_color_manual(values=c("navy", "red3", "deeppink1", "greenyellow", "gray20", "green4")) + 
      ggtitle("Mean with SD") +
      geom_line(aes(x=time,y=mean_value-se), alpha=0.2)+ 
      geom_line(aes(x=time,y=mean_value+se), alpha=0.2)+
      geom_errorbar(aes(ymin=mean_value-se, ymax=mean_value+se), alpha=0.2,
                    width=.2,                    # Width of the error bars
                    position=position_dodge(.9))
    lpmax2 <- lpmax2 + labs(x="", y="L/Lmax") + theme(legend.title=element_blank(),
                                                      panel.background = element_rect(fill = "white", colour = "grey90"),
                                                      panel.grid.major = element_line(color = "grey90"),
                                                      strip.background = element_rect(fill = "grey90", colour = NA),
                                                      panel.grid.minor = element_line(colour = "grey90", size = 0.25),
                                                      legend.key=element_rect(fill='white'),
                                                      strip.text.x = element_text(size=10))
    return(lpmax2)
  })
  
  output$plot <- renderPlot({
    plot <- well_plots()
    show(plot)
  })
  
  output$plate_layout <- renderPlot({
    layout <- plate_layout()
    show(layout)
  })
  
  output$bar_max <- renderPlot({
    barplot <- bar_plots_max()
    show(barplot)
  })
  
  output$graph_mean <- renderPlot({
    mean_graph <- mean_graphs()
    show(mean_graph)
  })
  
  output$ui.downloaddata <- renderUI({
    if (is.null(calculate())) 
      return(NULL)
    if (input$assay_type == 2)
      return(NULL)
    downloadButton('downloadData', 'Download Data')
  })
  output$downloadData <- downloadHandler(
    
    filename = function() {paste0("norm_", gsub(unlist(strsplit(input$data_file$name, "[.]")[1])[-1], "", input$data_file$name), "xlsx")},
    content = function(file) {
      data <- calculate()
      #writeWorksheetToFile(file, data = data$ms_normdata, sheet = "normalized data")
      fname <- paste0("norm_", gsub(unlist(strsplit(input$data_file$name, "[.]")[1])[-1], "", input$data_file$name), "xlsx")
      wb <- loadWorkbook(fname, create = TRUE)
      createSheet(wb, name = "raw data")
      writeWorksheet(wb, data$ms_rawdata, sheet = "raw data")
      createSheet(wb, name = "normalized data")
      writeWorksheet(wb, data$ms_normdata, sheet = "normalized data")
      saveWorkbook(wb)
      file.rename(fname,file)
    
    }
  )
  
  output$ui.downloadplot <- renderUI({
    if (is.null(calculate()))
      return(NULL)
    downloadButton('downloadPlot', 'Download Plot')
  })
  output$downloadPlot <- downloadHandler(
    filename = function() {paste0("plot_", gsub(unlist(strsplit(input$data_file$name, "[.]")[1])[-1], "", input$data_file$name), "pdf")},
    content = function(file) {
      #well_plots()
      #ggsave(file, width = 29.7, height = 21.0, units = "cm", limitsize = FALSE)
      wellcurves <- well_plots()
      platelayout <- plate_layout()
      bar_plots <- bar_plots_max()
      graph_mean <- mean_graphs()
      
      pdf(file, width = 29.7, height = 21.0, paper = "a4r")
      invisible(print(wellcurves))
      if (is.null(platelayout) == FALSE){invisible(print(platelayout))}
      if (is.null(bar_plots) == FALSE){invisible(print(bar_plots))}
      if (is.null(graph_mean) == FALSE){invisible(print(graph_mean))}
      dev.off()
      
    }
  )
  
  output$ui.settings1 <- renderUI({
    if (input$assay_type == 2){
      HTML("At the moment it is only possible to view and download raw data graphs for ROS measurements! <br/> <br/> <div><span style='color:red'>Keep in mind that the data has to be in the 'Luminoscan' format!</span></div>")
    } else{
    radioButtons(inputId="data_type",
                 label="Show graphs for:",
                 choices = list("raw data" = 1, "normalized data" = 2), selected = 2)
    } 
  })
  
  output$ui.settings2 <- renderUI({
    if (input$assay_type == 2){
      return(NULL)
    }
    sliderInput("dc", "Discharge measurement points",min=0, max=20, value=15)
  })
  
  
  output$menu <- renderMenu({
    
    inputlayout <- input$layout_file
    
    sidebarMenu(
      menuItem("Raw Data", tabName = "rawdata", icon = icon("calculator"), selected = TRUE),
      if (is.null(inputlayout)==FALSE){
        menuItem("Data Summary", tabName = "data_summary", icon = icon("bar-chart"),
                menuSubItem("Mean Maxima", tabName = "norm_data1"),
                menuSubItem("Mean Kinetics", tabName = "norm_data2"))
      }
    )
  })
})