library(shiny)
library(ggplot2)
library(XLConnect)
library(plyr)
library(reshape2)
library(grid)
library(shinydashboard)
library(gridExtra)
library(gtools)
library(ggsci)
options(java.parameters = "-Xss2560k")

##################### functions #####################

#normalising one well for calcium measurements

norm.well <- function(x){
  if(sum(x) == 0) return(x) #avoids errors in empty sets
  cs <- rev(cumsum(rev(x))); #cumulative sum
  xn <- x/(10*(cs))*1000; #correction factor 10 because of timepoints
  return(xn)
}

#converting raw .txt files of Luminoscan to complete dataframes

curateQ <- function(QE = NULL, QD = NULL){
  if(is.null(QE)) return(0)
  if(is.null(QD)) return(0)
  
  QE <- read.table(QE, sep = "\t", dec = ",", header = FALSE, na.strings = "", as.is = TRUE)
  QE <- QE[,!apply(QE, 2, function(x){all(is.na(x))})] 
  QE <- QE[!apply(QE, 1, function(x){all(is.na(x))}),]
  QE <- QE[,2:(ncol(QE))]
  
  QE_curated <- NULL
  index = 1
  for(i in seq(nrow(QE)/73)){
    if(is.null(QE_curated)) {
      QE_curated <- QE[index:(index+72),1:ncol(QE)]
    } else {
      QE_curated <- cbind(QE_curated, QE[index:(index+72),1:ncol(QE)])
    }
    index = index + 73
  }
  
  #QE_curated[1,] <- gsub("M.", "", QE_curated[1,])
  colnames(QE_curated) <- QE_curated[1,]
  newcolnames <- unlist(strsplit(colnames(QE_curated), ".", fixed = TRUE))
  colnames(QE_curated) <- newcolnames[!newcolnames=="M"]
  
  QE_curated <- QE_curated[2:nrow(QE_curated),]
  QE_curated <- apply(QE_curated, 2, type.convert, dec = ",")
  
  QD <- read.table(QD, sep = "\t", dec = ",", header = FALSE, na.strings = "", as.is = TRUE)
  QD <- QD[,!apply(QD, 2, function(x){all(is.na(x))})]
  QD <- QD[!apply(QD, 1, function(x){all(is.na(x))}),]
  QD <- QD[,2:(ncol(QD))]
  
  QD_curated <- NULL
  index = 1
  for(i in seq(nrow(QD)/7)){
    if(is.null(QD_curated)) {
      QD_curated <- QD[index:(index+6),1:ncol(QD)]
    } else {
      QD_curated <- cbind(QD_curated, QD[index:(index+6),1:ncol(QD)])
    }
    index = index + 7
  }

  
#  QD_curated[1,] <- gsub("M.", "", QD_curated[1,])
  
  colnames(QD_curated) <- QD_curated[1,]
  newcolnames <- unlist(strsplit(colnames(QD_curated), ".", fixed = TRUE))
  colnames(QD_curated) <- newcolnames[!newcolnames=="M"]
  
  QD_curated <- QD_curated[2:nrow(QD_curated),]
  QD_curated <- apply(QD_curated, 2, type.convert, dec = ",")
  
  Q_curated <- rbind(QE_curated, QD_curated)
  Q_curated_df <- data.frame(Q_curated, row.names = NULL)
  Q_curated_df[is.na(Q_curated_df)] = 0
  newcolnames <- unlist(strsplit(colnames(Q_curated_df), ".", fixed = TRUE))
  colnames(Q_curated_df) <- newcolnames[!newcolnames=="M"]
  
  return(Q_curated_df)
}

#get the plate layout out of the raw plate layout file

layout <- function(inputlayout){

  if (is.null(inputlayout))
    return(NULL)

  raw_plate_layout <- readWorksheetFromFile(inputlayout$datapath, sheet=1)

  plate_layout <- mutate(raw_plate_layout,
                        Row=as.numeric(match(toupper(substr(well, 1, 1)), LETTERS)),
                        Column=as.numeric(substr(well, 2, 5)))

  empty_wells <- raw_plate_layout$well[is.na(raw_plate_layout$elicitor) & is.na(raw_plate_layout$genotype)]

  output <- list("empty_wells" = empty_wells,
                "plate_layout" = plate_layout)
  return(output)
}

# calculate one dataset

calculate_data <- function(inputfile, dc, data_type, assay_type){

    if((tools::file_ext(inputfile[1]) == "xlsx")||(tools::file_ext(inputfile[1]) == "xls")){ #check the file format to call right function to read it
    rawdata <- readWorksheetFromFile(inputfile$datapath, sheet=1)
  } else if (tools::file_ext(inputfile[1]) == "csv"){
    rawdata <- read.csv(inputfile$datapath)
  }
  
  #colnames(rawdata) <- gsub("M.", "", colnames(rawdata))
  newcolnames <- unlist(strsplit(colnames(rawdata), ".", fixed = TRUE))
  colnames(rawdata) <- newcolnames[!newcolnames=="M"]
  rawdata[is.na(rawdata)]=0
  
  rawdata <- rawdata[ , mixedsort(names(rawdata))] #this sorts the rawdata for well names, requires gtools
  

  #normdata <- apply(rawdata, 2, norm.well) #w/o dependencies 
  normdata <- colwise(norm.well)(rawdata) # colwise needs library(plyr)

  ms_rawdata <- rawdata[1:(nrow(rawdata)-dc),] # isolates MS values of raw data
  ms_normdata <- normdata[1:(nrow(normdata)-dc),] # isolates MS values of normalized data
  ms_normdata[is.nan(colSums(ms_normdata))] <- 0 #first replace all colSums resulting in NaNs with zeros

  ms_rawdata$time <- seq(0, (nrow(ms_rawdata)-1)*10, 10) # adds time values in seconds
  ms_normdata$time <- seq(0, (nrow(ms_normdata)-1)*10, 10) # adds time values in seconds

  ms_rd_melted <- melt(ms_rawdata, id.vars = "time", variable.name = "well") # formatting for ggplot
  ms_nd_melted <- melt(ms_normdata, id.vars = "time", variable.name = "well") # formatting for ggplot

  ms_normdata_copy <- ms_normdata

  #ms_normdata_copy[is.nan(colSums(ms_normdata_copy))] <- 0 #first replace all colSums resulting in NaNs with zeros
  ms_normdata_copy[colSums(ms_normdata_copy)==0] <- NA #then set every 0 column to NA 

  plate_mean <- matrix(ncol = 2, nrow = nrow(ms_normdata_copy))
  plate_mean[,1] <- ms_normdata_copy$time
  plate_mean[,2] <- rowMeans(ms_normdata_copy[!ms_normdata_copy$time], na.rm=TRUE)
  plate_mean <- data.frame(plate_mean)
  colnames(plate_mean) <- c("time", "values")

  if (data_type == 1){
    chosen_set <- ms_rd_melted
  } else if (assay_type == 2){
    chosen_set <- ms_rd_melted
  }else {
    chosen_set <- ms_nd_melted
  }

  if (ncol(rawdata) > 96){
    well_plate <- 384
  } else {
    well_plate <- 96
  }

  calc_output <- list("ms_rawdata" = ms_rawdata, 
                      "ms_normdata" = ms_normdata, 
                      "ms_rd_melted" = ms_rd_melted, 
                      "ms_nd_melted" = ms_nd_melted,
                      "chosen_set" = chosen_set,
                      "yrange" = c(min(chosen_set$value), max(chosen_set$value)),
                      "xrange" = c(min(chosen_set$time), max(chosen_set$time)),
                      "well_plate" = well_plate,
                      "plate_mean" = plate_mean)
  
  return(calc_output)
}
# draw plate layout

draw_layout <- function(plate_layout, well_plate, layoutfile){

  glayout <- ggplot(data=plate_layout, aes(x=Column, y=Row)) +
    scale_alpha_discrete(range = c(0.1, 0.005)) +
    labs(title= paste("Plate Layout:", layoutfile))

  if(well_plate == 96){
    glayout <- glayout +
      geom_point(size=9) +
      geom_point(size=7, aes(colour = elicitor)) +
      geom_point(size=5, aes(shape = genotype)) +
      coord_fixed(ratio=(13/12)/(9/8), xlim=c(0.8, 12.2), ylim=c(0.6, 8.4)) +
      scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) +
      scale_x_continuous(breaks=seq(1, 12))
  } else if (well_plate == 384){
    glayout <- glayout +
      geom_point(size=6) +
      geom_point(size=4, aes(colour = genotype)) +
      geom_point(size=5, aes(shape = elicitor)) +
      coord_fixed(ratio=(25/24)/(17/16), xlim=c(0.8, 24.2), ylim=c(0.6, 16.4)) +
      scale_y_reverse(breaks=seq(1, 16), labels=LETTERS[1:16]) +
      scale_x_continuous(breaks=seq(1, 24))
  }
  glayout <- glayout + theme_bw() + guides(colour = guide_legend(title = "Elicitors", ncol = 2, byrow = TRUE), alpha = guide_legend(title = "Genotypes", ncol = 2, byrow = TRUE)) + #scale_color_npg() + 
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
}


# draw well plate

draw_well_plate <- function(data, empty_wells, mean_overlay, mean_overlay_plate, mean_overlay_WT, file_name, data_file, xlim, ylim){
  chosen_set <- data$chosen_set
  
  if(is.null(empty_wells) == FALSE){
    chosen_set$value[chosen_set$well %in% empty_wells] <- 0
  }
  
  g <- ggplot(chosen_set, aes(time, value)) 

  if(data$well_plate == 384){
    g <- g + facet_wrap(~well, ncol = 24)
  } else {
    g <- g + facet_wrap(~well, ncol = 12)
  }


  if((is.null(mean_overlay) == FALSE) && (mean_overlay)){
    mean_melted <- mean_max_calculate()$mean_melted
    g <- g + geom_line(data=mean_melted, aes(time, value), color = "red")
  }


  if((is.null(mean_overlay_plate) == FALSE) && (mean_overlay_plate)){
    g <- g + geom_line(data=data$plate_mean, aes(time, values), color = "limegreen")
  }

  if((is.null(mean_overlay_WT) == FALSE) && (mean_overlay_WT)){
    WTplate_mean <- calculate_plate_mean()
    g <- g + geom_line(data=WTplate_mean, aes(time, values), color = "blue")
  }

  if(file_name){g <- g + labs(title = paste("rawdata file:", data_file))}
  
  g <- g + geom_line() + theme_bw() + #scale_color_npg() + 
    theme(
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
                    x = (xlim[2]*0.85),
                    y = (ylim*0.9),
                   label = colnames(data$ms_normdata)[1:(length(data$ms_normdata)-1)],
                   size = 2.5)

  g <- g + ylim(0, ylim) + xlim(xlim[1], xlim[2])



  #if (data$yaxis_max != input$ylim){
  #  g <- g + ylim(0, input$ylim)
  #}
}
####################################################

shinyServer(function(input, output){
  
  calculate_plate_mean <- reactive({
    
    dc <- input$dc
    
    inputfile <- input$WT_file
    
    if (is.null(inputfile)) return(NULL)
    
    if((tools::file_ext(inputfile[1]) == "xlsx")||(tools::file_ext(inputfile[1]) == "xls")){
      rawdata <- readWorksheetFromFile(inputfile$datapath, sheet=1)
    } else if (tools::file_ext(inputfile[1]) == "csv"){
      rawdata <- read.csv(inputfile$datapath)
    }
    
    #colnames(rawdata) <- gsub("M.", "", colnames(rawdata))
    newcolnames <- unlist(strsplit(colnames(rawdata), ".", fixed = TRUE))
    colnames(rawdata) <- newcolnames[!newcolnames=="M"]
    rawdata[is.na(rawdata)]=0
    
    
    #normdata <- apply(rawdata, 2, norm.well) #w/o dependencies 
    normdata <- colwise(norm.well)(rawdata) # colwise needs library(plyr)
    
    ms_normdata <- normdata[1:(nrow(normdata)-dc),] # isolates MS values of normalized data
    
    ms_normdata$time <- seq(0, (nrow(ms_normdata)-1)*10, 10) # adds time values in seconds
    
    ms_normdata[colSums(ms_normdata)==0] <- NA
    
    plate_mean <- matrix(ncol = 2, nrow = nrow(ms_normdata))
    plate_mean[,1] <- ms_normdata$time
    plate_mean[,2] <- rowMeans(ms_normdata[!ms_normdata$time], na.rm = TRUE)
    plate_mean <- data.frame(plate_mean)
    colnames(plate_mean) <- c("time", "values")
    
    return(plate_mean)
    
  })
  
  calculate <- reactive({

    if (input$assay_type == 2){
      dc <- 0
    } else { 
      dc <- input$dc
    }
    
    inputfile <- input$data_file
    if (is.null(inputfile)) return(NULL)
    
    calc_output <- calculate_data(inputfile, dc, input$data_type, input$assay_type)
    
    return(calc_output)
    
  })
  
  multiple_calulate <- reactive({
    if(is.null(input$data2)&is.null(input$data3)&is.null(input$data4)){
      return(NULL)
    }
    
    if (input$assay_type == 2){
      dc <- 0
    } else { 
      dc <- input$dc
    }
    
    output <- list()
    data1 <- calculate()
    output$data1 <- data1
    
    if(!is.null(input$data2)){
      data2 <- calculate_data(input$data2, dc, input$data_type, input$assay_type)
      output$data2 <- data2
    }
    if(!is.null(input$data3)){
      data3 <- calculate_data(input$data3, dc, input$data_type, input$assay_type)
      output$data3 <- data3
    }
    if(!is.null(input$data4)){
      data4 <- calculate_data(input$data4, dc, input$data_type, input$assay_type)
      output$data4 <- data4
    }
      
    return(output)
  })
  

  get_layout <- reactive({
    inputlayout <- input$layout_file
    
    output <- layout(inputlayout)
    return(output)
    
  })
  
  
  plate_layout <- reactive({
    
    plate_layout <- get_layout()$plate_layout
    well_plate <- calculate()$well_plate
    
    if (is.null(plate_layout))
      return(NULL)
    if (is.null(well_plate))
      well_plate = 96
    layoutfile <- input$layoutfile
    glayout <- draw_layout(plate_layout, well_plate, layoutfile)
    
    return(glayout)
    
  })
  
  mean_max_calculate <- reactive({
    
    plate_layout <- get_layout()$plate_layout
    
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
    ms_normdata[colSums(ms_normdata)==0] <- NA #fill empty columns with NA so they dont count into the mean
    ms_normdata_names <- ms_normdata
    colnames(ms_normdata_names) <- paste(plate_layout[plate_layout$well==colnames(ms_normdata_names),]$genotype, plate_layout[plate_layout$well==colnames(ms_normdata_names),]$elicitor)
    all_means <- matrix(nrow = length(elicitors)*length(genotypes), ncol = 4, dimnames = NULL)
    all_means_graph <- matrix(nrow = length(elicitors)*length(genotypes)*length(ms_normdata[,1]), ncol = 5, dimnames = NULL)
    
    means <- matrix(nrow = nrow(data$ms_normdata), ncol = ncol(data$ms_normdata)) # last column stores time
    colnames(means) <- colnames(data$ms_normdata)[1:(ncol(data$ms_normdata))]
    means[,"time"] <- seq(0, (nrow(data$ms_rawdata)-1)*10, 10) # adds time values in seconds
    
    all_means_shaped <- matrix(nrow = 2, ncol = (length(elicitors)*length(genotypes)), dimnames = NULL)
    sd_shaped <- all_means_shaped
    graphs_shaped <- matrix(nrow = nrow(data$ms_normdata), ncol = (length(elicitors)*length(genotypes)), dimnames = NULL)
    graphs_shaped_sd <- graphs_shaped
    info <- c()
    
    i = 1
    i2 = 1
    
    for(eli in elicitors){
      for(gen in genotypes){
        wells <- raw_plate_layout$well[raw_plate_layout$elicitor == eli & raw_plate_layout$genotype == gen]
        current_set <- ms_normdata[colnames(ms_normdata) %in% wells]
        
        current_rowmeans <- rowMeans(current_set, na.rm = TRUE)
        current_rowsds <- apply(current_set, 1, sd, na.rm = TRUE)
        current_max <- max(current_rowmeans[(12+input$exclude4max):length(current_rowmeans)], na.rm = TRUE)
        current_sd <- current_rowsds[current_rowmeans==current_max]
        #print(length(c(eli, gen, current_max, current_sd)))
        #print(c(paste(gen, eli, sep = " "), current_max, current_sd))
        all_means[i,] <- c(eli, gen, current_max, current_sd)
        all_means_shaped[1,i] <- current_max
        sd_shaped[1,i] <- current_sd
        all_means_graph[i2:(i2+length(current_rowmeans)-1),1] <- rep(eli, length(current_rowmeans))
        all_means_graph[i2:(i2+length(current_rowmeans)-1),2] <- rep(gen, length(current_rowmeans))
        all_means_graph[i2:(i2+length(current_rowmeans)-1),3] <- current_rowmeans
        all_means_graph[i2:(i2+length(current_rowmeans)-1),4] <- current_rowsds
        all_means_graph[i2:(i2+length(current_rowmeans)-1),5] <- seq(0, (nrow(ms_normdata)-1)*10, 10)
        graphs_shaped[,i] <- current_rowmeans
        graphs_shaped_sd[,i] <- current_rowsds
        
        means[,colnames(means) %in% wells] <- current_rowmeans
        
        info <- append(info, paste(gen, eli, sep = " "))
        i = i+1
        i2 = i2 + length(current_rowmeans)
      }
    }
    
    df_mean <- data.frame(means)
    mean_melted <- melt(df_mean, id.vars = "time", variable.name = "well") # formatting for ggplot
    
    colnames(graphs_shaped) <- info
    colnames(graphs_shaped_sd) <- rep("SD", ncol(graphs_shaped_sd))
    colnames(all_means_shaped) <- info
    colnames(sd_shaped) <- rep("SD", ncol(sd_shaped))
    graphs_shaped_w_sd <- cbind(graphs_shaped, graphs_shaped_sd)
    graphs_shaped_w_sd <- graphs_shaped_w_sd[,c(matrix(1:ncol(graphs_shaped_w_sd), nrow = 2, byrow = TRUE))]
    all_means_shaped_w_sd <- cbind(all_means_shaped, sd_shaped)
    all_means_shaped_w_sd <- all_means_shaped_w_sd[,c(matrix(1:ncol(all_means_shaped_w_sd), nrow = 2, byrow = TRUE))]
    graphs_shaped_w_sd <- cbind(seq(0, (nrow(graphs_shaped_w_sd)-1)*10, 10), graphs_shaped_w_sd)
    colnames(graphs_shaped_w_sd)[1] <- "time"
    graphs_shaped <- cbind(seq(0, (nrow(graphs_shaped)-1)*10, 10), graphs_shaped)
    colnames(graphs_shaped)[1] <- "time"
    total_output <- list("all_means" = all_means,
                         "all_means_shaped" = all_means_shaped,
                         "all_means_shaped_w_sd" = all_means_shaped_w_sd,
                         "all_means_graph" = all_means_graph,
                         "graphs_shaped" = graphs_shaped,
                         "graphs_shaped_w_sd" = graphs_shaped_w_sd,
                         "mean_melted" = mean_melted,
                         "ms_normdata_names" = ms_normdata_names)
    
    return(total_output)
    
  })
  
  ROS_calculate <- reactive({
    
    if(input$assay_type==1)
      return(NULL)
    plate_layout <- get_layout()$plate_layout
    if (is.null(plate_layout))
      return(NULL)
    
    raw_plate_layout <- plate_layout[,1:4]
    raw_plate_layout <- raw_plate_layout[complete.cases(raw_plate_layout),] #delete rows with NAs
    
    elicitors <- unique(raw_plate_layout$elicitor)
    genotypes <- unique(raw_plate_layout$genotype)
    
    data <- calculate()
    if (is.null(data))
      return(NULL)
    
    bg_values <- 10
    how_many_bg_values <- 5
    
    #################################functions
    sdom.calculate <- function(x){
      xsdom <- sd(x)/sqrt(length(x));
      return(xsdom)
    }
    
    #calculate means and normalize to control
    calculate.sets <- function(dataset, raw_plate_layout, sdom = TRUE){
      
      eli <- unique(raw_plate_layout$elicitor)
      gen <- unique(raw_plate_layout$genotype)
      
      dataset_mean <- matrix(nrow = nrow(dataset), ncol = length(elicitors)*length(genotypes), dimnames = NULL)
      dataset_sd <- matrix(nrow = nrow(dataset), ncol = length(elicitors)*length(genotypes), dimnames = NULL)
      dataset_info <- matrix(nrow = 2, ncol = length(elicitors)*length(genotypes), dimnames = NULL)
      
      columnames <- c()
      
      i = 1
      
      for(eli in elicitors){
        for(gen in genotypes){
          wells <- raw_plate_layout$well[raw_plate_layout$elicitor == eli & raw_plate_layout$genotype == gen]
          current_set <- dataset[colnames(dataset) %in% wells]
          current_rowmeans <- rowMeans(current_set)
          if(sdom){
            current_rowsds <- apply(current_set, 1, sdom.calculate)
          } else{
            current_rowsds <- apply(current_set, 1, sd)
          }
          
          dataset_mean[,i] <- current_rowmeans
          dataset_sd[,i] <- current_rowsds
          
          dataset_info[1,i] <- gen
          
          dataset_info[2,i] <- eli
          
          
          i = i+1
          columnames <- append(columnames, paste(gen, eli, sep = " "))
          
        }
      }
      
      dataset_mean <- rbind(dataset_info, dataset_mean)
      #  dataset_sd <- rbind(dataset_info, dataset_sd)
      
      colnames(dataset_mean) <- columnames
      colnames(dataset_sd) <- columnames
      colnames(dataset_info) <- columnames
      rownames(dataset_info) <- c("genotype", "elicitor")
      
      #  normed_means <- dataset_mean #dublicate dataset_mean to use frost two rows
      normed_mean_values <- apply(dataset_mean, 2, subtract.control, dataset_mean = dataset_mean)
      
      #  normed_means[3:nrow(dataset_mean),] <- normed_mean_values
      
      output <- list("normed_means" = normed_mean_values,
                     "dataset_means" = dataset_mean,
                     "dataset_info" = dataset_info,
                     "sd" = dataset_sd)
      
      return(output)
    }
    
    #subtract background
    subtract.background <- function(x, bg_values, how_many_bg_values, with_info = FALSE){
      if (with_info){
        backgroundvalues <- as.numeric(x[(3 + bg_values - how_many_bg_values):(3 + bg_values-1)])
        mean_background <- mean(backgroundvalues)
        normalized_values <- as.numeric(x[3:length(x)]) - mean_background
      } else{
        backgroundvalues <- as.numeric(x[(bg_values - how_many_bg_values):(bg_values-1)])
        mean_background <- mean(backgroundvalues)
        normalized_values <- as.numeric(x) - mean_background
      }
      
      return(normalized_values)
    }
    
    #subtract control (needed for calculate.sets)
    subtract.control <- function(x, dataset_mean){
      values <- as.numeric(x[3:nrow(dataset_mean)])
      genotype <- x[1]
      blank <- as.numeric(dataset_mean[,dataset_mean[1,]==genotype & dataset_mean[2,]=="control"][3:nrow(dataset_mean)])
      normed_mean = values - blank
      return(normed_mean)
    }
    
    #ggplots formatting
    annotate_type <- function(x, data_info, type){
      type <- data_info[,colnames(data_info) == x["name"]][rownames(data_info) == type ]
      return(type)
    }
    
    #ggplots formatting
    add_sd <- function(x, data_info, type){
      type <- data_info[,colnames(data_info) == x["name"]][rownames(data_info) == type ]
      return(type)
    }
    
    ##################### calculations #####################
    
    plate_layout <- get_layout()$plate_layout
    
    if (is.null(plate_layout))
      return(NULL)
    
    raw_plate_layout <- plate_layout[,1:4]
    raw_plate_layout <- raw_plate_layout[complete.cases(raw_plate_layout),] #delete rows with NAs
    
    rawdata <- calculate()$ms_rawdata
    
    
    #1 subtract control from measurements
    meandata <- calculate.sets(rawdata, raw_plate_layout)
    
    #2 subtract background
    normdata <- meandata$normed_means #duplicate to overwrite
    normdata <- apply(normdata, 2, subtract.background, bg_values = bg_values, how_many_bg_values = how_many_bg_values)
    
    #3 prepare for plotting with ggplot2
    rownames(normdata) <- seq(-bg_values, (nrow(normdata)-bg_values-1), 1)
    
    normdata_melted <- melt(normdata, varnames = c("time", "name"), value.name = "values", as.is = TRUE)
    normdata_melted$genotype <- apply(normdata_melted, 1, annotate_type, data_info = meandata$dataset_info, type = "genotype")
    normdata_melted$elicitor <- apply(normdata_melted, 1, annotate_type, data_info = meandata$dataset_info, type = "elicitor")
    
    sd <- meandata$sd
    rownames(sd) <- seq(-bg_values, (nrow(sd)-bg_values-1), 1)
    sd_melted <-melt(sd, varnames = c("time", "name"), value.name = "sd", as.is = TRUE)
    normdata_melted$sd <- sd_melted$sd
    
    class(normdata_melted$time) <- "numeric"
    
    maxima <- apply(normdata[(bg_values+input$exclude4max+1):nrow(normdata),1:ncol(normdata)], 2, max) # change here for maxima settings +2 = skip first vaulue
    maxima <- maxima[!grepl("control", names(maxima))] # removes control values
    maxima_sd <- sd[normdata %in% maxima]
    
    print(normdata[normdata %in% maxima])
    
    
    maxima_melted <- melt(maxima, value.name = "values")
    maxima_melted$name <- rownames(maxima_melted)
    rownames(maxima_melted) <- NULL
    maxima_sd_melted <- melt(maxima_sd, value.name = "sd")
    
    maxima_melted$genotype <- apply(maxima_melted, 1, annotate_type, data_info = meandata$dataset_info, type = "genotype")
    maxima_melted$elicitor <- apply(maxima_melted, 1, annotate_type, data_info = meandata$dataset_info, type = "elicitor")
    
    maxima_melted$sd <- maxima_sd_melted$sd #Error here if a maximum value appears twice in the dataset 
    
    sd2 <- sd
    colnames(sd2) <- replicate(ncol(sd2), "SDOM")
    normdata_w_sd <- cbind(normdata, sd2)
    normdata_w_sd <- normdata_w_sd[,c(matrix(1:ncol(normdata_w_sd), nrow = 2, byrow = TRUE))]
    normdata_w_sd <- cbind((seq(-bg_values, (nrow(normdata_w_sd)-bg_values-1), 1)), normdata_w_sd)
    normdata <- cbind((seq(-bg_values, (nrow(normdata)-bg_values-1), 1)), normdata)
    colnames(normdata)[1] <- "time"
    colnames(normdata_w_sd)[1] <- "time"
    
    output <- list("normdata_melted" = normdata_melted,
                   "normdata" = normdata,
                   "normdata_sd" = sd,
                   "normdata_w_sd" = normdata_w_sd,
                   "maxima_melted" = maxima_melted,
                   "maxima" = maxima,
                   "maxima_sd" = maxima_sd,
                   "maxima_w_sd" = cbind(maxima, maxima_sd))
    
    return(output)
    
  })
  
  well_plots <- reactive({
    
    data <- calculate()
    
    if ((is.null(data))||(is.null(input$ylim))||(is.null(input$xlim)))
      return(NULL)
    
    empty_wells <- get_layout()$empty_wells
    
    g <- draw_well_plate(data, empty_wells, input$mean_overlay, input$mean_overlay_plate, 
                         input$mean_overlay_WT, input$file_name, input$data_file, input$xlim, input$ylim)

    return(g)

  })

  
  
  bar_plots_max <- reactive({
    
    plate_layout <- get_layout()$plate_layout
    
    if (is.null(plate_layout))
      return(NULL)
    
    raw_plate_layout <- plate_layout[,1:4]
    raw_plate_layout <- raw_plate_layout[complete.cases(raw_plate_layout),] #delete rows with NAs
    elicitors <- unique(raw_plate_layout$elicitor)
    genotypes <- unique(raw_plate_layout$genotype)
    
    if (input$assay_type == 1){
      complete_data <- mean_max_calculate()
      all_means <- complete_data$all_means
      
      header <- "Maxima with SD"
      ylabel <- "L/Lmax"
      
      num_all_means <- all_means[,3:4]
      class(num_all_means) <- "numeric" # sonst werden max und se nicht alszahlen erkannt...
      max_barplot <- data.frame(
        elicitor = all_means[,1],
        genotype = factor(all_means[,2], levels = genotypes),
        values = num_all_means[,1],
        sd = num_all_means[,2]
      )
    } else if (input$assay_type == 2){
      complete_data <- ROS_calculate()
      header <- "Maxima with SDOM"
      ylabel <- "Relative Luminescence"
      max_barplot <- complete_data$maxima_melted
      
    }
    
    if (input$bar_columns == 2){
      bpmax <- ggplot(max_barplot, aes(fill=elicitor, y=values, x=elicitor)) + facet_wrap(~genotype, ncol = 2) + ggtitle(header) + theme(legend.position = "none")
    } else {
      bpmax <- ggplot(max_barplot, aes(fill=elicitor, y=values, x=elicitor)) + facet_wrap(~genotype, ncol = 1) + ggtitle(header) + theme(legend.position = "none")
      }
    
    if(input$bar_rotation == 2){bpmax <- bpmax + coord_flip() + theme(legend.position = "none")}
    
    bpmax <- bpmax + geom_bar(position="dodge", stat="identity")
    bpmax <- bpmax + geom_errorbar(aes(ymin=values-sd, ymax=values+sd),
                                   width=.2,                    # Width of the error bars
                                   position=position_dodge(.9))
    
    
    bpmax <- bpmax + labs(x="", y=ylabel) + #scale_fill_npg() + 
      theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5),
                                                    legend.title=element_blank(),
                                                    panel.background = element_rect(fill = "white", colour = "grey90"),
                                                    panel.grid.major = element_line(color = "grey90"),
                                                    strip.background = element_rect(fill = "grey90", colour = NA),
                                                    panel.grid.minor = element_line(colour = "grey90", size = 0.25))
    
    return(bpmax)
    
  })
  
  
  mean_graphs <- reactive({
    
    plate_layout <- get_layout()$plate_layout
    
    if (is.null(plate_layout))
      return(NULL)
    
    raw_plate_layout <- plate_layout[,1:4]
    raw_plate_layout <- raw_plate_layout[complete.cases(raw_plate_layout),] #delete rows with NAs
    
    elicitors <- unique(raw_plate_layout$elicitor)
    genotypes <- unique(raw_plate_layout$genotype)
    
    if (input$assay_type == 1){
      complete_data <- mean_max_calculate()
      header <- "Mean with SD"
      ylabel <- "L/Lmax"
      all_means_graph <- complete_data$all_means_graph
    
      num_all_means_graph <- all_means_graph[,3:5]
      class(num_all_means_graph) <- "numeric" # sonst werden max und se nicht als zahlen erkannt...
      max_lineplot <- data.frame(
        elicitor = all_means_graph[,1],
        genotype = factor(all_means_graph[,2], levels = genotypes),
        values = num_all_means_graph[,1],
        sd = num_all_means_graph[,2],
        time = num_all_means_graph[,3]
      )
    } else if(input$assay_type == 2){
      complete_data <- ROS_calculate()
      header <- "Mean with SDOM"
      ylabel <- "Relative Luminescence"
      max_lineplot <- complete_data$normdata_melted
      
    }
    
    lpmax2 <- ggplot(max_lineplot, aes(x=time, y=values)) #+ 
      #scale_color_manual(values=c("navy", "red3", "deeppink1", "greenyellow", "gray20", "green4", "moccasin", "blue2", "maroon4", "turquoise4", "purple3", "pink3", "lemonchiffon2", "gray70", "tomato2", "limegreen", "salmon1", "hotpink"))
    
    if(input$graph_sorting == 2){
      lpmax2 <- lpmax2 + geom_line(aes(color = elicitor)) + facet_wrap(~genotype)
      
      if(1 %in% input$settings_mean){
        lpmax2 <- lpmax2 + ggtitle(header) +
          geom_line(aes(x=time,y=values-sd, color = elicitor), alpha=0.2) + 
          geom_line(aes(x=time,y=values+sd, color = elicitor), alpha=0.2) +
          geom_errorbar(aes(ymin=values-sd, ymax=values+sd, color = elicitor), alpha=0.2,
                        width=0,                    # Width of the error bars
                        position=position_dodge(.9))
      } else{
        lpmax2 <- lpmax2 + ggtitle("Mean")
      }
    } else {
      lpmax2 <- lpmax2 + geom_line(aes(color = genotype)) + facet_wrap(~elicitor)
      
      if(1 %in% input$settings_mean){
        lpmax2 <- lpmax2 + ggtitle(header) +
          geom_line(aes(x=time,y=values-sd, color = genotype), alpha=0.2) + 
          geom_line(aes(x=time,y=values+sd, color = genotype), alpha=0.2) +
          geom_errorbar(aes(ymin=values-sd, ymax=values+sd, color = genotype), alpha=0.2,
                        width=0,                    # Width of the error bars
                        position=position_dodge(.9))
      } else{
        lpmax2 <- lpmax2 + ggtitle("Mean")
      }
    }
    
    
    lpmax2 <- lpmax2 + labs(x="", y=ylabel) + #scale_color_npg() +
      theme(legend.title=element_blank(),
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
  
  output$multiple_wellplots <- renderPlot({
    if ((is.null(data))||(is.null(input$ylim))||(is.null(input$xlim)))
      return(NULL)
    empty_wells <- get_layout()$empty_wells
    complete_data <- multiple_calulate()
    if (is.null(complete_data))
      return(NULL)
    plots <- list()
    i=1
    empty_wells <- get_layout()$empty_wells
    for (dataset in complete_data){
      g <- draw_well_plate(dataset, empty_wells, input$mean_overlay, input$mean_overlay_plate, 
                           input$mean_overlay_WT, input$file_name, input$data_file, input$xlim, input$ylim)
      plots[[i]] <- g
      i = i+1
    }
    
    do.call("grid.arrange", c(plots, ncol=2))
    
  })
  
  output$mappingplot <- renderPlot({
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
    downloadButton('downloadData', 'Download Data')
  })
  output$downloadData <- downloadHandler(
    
    filename = function() {paste0("norm_", gsub(unlist(strsplit(input$data_file$name, "[.]")[1])[-1], "", input$data_file$name), "xlsx")},
    content = function(file) {
      data <- calculate()
      ROSdata <- NULL
      CaMean <- NULL
      if(input$assay_type==1){CaMean <- mean_max_calculate()}
      if(input$assay_type==2){ROSdata <- ROS_calculate()}
      wellplot <- well_plots()
      layout <- plate_layout()
      barplot <- bar_plots_max()
      mean_graph <- mean_graphs()
      #writeWorksheetToFile(file, data = data$ms_normdata, sheet = "normalized data")
      fname <- paste0("norm_", gsub(unlist(strsplit(input$data_file$name, "[.]")[1])[-1], "", input$data_file$name), "xlsx")
      wb <- loadWorkbook(fname, create = TRUE)
      createSheet(wb, name = "raw data")
      writeWorksheet(wb, data$ms_rawdata, sheet = "raw data")
      if(input$assay_type== 1){
        createSheet(wb, name = "normalized data")
        writeWorksheet(wb, data$ms_normdata, sheet = "normalized data")}
      if((is.null(ROSdata)==FALSE)||(is.null(CaMean)==FALSE)){
        createSheet(wb, name = "data with names")
        createSheet(wb, name = "mean data")
        createSheet(wb, name = "maxima")
        if(input$assay_type==1){
          writeWorksheet(wb, CaMean$ms_normdata_names, sheet = "data with names")
          writeWorksheet(wb, CaMean$graphs_shaped, sheet = "mean data")
          writeWorksheet(wb, CaMean$all_means_shaped, sheet = "maxima")
          createSheet(wb, name = "mean data with sd")
          createSheet(wb, name = "maxima with sd")
          writeWorksheet(wb, CaMean$graphs_shaped_w_sd, sheet = "mean data with sd")
          writeWorksheet(wb, CaMean$all_means_shaped_w_sd, sheet = "maxima with sd")
          }
        if(input$assay_type==2){
          writeWorksheet(wb, ROSdata$normdata, sheet = "mean data")
          writeWorksheet(wb, ROSdata$maxima, sheet = "maxima")
          createSheet(wb, name = "mean data with sdom")
          createSheet(wb, name = "maxima with sdom")
          writeWorksheet(wb, ROSdata$normdata_w_sd, sheet = "mean data with sdom")
          writeWorksheet(wb, ROSdata$maxima_w_sd, sheet = "maxima with sdom")
        }
        if(1 %in% input$settings_data_download1){
          png("wellplot.png", width = 1200, height = 1100)
          invisible(print(wellplot))
          dev.off()
          createSheet(wb, name = "wellplots")
          createName(wb, name = "well_plots", formula = "wellplots!$B$2")
          addImage(wb, filename = "wellplot.png", name = "well_plots", originalSize = TRUE)
          file.remove("wellplot.png")
          png("layout.png", width = 600, height = 600)
          invisible(print(layout))
          dev.off()
          createSheet(wb, name = "platelayout")
          createName(wb, name = "plate_layout", formula = "platelayout!$B$2")
          addImage(wb, filename = "layout.png", name = "plate_layout", originalSize = TRUE)
          file.remove("layout.png")
          png("barplot.png", width = 900, height = 900)
          invisible(print(barplot))
          dev.off()
          createSheet(wb, name = "maximabarplots")
          createName(wb, name = "maxima_bar_plots", formula = "maximabarplots!$B$2")
          addImage(wb, filename = "barplot.png", name = "maxima_bar_plots", originalSize = TRUE)
          file.remove("barplot.png")
          png("mean_graph.png", width = 1200, height = 1200)
          invisible(print(mean_graph))
          dev.off()
          createSheet(wb, name = "meangraphs")
          createName(wb, name = "mean_graphs", formula = "meangraphs!$B$2")
          addImage(wb, filename = "mean_graph.png", name = "mean_graphs", originalSize = TRUE)
          file.remove("mean_graph.png")
        }
      }
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
  

  output$downloadCombined <- downloadHandler(
    filename <- function(){
      paste("combined_quadrants", "csv", sep = ".")
    },
    content = function(file) {
      
      Q1 <- curateQ(input$Q1E_file$datapath, input$Q1D_file$datapath)
      Q2 <- curateQ(input$Q2E_file$datapath, input$Q2D_file$datapath)
      Q3 <- curateQ(input$Q3E_file$datapath, input$Q3D_file$datapath)
      Q4 <- curateQ(input$Q4E_file$datapath, input$Q4D_file$datapath)
      
      measurement <- Q1 + Q2 + Q3 + Q4
      
      write.csv(measurement, file = file, row.names = FALSE)
      
      
      
    }
  )
  
  output$downloadLicense <- downloadHandler(
    filename <- function(){
      paste("license", "txt", sep = ".")
    },
    content <-function(file){
      file.copy("license.txt", file)
    }
  )
  
  output$downloadCaExample <- downloadHandler(
    filename <- function(){
      paste("CaExample", "xlsx", sep = ".")
    },
    content <-function(file){
      file.copy("CaExample.xlsx", file)
    }
  )
  
  output$downloadCaExample_layout <- downloadHandler(
    filename <- function(){
      paste("CaExample_layout", "xlsx", sep = ".")
    },
    content <-function(file){
      file.copy("CaExample_layout.xlsx", file)
    }
  )
  
  output$downloadROSExample <- downloadHandler(
    filename <- function(){
      paste("ROSExample", "xlsx", sep = ".")
    },
    content <-function(file){
      file.copy("ROSExample.xlsx", file)
    }
  )
  
  output$downloadROSExample_layout <- downloadHandler(
    filename <- function(){
      paste("ROSExample_layout", "xlsx", sep = ".")
    },
    content <-function(file){
      file.copy("ROSExample_layout.xlsx", file)
    }
  )
  
  output$downloadDatatemplate <- downloadHandler(
    filename <- function(){
      paste("Datatemplate", "xlsx", sep = ".")
    },
    content <-function(file){
      file.copy("Datatemplate.xlsx", file)
    }
  )
  
  output$downloadLayouttemplate <- downloadHandler(
    filename <- function(){
      paste("Layouttemplate", "xlsx", sep = ".")
    },
    content <-function(file){
      file.copy("Layouttemplate.xlsx", file)
    }
  )
  
  output$ui.settings1 <- renderUI({
    if (input$assay_type == 2){
      HTML("<div><span style='color:red'>Keep in mind that the data has to be in the 'Luminoscan' format! </br> (for now, soon to be updated!)</span></div>")
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
  output$ui.settings3 <- renderUI({
    yrange <- calculate()$yrange
    ylim <- yrange[2]
    ymin <- yrange[1]
    if (is.null(ylim)){
      return(NULL)
    }
    ylim.steps = round(ylim/20, digits = 2)
    
#    numericInput("ylim", "y-axis limit", value = ylim, min = ymin, max = 2*ylim, step = ylim.steps)
    sliderInput("ylim", "y-axis limit", min = floor(ymin), max = ceiling(ylim*2), value = ylim, step = ylim.steps)
  })
  
  output$ui.settings4 <- renderUI({
    xrange <- calculate()$xrange
    xlim <- xrange[2]
    xmin <- xrange[1]
    if (is.null(xlim)){
      return(NULL)
    }
    sliderInput("xlim", "x-axis limit", min = xmin, max = xlim, value = c(xmin, xlim))
  })
  
  output$ui.settings5 <- renderUI({
    plate_layout <- get_layout()$plate_layout
    if ((is.null(plate_layout))||(input$assay_type == 2))
      return(NULL)
    checkboxInput("mean_overlay", label = "plot respective means", value = FALSE)
  })
  
  output$menu1 <- renderMenu({
    inputlayout <- input$layout_file
    
    if (is.null(inputlayout)==FALSE){
      menuItem("Data Summary", tabName = "data_summary", icon = icon("area-chart"),
              menuSubItem("Mean Maxima", tabName = "norm_data1", icon = icon("bar-chart")),
              menuSubItem("Mean Kinetics", tabName = "norm_data2", icon = icon("line-chart"))
      )
    } 
  })
  
  output$menu2 <- renderMenu({
    inputfile <- input$data_file
    if (is.null(inputfile)==FALSE){
      menuItem("Download", tabName = "download", icon = icon("cloud-download"))
    }
  })
    
  output$ui.plate_layout<-renderUI({
    if (is.null(get_layout()))
      return(NULL)
    
    box(title = "Plate Layout",
        solidHeader = TRUE,
        status = "success",
        width = 12,
        collapsible = TRUE,
        plotOutput("plate_layout", height = 300))
  })
  
})
