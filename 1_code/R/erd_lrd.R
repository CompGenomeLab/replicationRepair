#### library ####

library(ggplot2)
library(stringr)

#### prepare data ####

sourcePath <- paste("~/Documents/myprojects/replicationRepair/1_code/R/")

date <- "[2020.02.01]"

fr_name <- "InZones_LRD_windows"

source(paste(sourcePath, "2_report_sub_dfs.R", sep = ""))

fr_LRD <- fr

fr_name <- "InZones_ERD_windows"

source(paste(sourcePath, "2_report_sub_dfs.R", sep = ""))

fr_ERD <- fr

#### set variables and rearrange ####

for (rearrange in 1) {  
  
  #### set variables ####
  
  dateout <- Sys.Date()
  dateout <- format(dateout, format = "[%Y.%m.%d]")
  
  window_number <- 201
  if (window_number %% 2 == 0) {
    half_window <- window_number / 2
  } else {
    half_window <- (window_number - 1) / 2 + 1
  }
  
  #### rearrange ####
  
  # window numbers separated from dataset names
  window <- data.frame(str_split_fixed(fr_ERD$dataset, "_", -1))
  
  fr_ERD$windows <- as.numeric(levels(
    window[ , ncol(window)]))[window[ , ncol(window)]] - half_window
  
  fr_ERD$dataset <- "ERD"
  
  rm(window, rearrange)
  
  # window numbers separated from dataset names
  window <- data.frame(str_split_fixed(fr_LRD$dataset, "_", -1))
  
  fr_LRD$windows <- as.numeric(levels(
    window[ , ncol(window)]))[window[ , ncol(window)]] - half_window
  
  fr_LRD$dataset <- "LRD"
  
  rm(window, rearrange)
  
}

df <- rbind(fr_ERD, fr_LRD)

df <- dcast(df, chromosomes + start_position + end_position + 
                score + dataset_strand + product + phase + 
                time_after_exposure + sample_strand + method ~ dataset, 
              value.var = "RPKM")
