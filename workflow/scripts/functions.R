#### libraries ####

library(stringr)

#### functions ####

rmv_low_counted_regions <- function ( df, count_threshold = 1 ){
  trash = filter(df, counts < count_threshold)
  df = df[!(df$chromosomes %in% trash$chromosomes & 
            df$start_position %in% trash$start_position &
            df$end_position %in% trash$end_position), ]
  
  return(df)
}

repair_rate <- function ( df, ds_rep = "A" ){
  ds <- filter(df, method == "Damage_seq" & replicate == ds_rep)
  xr <- filter(df, method == "XR_seq")
  xr_list <- split(xr, xr$replicate)
  df_xr_ds <- xr_list[[1]][0, ]
  
  for ( xr_rep in 1:length(xr_list) ){
    temp <- xr_list[[xr_rep]]
    temp <- rbind(temp, ds) 
    temp <- dcast(temp, chromosomes + start_position + end_position + 
                    dataset + score + dataset_strand + product + phase + 
                    time_after_exposure + sample_strand ~ method, 
                  value.var = "RPKM")
    temp <- cbind(temp, xr_list[[xr_rep]]["replicate"])
    df_xr_ds <- rbind(df_xr_ds, temp)
  }
  
  df_xr_ds$xr_ds <- df_xr_ds$XR_seq / df_xr_ds$Damage_seq
  df_xr_ds <- select(df_xr_ds, -c("Damage_seq", "XR_seq"))
  
  return(df_xr_ds)
}

rrEarly_rrLate <- function ( df_rr ){
  df_rr_filt <- filter(df_rr, phase != "async")
  df_rr_ear_la <- dcast(df_rr_filt, chromosomes + start_position + 
                        end_position + dataset + score + dataset_strand + 
                        product + time_after_exposure + replicate + 
                        sample_strand ~ phase, value.var = "xr_ds")
  df_rr_ear_la$ear_la <- df_rr_ear_la$early / df_rr_ear_la$late
  df_rr_ear_la <- select(df_rr_ear_la, -c(early, late))

  return(df_rr_ear_la)
}

chrState_naming <- function (df, chr_states, general_states, 
                             chrState2generalState){
  df$dataset <- factor(df$dataset, levels = chr_states)
  df$states <- NA
  
  for (i in 1:length(chr_states)) {
    df <- within(df, states[dataset == chr_states[i]] <- 
                   chrState2generalState[i])
  }
  
  df$states <- as.factor(df$states)
  df$states <- factor(
    df$states, levels = general_states)
  
  return(df)
}

rr_boxplot <- function (df){
  mut_plus <- filter(df, windows > 0)
  mut_plus_agg <- aggregate(x = mut_plus$xr_ds, by = 
                              list(mut_plus$dataset, mut_plus$repdomains, 
                                   mut_plus$sample_strand), FUN = "mean")
  mut_plus_agg$direction <- "Right Replicating"
  mut_minus <- filter(df, windows < 0)
  mut_minus_agg <- aggregate(x = mut_minus$xr_ds, by = 
                               list(mut_minus$dataset, mut_minus$repdomains, 
                                    mut_minus$sample_strand), FUN = "mean")
  mut_minus_agg$direction <- "Left Replicating"
  mut_agg <- rbind(mut_plus_agg, mut_minus_agg)
  
  return(mut_agg)
}

rr_boxplot_plus_minus <- function (df){
  mut_casted <- dcast(df, dataset + repdomains + windows ~ sample_strand, 
                      value.var = "xr_ds")
  mut_casted$plus_min <- mut_casted$"+" - mut_casted$"-" 
  mut_plus <- filter(mut_casted, windows > 0)
  mut_plus_agg <- aggregate(x = mut_plus$plus_min, by = 
                              list(mut_plus$dataset, mut_plus$repdomains), 
                            FUN = "mean")
  mut_plus_agg$direction <- "Right Replicating"
  mut_minus <- filter(mut_casted, windows < 0)
  mut_minus_agg <- aggregate(x = mut_minus$plus_min, by = 
                               list(mut_minus$dataset, mut_minus$repdomains), 
                             FUN = "mean")
  mut_minus_agg$direction <- "Left Replicating"
  mut_agg <- rbind(mut_plus_agg, mut_minus_agg)
  
  return(mut_agg)
}

output_date <- function (){
  dateout <- Sys.Date()
  dateout <- format(dateout, format = "[%Y.%m.%d]")
  return(dateout)
}


get_sub_str <- function ( whole_str, substr1, substr2 ){
  target_str <- str_match(whole_str, paste(substr1, "(.*?)", substr2, sep = ""))
  return(target_str[2])
}


middle <- function( window_number ){
  if (window_number %% 2 == 0) {
    middle <- window_number / 2
  } else {
    middle <- (window_number - 1) / 2 + 1
  }
  return(middle)
}

region_length <- function( window_number, window_length ){
  rlength <- (window_number - 1) * window_length / 2000 # in kb and half of it. 
  return(rlength)
}

window_numbering <- function( mydata, name_column, middle ){
  window <- cbind(mydata[ , name_column], 
                  data.frame(str_split_fixed(mydata[ , name_column], "_", -1)))
  window <- window[ , c(1, ncol(window))]
  names(window) <- c("dataset", "windows")
  window$windows <- as.numeric(as.character(window$windows)) - middle
  newdata <- merge(mydata, window, by.x="dataset", by.y="dataset")
  newdata <- unique(newdata) # farklı samplelar ve 2 strand olduğu için gerekli
  return(newdata)
}

domain_name <- function( mydata, name_column ){
  window <- cbind(mydata[ , name_column], 
                  data.frame(str_split_fixed(mydata[ , name_column], "_", -1)))
  window <- window[ , c(1, ncol(window)-1)]
  names(window) <- c("dataset", "repdomains")
  mydata$repdomains <- window[,2]
  mydata$repdomains <- factor(mydata$repdomains, 
                              levels = c("UTZ", "ERD", "DTZ", "LRD"))
  return(mydata)
}

fig_save <- function( figpath, figname ) {
  ggsave(path = figurePath, filename = figureName, 
         width = 297, height = 210, units = "mm")
  figurePNG <- sub(".pdf", ".png", figureName)
  ggsave(path = figurePath, filename = figurePNG, 
         width = 297, height = 210, units = "mm")
}

