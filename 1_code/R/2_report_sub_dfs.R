#### library #####

library(dplyr)
library(reshape2)

#### variable set ####

dataInfoPath <- paste("~/Documents/myprojects/replicationRepair/4_output/",
                  "gitignore/1_TextforPlotting/", sep = "")

temp <- list.files(path = dataInfoPath, pattern = ".txt")

temp <- temp[grepl("*final_report.*\\.txt", temp)]

# date and fr_name will be delivered by the main plot scripts (3_*)

for (run in 1) {
  #### import data ####
  
  fr <- read.csv(paste(dataInfoPath, date, "final_report_", fr_name, 
                       "_ready.csv", sep = ""))
  
  #### filter zeros ####
  
  if (grepl("windows", fr_name)) { } else {
    
    trash = filter(fr, counts < 1)
    fr = fr[!(fr$chromosomes %in% trash$chromosomes & 
            fr$start_position %in% trash$start_position &
            fr$end_position %in% trash$end_position), ]
    
    rm(trash)
  }
  
  #### xr / ds ####
  
  ds_a <- filter(fr, method == "Damage_seq" & replicate == "A")
  
  fr_xr <- filter(fr, method == "XR_seq")
  
  fr_list <- split(fr_xr, fr_xr$replicate)
  
  fr_xr_ds <- fr_list[[1]][0, ]
  
  for (xr in 1:length(fr_list)) {
    
    temp <- fr_list[[xr]]
    
    temp <- rbind(temp, ds_a)
    
    temp <- dcast(temp, chromosomes + start_position + end_position + 
                    dataset + score + dataset_strand + product + phase + 
                    time_after_exposure + sample_strand ~ method, 
                  value.var = "RPKM")
    
    temp <- cbind(temp, fr_list[[xr]]["replicate"])
  
    fr_xr_ds <- rbind(fr_xr_ds, temp)
  }
  
  fr_xr_ds$xr_ds <- fr_xr_ds$XR_seq / fr_xr_ds$Damage_seq
  
  fr_xr_ds = select(fr_xr_ds, -c("Damage_seq", "XR_seq"))
  
  rm(ds_a, fr_xr, fr_list, xr, temp)
  
  
  #### early / late after xr / ds ####
  
  fr_xr_ds_ear_la <- filter(fr_xr_ds, phase != "async")
  
  fr_xr_ds_ear_la <- dcast(fr_xr_ds_ear_la, chromosomes + start_position + 
                             end_position + dataset + score + dataset_strand + 
                             product + time_after_exposure + replicate + 
                             sample_strand ~ phase, value.var = "xr_ds")
  
  fr_xr_ds_ear_la$ear_la <- fr_xr_ds_ear_la$early / fr_xr_ds_ear_la$late
  
  fr_xr_ds_ear_la = select(fr_xr_ds_ear_la, -c(early, late))
  
  
  #### (+ / -) after xr / ds ####
  
  fr_xr_ds_plus_min <- dcast(fr_xr_ds, chromosomes + start_position + 
                             end_position + dataset + score + dataset_strand + 
                             product + phase + time_after_exposure + replicate ~ 
                             sample_strand, value.var = "xr_ds")
  
  fr_xr_ds_plus_min$plus_min <- fr_xr_ds_plus_min$`+` / fr_xr_ds_plus_min$`-`
  
  fr_xr_ds_plus_min = select(fr_xr_ds_plus_min, -c(`-`, `+`))
  
  
  #### (+ / -) after early / late after xr / ds ####
  
  fr_xr_ds_ear_la_plus_min <- dcast(fr_xr_ds_ear_la, chromosomes + 
                                      start_position + end_position + dataset + 
                                      score + dataset_strand + product + 
                                      time_after_exposure + replicate ~ 
                                      sample_strand, value.var = "ear_la")
  
  fr_xr_ds_ear_la_plus_min$plus_min <- fr_xr_ds_ear_la_plus_min$`+` / 
    fr_xr_ds_ear_la_plus_min$`-`
  
  fr_xr_ds_ear_la_plus_min = select(fr_xr_ds_ear_la_plus_min, -c(`-`, `+`))

    
  #### early / late ####
  
  fr_ear_la <- filter(fr, phase != "async")
  
  fr_ear_la <- dcast(fr_ear_la, chromosomes + start_position + 
                       end_position + dataset + score + dataset_strand + 
                       method + product + time_after_exposure + replicate + 
                       sample_strand ~ phase, value.var = "RPKM")
  
  fr_ear_la$ear_la <- fr_ear_la$early / fr_ear_la$late
  
  fr_ear_la = select(fr_ear_la, -c(early, late))
  
  #### (+ / -) after early / late ####
  
  fr_ear_la_plus_min <- dcast(fr_ear_la, chromosomes + start_position + 
                       end_position + dataset + score + dataset_strand + 
                       method + product + time_after_exposure + replicate ~ 
                       sample_strand, value.var = "ear_la")
  
  fr_ear_la_plus_min$plus_min <- fr_ear_la_plus_min$`+` / fr_ear_la_plus_min$`-`
  
  fr_ear_la_plus_min = select(fr_ear_la_plus_min, -c(`+`, `-`))
  
  
  #### xr / input ####
  
  dna_a <- filter(fr, method == "DNA_seq" & phase != "async")
  
  fr_xr = filter(fr, method == "XR_seq" & phase != "async")
  
  fr_list <- split(fr_xr, list(fr_xr$product, fr_xr$time_after_exposure, 
                               fr_xr$replicate), drop = TRUE)
  fr_xr_dna <- fr_list[[1]][0,]
  
  for (xr in 1:length(fr_list)) {
    
    temp <- fr_list[[xr]]
    
    temp <- rbind(temp, dna_a)
    
    temp <- dcast(temp, chromosomes + start_position + 
                     end_position + dataset + score + dataset_strand + 
                     phase +  sample_strand ~ method, 
                   value.var = "RPKM")
    
    temp <- cbind(temp, fr_list[[xr]]["replicate"])
    
    temp <- cbind(temp, fr_list[[xr]]["product"])
    
    temp <- cbind(temp, fr_list[[xr]]["time_after_exposure"])
    
    fr_xr_dna <- rbind(fr_xr_dna, temp)
  }
  
  fr_xr_dna$xr_dna <- fr_xr_dna$XR_seq / fr_xr_dna$DNA_seq
  
  fr_xr_dna = select(fr_xr_dna, -c("DNA_seq", "XR_seq"))
  
  rm(dna_a, fr_xr, fr_list, xr, temp)
  
  
  #### ds / input ####
  
  dna_a <- filter(fr, method == "DNA_seq" & phase != "async")
  
  fr_ds = filter(fr, method == "Damage_seq" & phase != "async")
  
  fr_list <- split(fr_ds, list(fr_ds$product, fr_ds$time_after_exposure, 
                               fr_ds$replicate), drop = TRUE)
  fr_ds_dna <- fr_list[[1]][0,]
  
  for (ds in 1:length(fr_list)) {
    
    temp <- fr_list[[ds]]
    
    temp <- rbind(temp, dna_a)
    
    temp <- dcast(temp, chromosomes + start_position + 
                    end_position + dataset + score + dataset_strand + 
                    phase +  sample_strand ~ method, 
                  value.var = "RPKM")
    
    temp <- cbind(temp, fr_list[[ds]]["replicate"])
    
    temp <- cbind(temp, fr_list[[ds]]["product"])
    
    temp <- cbind(temp, fr_list[[ds]]["time_after_exposure"])
    
    fr_ds_dna <- rbind(fr_ds_dna, temp)
  }
  
  fr_ds_dna$ds_dna <- fr_ds_dna$Damage_seq / fr_ds_dna$DNA_seq
  
  fr_ds_dna = select(fr_ds_dna, -c("DNA_seq", "Damage_seq"))
  
  rm(dna_a, fr_ds, fr_list, ds, temp)
  
  
  #### (+ / -) ####
  
  fr_plus_min <- dcast(fr, chromosomes + start_position + end_position + 
                         dataset + score + dataset_strand + method + phase + 
                         product + time_after_exposure + replicate ~ 
                         sample_strand, value.var = "RPKM")
  
  fr_plus_min$plus_min <- fr_plus_min$`+` / fr_plus_min$`-`
  
  fr_plus_min = select(fr_plus_min, -c(`-`, `+`))

  rm(run)  
}

