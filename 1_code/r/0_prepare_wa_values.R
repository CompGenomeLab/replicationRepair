#### library ####

library(reshape2)
library(ggplot2)
library(stringr)

#### set variables ####

dateout <- Sys.Date()
dateout <- format(dateout, format = "[%Y.%m.%d]")

setwd(paste("~/Documents/My_Projects/Project_Repair_Replication/Data/", 
            "S1toS2_regions/Original_files", sep = ""))

temp <- list.files(pattern = "PctSignalRep1.bed$")

#### import data and rearrange ####

all_phases <- data.frame()

for(counter in 1:length(temp)) {
  
  repliseq <- read.table(temp[[counter]], header = FALSE)
  
  phase <- substr(temp[[counter]], start = 25, stop = 26) 
  
  repliseq$phase <- phase
  
  all_phases <- rbind(all_phases, repliseq)
  
}

all_phases_dcasted = dcast(all_phases, V1 + V2 + V3 ~ phase, value.var = "V4")

all_phases_dcasted[is.na(all_phases_dcasted)] <- 0

# wa value calculation
all_phases_dcasted$WA = (0.917*all_phases_dcasted$G1) + 
  (0.750*all_phases_dcasted$S1) + (0.583*all_phases_dcasted$S2) + 
  (0.417*all_phases_dcasted$S3) + (0.250*all_phases_dcasted$S4) + 
  (0*all_phases_dcasted$G2)

all_phases_dcasted = all_phases_dcasted[ , c("V1", "V2", "V3", "G1", "S1", 
                                             "S2", "S3", "S4", "G2", "WA")]

# adding quartiles to wa values
all_phases_dcasted$quartiles <- NA

all_phases_dcasted <- within(all_phases_dcasted, quartiles
                             [81 <= WA & WA <= 100] <- "81-100%")

all_phases_dcasted <- within(all_phases_dcasted, quartiles
                             [61 <= WA & WA < 81] <- "61-80%")

all_phases_dcasted <- within(all_phases_dcasted, quartiles
                             [41 <= WA & WA < 61] <- "41-60%")

all_phases_dcasted <- within(all_phases_dcasted, quartiles
                             [21 <= WA & WA < 41] <- "21-40%")

all_phases_dcasted <- within(all_phases_dcasted, quartiles
                             [0 <= WA & WA < 12] <- "0-20%")

# write to a file
setwd(paste("~/Documents/My_Projects/Project_Repair_Replication/Data/", 
            "S1toS2_regions/WA_scores", sep = ""))

write.table(all_phases_dcasted, file = "WA_values.txt", sep = "\t", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)


# rearrange as bed file
all_phases_dcasted$name <- "weighted_average" 

all_phases_dcasted$empty <- "." 

all_phases_dcasted_bed <- all_phases_dcasted[ , c("V1", "V2", "V3", "name", 
                                                  "WA", "empty")]

# write to a file
write.table(all_phases_dcasted_bed, file = "WA_values.bed", sep = "\t", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)

#### plot histogram ####

all_phases_dcasted <- read.table("WA_values.txt", header = FALSE, sep = "\t")

ggplot(all_phases_dcasted) + geom_bar(aes(x=V10), colour="red") 

setwd(paste("~/Documents/My_Projects/Project_Repair_Replication/Results/", 
            "S1toS2_regions/", sep = ""))

ggsave(filename = paste(dateout, "histogram_WA.png", sep = ""), 
       width = 297, height = 210, units = "mm") 

rm(all_phases, all_phases_dcasted, all_phases_dcasted_bed, repliseq, counter, 
   dateout, phase, temp)


