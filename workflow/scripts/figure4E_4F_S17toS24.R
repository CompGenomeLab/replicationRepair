#### Packages and Libraries ####

library(argparser)
library(stringr)
library(ggpubr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(patchwork)
library(GGally)
set.seed(1) 

######## Arguments ##########
p <- arg_parser("producing figure 4E, 4F, S17-S24 - the sequence context of (mono to penta)-mer around initiation zones")
p <- add_argument(p, "-i", help="input")
p <- add_argument(p, "--data_prefix", help="name prefix of the dataframes that generate the plots")
p <- add_argument(p, "--o1", help="output 1-mer")
p <- add_argument(p, "--o2", help="output 2-mer")
p <- add_argument(p, "--o3", help="output 3-mer")
p <- add_argument(p, "--o4", help="output 4-mer")
p <- add_argument(p, "--o5", help="output 5-mer")
p <- add_argument(p, "--oT", help="output polyT tracks")
p <- add_argument(p, "--oTT", help="output polyTT tracks")
p <- add_argument(p, "--oC", help="output polyC tracks")
p <- add_argument(p, "--oCC", help="output polyCC tracks")
p <- add_argument(p, "--oTC", help="output polyTC tracks")

# Parse the command line arguments
argv <- parse_args(p)

counts_txt <- argv$i

#### Default Plot Format ####

source("workflow/scripts/plot_format.R")

#### Functions ####

source("workflow/scripts/functions.R")

calculateAsymmetry <- function(df, pattern1, pattern2){
  
  for (nuc_id in 1:length(df$variable)){
    if (str_count(df$variable[nuc_id], paste0("(?=", pattern2,")")) 
        > str_count(df$variable[nuc_id], paste0("(?=", pattern1,")")) ){ 
      
      temp_base <- df$variable[nuc_id]
      df$variable[nuc_id] <- df$complementary[nuc_id]
      df$complementary[nuc_id] <- as.character(temp_base)  
    }
  }
  
  df_dcast <- dcast(df, direction + repdomains + variable + complementary ~ strand, 
                    value.var = "value")
  
  df_dcast$diff <- abs(df_dcast$`+` - df_dcast$`-`)

  df_dcast$lag_lead = 0
  df_dcast$lag_lead[df_dcast$direction == "Left Replicating"] <- df_dcast$`-`[df_dcast$direction == "Left Replicating"] / df_dcast$`+`[df_dcast$direction == "Left Replicating"]
  df_dcast$lag_lead[df_dcast$direction == "Right Replicating"] <- df_dcast$`+`[df_dcast$direction == "Right Replicating"] / df_dcast$`-`[df_dcast$direction == "Right Replicating"]
  
  baseA <- 29.474993117004857
  baseC <- 20.500734650602116
  baseG <- 20.515808882191084
  baseT <- 29.508463350201946
  
  occurT <- baseA / baseT 
  occurC <- baseG / baseC
  occurA <- baseT / baseA 
  occurG <- baseC / baseG
  
  for (i in 1:length(df_dcast$variable)){
    
    count1 <- str_count(df_dcast$variable[i], paste0("(?=", pattern1,")"))
    count2 <- str_count(df_dcast$variable[i], paste0("(?=", pattern2,")"))
    
    df_dcast$count1[i] <- count1
    df_dcast$count2[i] <- count2
    df_dcast$diff[i] <- count1 - count2
    df_dcast$sum[i] <- count1 + count2
    
    if (nchar(as.character(df_dcast$variable[i])) > 1){ 
      
      prev_nuc <- str_sub(df_dcast$variable[i], 1, 
                          nchar(as.character(df_dcast$variable[i]))-1)
      #print(prev_nuc)
      rep <- df_dcast$repdomains[i]
      direc <- df_dcast$direction[i]
      #print(rep)
      prev_obs <- df_dcast$lag_lead[
        df_dcast$variable == prev_nuc & 
          df_dcast$repdomains == rep &
          df_dcast$direction == direc] 
      if (length(prev_obs) == 0) {
        prev_obs <- df_dcast$lag_lead[
          df_dcast$complementary == prev_nuc & 
            df_dcast$repdomains == rep &
            df_dcast$direction == direc] 
      }
      #print(prev_obs)
      
      if (str_sub(df_dcast$variable[i], 
                  nchar(as.character(df_dcast$variable[i])), -1) == "T") {
        
        df_dcast$expected[i] <- prev_obs * occurT
        
      } else if (str_sub(df_dcast$variable[i], 
                         nchar(as.character(df_dcast$variable[i])), -1) == "C") {
        
        df_dcast$expected[i] <- prev_obs * occurC
        
      } else if (str_sub(df_dcast$variable[i], 
                         nchar(as.character(df_dcast$variable[i])), -1) == "G") {
        
        df_dcast$expected[i] <- prev_obs * occurG
        
      } else if (str_sub(df_dcast$variable[i], 
                         nchar(as.character(df_dcast$variable[i])), -1) == "A") {
        
        df_dcast$expected[i] <- prev_obs * occurA
      }
      
    } else if (df_dcast$variable[i] == "T") { 
      df_dcast$expected[i] <- occurT
      
    } else if (df_dcast$variable[i] == "C") { 
      df_dcast$expected[i] <- occurC 
      
    } else if (df_dcast$variable[i] == "A") { 
      df_dcast$expected[i] <- occurA
      
    } else if (df_dcast$variable[i] == "G") { 
      df_dcast$expected[i] <- occurG }
    
  }
  
  df_dcast_arranged <- df_dcast %>% 
    group_by(diff) %>% 
    arrange(lag_lead)
  
  
  return(df_dcast_arranged)
}

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

#### Main ####

counts_df <- read.table(counts_txt, header = T)

windows <- data.frame(str_split_fixed(counts_df$id, "_", -1))
counts_df$windows <- as.numeric(
  as.character(windows[,3])) - 101
counts_df$repdomains <- windows[,2]
counts_df$strand <- windows[,4]
counts_df$strand <- gsub("\\(.*", "", counts_df$strand)

counts_total <- filter(counts_df, windows < -20 | windows > 20)

counts_total$direction[ counts_total$windows > 0 ] <- "Right Replicating"
counts_total$direction[ counts_total$windows < 0 ] <- "Left Replicating"
counts_total_molten <- melt(counts_total, id = c("direction", "windows", "repdomains", "strand", "id"))

counts_total_left <- filter(counts_total_molten, direction == "Left Replicating")
counts_total_right <- filter(counts_total_molten, direction == "Right Replicating")
counts_total_left$value <- counts_total_left$value*-1
counts_total_df <- rbind(counts_total_right, counts_total_left)

counts_total_df <- data_summary(counts_total_df, varname="value",
                            groupnames=c("direction", "repdomains", "strand", "variable"))

complementary_list <- rev(unique(as.array( counts_total_df$variable )))

for (nuc_id in 1:length(complementary_list)){
  comp_nuc <- ""
  for (idx in 1:nchar(as.character(complementary_list[nuc_id]))){ 
    if (substr(complementary_list[nuc_id], idx, idx) == "A"){
      base <- "T"
    } else if (substr(complementary_list[nuc_id], idx, idx) == "T"){
      base <- "A"
    } else if (substr(complementary_list[nuc_id], idx, idx) == "G"){
      base <- "C"
    } else if (substr(complementary_list[nuc_id], idx, idx) == "C") {
      base <- "G"
    }
    comp_nuc <- paste0(base, comp_nuc)
  }
  counts_total_df$complementary[ counts_total_df$variable == complementary_list[nuc_id] ] <- comp_nuc
  if (comp_nuc == complementary_list[nuc_id]) {
    next
  } else if (complementary_list[nuc_id] %in% unique(as.array( counts_total_df$variable ))){
    counts_total_df <- filter(counts_total_df, variable != comp_nuc)
  }
}

##### Asymmetry Plots #####

counts_total_T <- calculateAsymmetry(counts_total_df, "T", "A")

write.table(counts_total_T, file = paste0(argv$data_prefix, "4_E.csv"), 
            quote = FALSE, row.names = FALSE, sep = ",")

p1 <- ggplot(data=counts_total_T, aes(x=diff, y=lag_lead)) +
  facet_grid(~repdomains) +
  geom_point(position = position_dodge2(.9), size = .5) +
  scale_x_continuous(breaks = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)) +
  xlab("# of net Ts (T-A)") +
  ylab("Lagging/Leading strand\n occurrence of sequences")

p1 <- p1 + p_format +
  geom_stripped_cols(odd = "brown", even = "white", alpha = .1) +
  theme(panel.border=element_rect(size=1, fill = NA)) +
  stat_summary(fun = median,
               geom="crossbar",
               size = 0.25,
               width = 0.5,
               color="red") 

ggsave(argv$oT, width = 11, height = 9, units = "cm") 

counts_total_C <- calculateAsymmetry(counts_total_df, "C", "G")

write.table(counts_total_C, file = paste0(argv$data_prefix, "S23.csv"), 
            quote = FALSE, row.names = FALSE, sep = ",")

p2 <- ggplot(data=counts_total_C, aes(x=diff, y=lag_lead)) +
  facet_grid(~repdomains) +
  geom_point(position = position_dodge2(.9), size = .5) +
  scale_x_continuous(breaks = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)) +
  xlab("# of net Cs (C-G)") +
  ylab("Lagging/Leading strand\n occurrence of sequences")

p2 <- p2 + p_format +
  geom_stripped_cols(odd = "brown", even = "white", alpha = .1) +
  theme(panel.border=element_rect(size=1, fill = NA)) +
  stat_summary(fun = median,
               geom="crossbar",
               size = 0.25,
               width = 0.5,
               color="red") 

ggsave(argv$oC, width = 11, height = 9, units = "cm") 

counts_total_TT <- calculateAsymmetry(counts_total_df, "TT", "AA")


write.table(counts_total_TT, file = paste0(argv$data_prefix, "4_F.csv"), 
            quote = FALSE, row.names = FALSE, sep = ",")

p3 <- ggplot(data=counts_total_TT, aes(x=diff, y=lag_lead)) +
  facet_grid(~repdomains) +
  geom_point(position = position_dodge2(.9), size = .5) +
  scale_x_continuous(breaks = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)) +
  xlab("# of net TTs (TT-AA)") +
  ylab("Lagging/Leading strand\n occurrence of sequences")

p3 <- p3 + p_format +
  geom_stripped_cols(odd = "brown", even = "white", alpha = .1) +
  theme(panel.border=element_rect(size=1, fill = NA)) +
  stat_summary(fun = median,
               geom="crossbar",
               size = 0.25,
               width = 0.5,
               color="red") 

ggsave(argv$oTT, width = 11, height = 9, units = "cm") 

counts_total_CC <- calculateAsymmetry(counts_total_df, "CC", "GG")

write.table(counts_total_CC, file = paste0(argv$data_prefix, "S24.csv"), 
            quote = FALSE, row.names = FALSE, sep = ",")

p4 <- ggplot(data=counts_total_CC, aes(x=diff, y=lag_lead)) +
  facet_grid(~repdomains) +
  geom_point(position = position_dodge2(.9), size = .5) +
  scale_x_continuous(breaks = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)) +
  xlab("# of net CCs (CC-GG)") +
  ylab("Lagging/Leading strand\n occurrence of sequences")

p4 <- p4 + p_format +
  geom_stripped_cols(odd = "brown", even = "white", alpha = .1) +
  theme(panel.border=element_rect(size=1, fill = NA)) +
  stat_summary(fun = median,
               geom="crossbar",
               size = 0.25,
               width = 0.5,
               color="red") 

ggsave(argv$oCC, width = 11, height = 9, units = "cm") 

counts_total_TC <- calculateAsymmetry(counts_total_df, "TC", "GA")

write.table(counts_total_TC, file = paste0(argv$data_prefix, "S22.csv"), 
            quote = FALSE, row.names = FALSE, sep = ",")

p5 <- ggplot(data=counts_total_TC, aes(x=diff, y=lag_lead)) +
  facet_grid(~repdomains) +
  geom_point(position = position_dodge2(.9), size = .5) +
  scale_x_continuous(breaks = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)) +
  xlab("# of net TCs (TC-GA)") +
  ylab("Lagging/Leading strand\n occurrence of sequences")

p5 <- p5 + p_format +
  geom_stripped_cols(odd = "brown", even = "white", alpha = .1) +
  theme(panel.border=element_rect(size=1, fill = NA)) +
  stat_summary(fun = median,
               geom="crossbar",
               size = 0.25,
               width = 0.5,
               color="red") 

ggsave(argv$oTC, width = 11, height = 9, units = "cm") 


##### Mono- to Penta-mer Plots #####

for (nuc_id in 1:length(counts_total_df$variable)){
  if (str_count(counts_total_df$variable[nuc_id], "A") 
      >= str_count(counts_total_df$variable[nuc_id], "T") ){ 

    temp_base <- counts_total_df$variable[nuc_id]
    counts_total_df$variable[nuc_id] <- counts_total_df$complementary[nuc_id]
    counts_total_df$complementary[nuc_id] <- as.character(temp_base)  
    }
}

counts_total_dcast <- dcast(counts_total_df, direction + repdomains + variable + complementary ~ strand, 
  value.var = "value")

counts_total_dcast$diff <- abs(counts_total_dcast$`+` - counts_total_dcast$`-`)

counts_total_final <- merge(counts_total_df, counts_total_dcast, by=c("direction", "repdomains", 
                                                                      "variable", "complementary"))

counts_total_final <- counts_total_final[, c(1,2,3,4,5,6,7,10)]

colnames(counts_total_final) <- c("direction", "replication_domains", "sequence", "complementary", "strand", 
                                  "percentage", "sd", "diff_strands")

mer1 <- filter(counts_total_final, nchar(as.character(sequence)) == 1)
mer2 <- filter(counts_total_final, nchar(as.character(sequence)) == 2)
mer3 <- filter(counts_total_final, nchar(as.character(sequence)) == 3)
mer4 <- filter(counts_total_final, nchar(as.character(sequence)) == 4)
mer5 <- filter(counts_total_final, nchar(as.character(sequence)) == 5)

#### 1-mer ####

write.table(mer1, file = paste0(argv$data_prefix, "S17.csv"), 
            quote = FALSE, row.names = FALSE, sep = ",")

p.mer1 <- ggplot(data=mer1, aes(x=reorder(sequence, diff_strands), y=percentage, fill=strand)) +
  geom_bar(stat = "identity", size = 1.5,
           position=position_dodge2(padding = 0.05), width=.8)  +
  geom_errorbar(aes(ymin=percentage-sd, ymax=percentage+sd), width=0.2, position=position_dodge(.8)) +
  facet_grid(~replication_domains~direction, scales="free") +
  coord_flip() + 
  scale_y_continuous(breaks = c(-30, -20, -10, 10, 20, 30), 
                     labels = c("30", "20", "10", "10", "20", "30")) +  
  xlab("Sequences") + 
  ylab("Percentage (%)") +
  scale_fill_manual(values = c("+" = "#0571b0", 
                               "-" = "#ca0020")) 

# adding and overriding the default plot format
p.mer1 <- p.mer1 + p_format + 
  theme(panel.spacing.x = unit(-0.5, "lines"))

ggsave(argv$o1, width = 22, height = 18, units = "cm")

#### 2-mer ####

write.table(mer2, file = paste0(argv$data_prefix, "S18.csv"), 
            quote = FALSE, row.names = FALSE, sep = ",")

p.mer2 <- ggplot(data=mer2, aes(x=reorder(sequence, diff_strands), y=percentage, fill=strand)) +
  geom_bar(stat = "identity", size = 1.5,
           position=position_dodge2(padding = 0.05), width=.8)  +
  geom_errorbar(aes(ymin=percentage-sd, ymax=percentage+sd), width=0.2, position=position_dodge(.8)) +
  facet_grid(~replication_domains~direction, scales="free") +
  coord_flip() + 
  scale_y_continuous(breaks = c(-12, -8, -4, 4, 8, 12), 
                     labels = c("12", "8", "4", "4", "8", "12")) + 
  xlab("Sequences") + 
  ylab("Percentage (%)") +
  scale_fill_manual(values = c("+" = "#0571b0", 
                               "-" = "#ca0020")) 

# adding and overriding the default plot format
p.mer2 <- p.mer2 + p_format + 
  theme(panel.spacing.x = unit(-0.5, "lines"))

ggsave(argv$o2, width = 22, height = 18, units = "cm")

#### 3-mer ####

mer3_erd <- filter(mer3, direction == "Left Replicating", strand == "+", replication_domains == "ERD")
mer3_lrd <- filter(mer3, direction == "Left Replicating", strand == "+", replication_domains == "LRD")

mer3_erd$quants <- ntile(mer3_erd$diff_strands, 32)
mer3_lrd$quants <- ntile(mer3_lrd$diff_strands, 32)
mer3_reps <- rbind(mer3_erd, mer3_lrd)
mer3_reps <- mer3_reps[,c(2,3,9)]

mer3 <- merge(mer3, mer3_reps, by=c("replication_domains", "sequence") )

write.table(mer3, file = paste0(argv$data_prefix, "S19.csv"), 
            quote = FALSE, row.names = FALSE, sep = ",")

p.mer3 <- ggplot(data=subset(mer3, quants>22), aes(x=reorder(sequence, diff_strands), y=percentage, fill=strand)) +
  geom_bar(stat = "identity", size = 1.5,
           position=position_dodge2(padding = 0.05), width=.8)  +
  geom_errorbar(aes(ymin=percentage-sd, ymax=percentage+sd), width=0.2, position=position_dodge(.8)) +
  facet_grid(~replication_domains~direction, scales="free") +
  coord_flip() + 
  scale_y_continuous(breaks = c(-5, -4, -3, -2, -1, 1, 2, 3, 4, 5), 
                     labels = c("5", "4", "3", "2", "1", "1", "2", "3", "4", "5")) +
  xlab("Sequences") + 
  ylab("Percentage (%)") +
  scale_fill_manual(values = c("+" = "#0571b0", 
                               "-" = "#ca0020")) 

# adding and overriding the default plot format
p.mer3 <- p.mer3 + p_format + 
  theme(panel.spacing.x = unit(-0.5, "lines"))

ggsave(argv$o3, width = 22, height = 18, units = "cm")

#### 4-mer ####

mer4_erd <- filter(mer4, direction == "Left Replicating", strand == "+", replication_domains == "ERD")
mer4_lrd <- filter(mer4, direction == "Left Replicating", strand == "+", replication_domains == "LRD")

mer4_erd$quants <- ntile(mer4_erd$diff_strands, 136)
mer4_lrd$quants <- ntile(mer4_lrd$diff_strands, 136)
mer4_reps <- rbind(mer4_erd, mer4_lrd)
mer4_reps <- mer4_reps[,c(2,3,9)]

mer4 <- merge(mer4, mer4_reps, by=c("replication_domains", "sequence") )

write.table(mer4, file = paste0(argv$data_prefix, "S20.csv"), 
            quote = FALSE, row.names = FALSE, sep = ",")

p.mer4 <- ggplot(data=subset(mer4, quants>126), aes(x=reorder(sequence, diff_strands), y=percentage, fill=strand)) +
  geom_bar(stat = "identity", size = 1.5,
           position=position_dodge2(padding = 0.05), width=.8)  +
  geom_errorbar(aes(ymin=percentage-sd, ymax=percentage+sd), width=0.2, position=position_dodge(.8)) +
  facet_grid(~replication_domains~direction, scales="free") +
  coord_flip() + 
  scale_y_continuous(breaks = c(-2, -1.5, -1, -0.5, .5, 1, 1.5, 2), 
                     labels = c("2", "1.5", "1", "0.5", "0.5", "1", "1.5", "2")) + 
  xlab("Sequences") + 
  ylab("Percentage (%)") +
  scale_fill_manual(values = c("+" = "#0571b0", 
                               "-" = "#ca0020")) 

# adding and overriding the default plot format
p.mer4 <- p.mer4 + p_format + 
  theme(panel.spacing.x = unit(-0.5, "lines"))

ggsave(argv$o4, width = 22, height = 18, units = "cm")

#### 5-mer ####

mer5_erd <- filter(mer5, direction == "Left Replicating", strand == "+", replication_domains == "ERD")
mer5_lrd <- filter(mer5, direction == "Left Replicating", strand == "+", replication_domains == "LRD")

mer5_erd$quants <- ntile(mer5_erd$diff_strands, 512)
mer5_lrd$quants <- ntile(mer5_lrd$diff_strands, 512)
mer5_reps <- rbind(mer5_erd, mer5_lrd)
mer5_reps <- mer5_reps[,c(2,3,9)]

mer5 <- merge(mer5, mer5_reps, by=c("replication_domains", "sequence") )

write.table(mer5, file = paste0(argv$data_prefix, "S21.csv"), 
            quote = FALSE, row.names = FALSE, sep = ",")

p.mer5 <- ggplot(data=subset(mer5, quants>502), aes(x=reorder(sequence, diff_strands), y=percentage, fill=strand)) +
  #geom_boxplot(outlier.shape = NA) +
  geom_bar(stat = "identity", size = 1.5,
           position=position_dodge2(padding = 0.05), width=.8)  +
  geom_errorbar(aes(ymin=percentage-sd, ymax=percentage+sd), width=0.2, position=position_dodge(.8)) +
  facet_grid(~replication_domains~direction, scales="free") +
  coord_flip() + 
  scale_y_continuous(breaks = c(-0.75, -0.50, -0.25, .25, .5, .75), 
                     labels = c("0.75", "0.5", "0.25", "0.25", "0.5", "0.75")) + 
  xlab("Sequences") + 
  ylab("Percentage (%)") +
  scale_fill_manual(values = c("+" = "#0571b0", 
                               "-" = "#ca0020")) 

# adding and overriding the default plot format
p.mer5 <- p.mer5 + p_format + 
  theme(panel.spacing.x = unit(-0.5, "lines"))

ggsave(argv$o5, width = 22, height = 18, units = "cm")