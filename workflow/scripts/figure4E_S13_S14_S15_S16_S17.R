#### Packages and Libraries ####

library(argparser)
library(stringr)
library(ggpubr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(patchwork)
library(GGally)

######## Arguments ##########
p <- arg_parser("producing figure 4E, S13-S17 - the sequence context of (mono to penta)-mer around initiation zones")
p <- add_argument(p, "-i", help="input")
p <- add_argument(p, "--o1", help="output 1-mer")
p <- add_argument(p, "--o2", help="output 2-mer")
p <- add_argument(p, "--o3", help="output 3-mer")
p <- add_argument(p, "--o4", help="output 4-mer")
p <- add_argument(p, "--o5", help="output 5-mer")
p <- add_argument(p, "--comb", help="output combined")
p <- add_argument(p, "--exp_obs", help="table that contains expected and observed lagging/leading")

# Parse the command line arguments
argv <- parse_args(p)

counts_txt <- argv$i
#counts_txt <-"/Users/azgarian/Desktop/iz_hela_repdomains_uv_mean0.5_windows_201_100_kmer_counts.txt"

#### Default Plot Format ####

#source("/Users/azgarian/Documents/myprojects/replicationRepair/workflow/scripts/plot_format.R")
source("workflow/scripts/plot_format.R")

#### Functions ####

#source("/Users/azgarian/Documents/myprojects/replicationRepair/workflow/scripts/functions.R")
source("workflow/scripts/functions.R")

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

counts_total_dcast_TC <- filter(counts_total_dcast 
                                , direction == "Left Replicating"
                                #, !grepl("A|G",counts_total_dcast$variable)
                                )

counts_total_dcast_TC$lag_lead <- counts_total_dcast_TC$`-` / counts_total_dcast_TC$`+`

baseA <- 29.474993117004857
baseC <- 20.500734650602116
baseG <- 20.515808882191084
baseT <- 29.508463350201946

occurT <- baseA / baseT 
occurC <- baseG / baseC
occurA <- baseT / baseA 
occurG <- baseC / baseG

for (i in 1:length(counts_total_dcast_TC$variable)){

  countA <- str_count(counts_total_dcast_TC$variable[i], "A")
  countT <- str_count(counts_total_dcast_TC$variable[i], "T")
  countC <- str_count(counts_total_dcast_TC$variable[i], "C")
  countG <- str_count(counts_total_dcast_TC$variable[i], "G")
  
  counts_total_dcast_TC$countT[i] <- countT
  counts_total_dcast_TC$`countT-A`[i] <- countT - countA
  counts_total_dcast_TC$`countT+A`[i] <- countT + countA
  counts_total_dcast_TC$`countT+A-G-C`[i] <- countT + countA - countG - countC
  counts_total_dcast_TC$`countT+C`[i] <- countT + countC 
  counts_total_dcast_TC$`countT+C-G-A`[i] <- countT + countC - countG - countA
  
  if (countG > countC & countG > countA & countG > countT) { 
    counts_total_dcast_TC$GorC[i] <- "G-rich" 
  } else if (countC > countG & countC > countA & countC > countT) { 
    counts_total_dcast_TC$GorC[i] <- "C-rich" 
  } else if (countT > countC & countT > countA & countT > countG) { 
    counts_total_dcast_TC$GorC[i] <- "T-rich" 
  } else if (countA > countC & countA > countG & countA > countT) { 
    counts_total_dcast_TC$GorC[i] <- "A-rich" 
  } else { counts_total_dcast_TC$GorC[i] <- "neutral"} 
  
  if (nchar(as.character(counts_total_dcast_TC$variable[i])) > 1){ 
    
    prev_nuc <- str_sub(counts_total_dcast_TC$variable[i], 1, 
                        nchar(as.character(counts_total_dcast_TC$variable[i]))-1)
    #print(prev_nuc)
    rep <- counts_total_dcast_TC$repdomains[i]
    #print(rep)
    prev_obs <- counts_total_dcast_TC$lag_lead[
      counts_total_dcast_TC$variable == prev_nuc & 
        counts_total_dcast_TC$repdomains == rep] 
    if (length(prev_obs) == 0) {
      prev_obs <- counts_total_dcast_TC$lag_lead[
        counts_total_dcast_TC$complementary == prev_nuc & 
          counts_total_dcast_TC$repdomains == rep] 
    }
    #print(prev_obs)
    
    if (str_sub(counts_total_dcast_TC$variable[i], 
                nchar(as.character(counts_total_dcast_TC$variable[i])), -1) == "T") {
      
      counts_total_dcast_TC$expected[i] <- prev_obs * occurT
      
    } else if (str_sub(counts_total_dcast_TC$variable[i], 
                       nchar(as.character(counts_total_dcast_TC$variable[i])), -1) == "C") {
      
      counts_total_dcast_TC$expected[i] <- prev_obs * occurC

    } else if (str_sub(counts_total_dcast_TC$variable[i], 
                       nchar(as.character(counts_total_dcast_TC$variable[i])), -1) == "G") {
      
      counts_total_dcast_TC$expected[i] <- prev_obs * occurG
      
    } else if (str_sub(counts_total_dcast_TC$variable[i], 
                       nchar(as.character(counts_total_dcast_TC$variable[i])), -1) == "A") {
      
      counts_total_dcast_TC$expected[i] <- prev_obs * occurA
    }
     
  } else if (counts_total_dcast_TC$variable[i] == "T") { 
    counts_total_dcast_TC$expected[i] <- occurT
  
  } else if (counts_total_dcast_TC$variable[i] == "C") { 
    counts_total_dcast_TC$expected[i] <- occurC 

  } else if (counts_total_dcast_TC$variable[i] == "A") { 
    counts_total_dcast_TC$expected[i] <- occurA
  
  } else if (counts_total_dcast_TC$variable[i] == "G") { 
    counts_total_dcast_TC$expected[i] <- occurG }

}

counts_total_dcast_TC <- counts_total_dcast_TC[, c(2,3,4,16,8,9,10,11,12,13,14,15)]

colnames(counts_total_dcast_TC) <- c("replication_domains", "sequence", 
                                     "complementary", "expected", "observed",
                                     "countT", "count(T-A)", "count(T+A)",
                                     "count(T+A-G-C)", "count(T+C)", 
                                     "count(T+C-A-G)", "GorC")

write.table(counts_total_dcast_TC, argv$exp_obs, quote = F, row.names = F, sep=",")

counts_total_dcast_TC <- counts_total_dcast_TC %>% 
  group_by(`count(T-A)`) %>% 
  arrange(observed)

p1 <- ggplot(data=counts_total_dcast_TC, aes(x=`count(T-A)`, y=observed)) +
  facet_grid(~replication_domains) +
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

ggsave(argv$comb, width = 11, height = 9, units = "cm") 

counts_total_final <- merge(counts_total_df, counts_total_dcast, by=c("direction", "repdomains", 
                                                                      "variable", "complementary"))

counts_total_final <- counts_total_final[, c(1,2,3,4,5,6,7,10)]

colnames(counts_total_final) <- c("direction", "replication_domains", "sequence", "complementary", "strand", 
                                  "percentage", "sd", "diff_strands")

#write.table(counts_total_final, argv$table, quote = F, row.names = F, sep=",")

mer1 <- filter(counts_total_final, nchar(as.character(sequence)) == 1)
mer2 <- filter(counts_total_final, nchar(as.character(sequence)) == 2)
mer3 <- filter(counts_total_final, nchar(as.character(sequence)) == 3)
mer4 <- filter(counts_total_final, nchar(as.character(sequence)) == 4)
mer5 <- filter(counts_total_final, nchar(as.character(sequence)) == 5)

#### 1-mer ####

p.mer1 <- ggplot(data=mer1, aes(x=reorder(sequence, diff_strands), y=percentage, fill=strand)) +
  #geom_boxplot(outlier.shape = NA) +
  geom_bar(stat = "identity", size = 1.5,
           position=position_dodge2(padding = 0.05), width=.8)  +
  geom_errorbar(aes(ymin=percentage-sd, ymax=percentage+sd), width=0.2, position=position_dodge(.8)) +
  #geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
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

p.mer2 <- ggplot(data=mer2, aes(x=reorder(sequence, diff_strands), y=percentage, fill=strand)) +
  #geom_boxplot(outlier.shape = NA) +
  geom_bar(stat = "identity", size = 1.5,
           position=position_dodge2(padding = 0.05), width=.8)  +
  geom_errorbar(aes(ymin=percentage-sd, ymax=percentage+sd), width=0.2, position=position_dodge(.8)) +
  #geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
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

p.mer3 <- ggplot(data=subset(mer3, quants>22), aes(x=reorder(sequence, diff_strands), y=percentage, fill=strand)) +
  #geom_boxplot(outlier.shape = NA) +
  geom_bar(stat = "identity", size = 1.5,
           position=position_dodge2(padding = 0.05), width=.8)  +
  geom_errorbar(aes(ymin=percentage-sd, ymax=percentage+sd), width=0.2, position=position_dodge(.8)) +
  #geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
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

p.mer4 <- ggplot(data=subset(mer4, quants>126), aes(x=reorder(sequence, diff_strands), y=percentage, fill=strand)) +
  #geom_boxplot(outlier.shape = NA) +
  geom_bar(stat = "identity", size = 1.5,
           position=position_dodge2(padding = 0.05), width=.8)  +
  geom_errorbar(aes(ymin=percentage-sd, ymax=percentage+sd), width=0.2, position=position_dodge(.8)) +
  #geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
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

p.mer5 <- ggplot(data=subset(mer5, quants>502), aes(x=reorder(sequence, diff_strands), y=percentage, fill=strand)) +
  #geom_boxplot(outlier.shape = NA) +
  geom_bar(stat = "identity", size = 1.5,
           position=position_dodge2(padding = 0.05), width=.8)  +
  geom_errorbar(aes(ymin=percentage-sd, ymax=percentage+sd), width=0.2, position=position_dodge(.8)) +
  #geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
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