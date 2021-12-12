#### Packages and Libraries ####

library(argparser)
library(stringr)
library(ggpubr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(patchwork)

######## Arguments ##########
p <- arg_parser("producing the sequence context of (1 to 5)-mer on initiation zones")
p <- add_argument(p, "-i", help="input")
p <- add_argument(p, "--o1", help="output 1-mer")
p <- add_argument(p, "--o2", help="output 2-mer")
p <- add_argument(p, "--o3", help="output 3-mer")
p <- add_argument(p, "--o4", help="output 4-mer")
p <- add_argument(p, "--o5", help="output 5-mer")
p <- add_argument(p, "--comb", help="output combined")
p <- add_argument(p, "--table", help="table that contains all the values")

# Parse the command line arguments
argv <- parse_args(p)

counts_txt <- argv$i
#counts_txt <-"/Users/azgarian/Desktop/iz_hela_repdomains_uv_mean0.5_windows_201_100_kmer_counts.txt"

#### Default Plot Format ####

source("workflow/scripts/plot_format.R")

#### Functions ####

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

counts_total_dcast <- dcast(counts_total_df, direction + repdomains + variable + complementary ~ strand, 
  value.var = "value")

counts_total_dcast$diff <- abs(counts_total_dcast$`+` - counts_total_dcast$`-`)

counts_total_final <- merge(counts_total_df, counts_total_dcast, by=c("direction", "repdomains", 
  "variable", "complementary"))

counts_total_final <- counts_total_final[, c(1,2,3,4,5,6,7,10)]

colnames(counts_total_final) <- c("direction", "replication_domains", "sequence", "complementary", "strand", 
  "percentage", "sd", "diff_strands")

write.table(counts_total_final, argv$table, quote = F, row.names = F, sep=",")

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


#### plot C #####
# plot C will be a drawing
# we are creating an empty text for that part
p.C <- wrap_elements(grid::textGrob(''))

#### Combining Plots with Patchwork ####

layout <- "
AABB
AACC
"

p.mer3 + p.mer1 + p.C +
  plot_annotation(tag_levels = 'A') +
  plot_layout(design = layout, guides="collect") & 
  theme(plot.tag = element_text(hjust = -0.2, size = 12, face="bold"),
    legend.position = "bottom")

ggsave(argv$comb, width = 22, height = 18, units = "cm") 
