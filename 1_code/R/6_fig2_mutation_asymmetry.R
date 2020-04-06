#### p.A.1.1 Initiation Zones ####

library(ggplot2)
library(stringr)
library(dplyr)
library(reshape2)
library(patchwork)

#### prepare Initiation Zone data

sourcePath <- paste("~/Documents/myprojects/replicationRepair/", 
                    "1_code/R/", sep = "") 

date <- "[2020.02.11.11:20]"

fr_name <- "Initiation_Zones_repdomains_mutations_combined_201_100"

#### set variables 

dateout <- Sys.Date()
dateout <- format(dateout, format = "[%Y.%m.%d]")

mut_dir <- paste("~/Documents/myprojects/replicationRepair/", 
                 "4_output/gitignore/1_TextforPlotting/", sep = "")

#### import and rearrange data 

mut <- read.table(paste(mut_dir, date, fr_name, ".bed", sep = ""))

names(mut) <- c("chr", "start_pos", "end_pos", "name", "score", "strands", "count")

Ccontent <- read.table(paste(mut_dir, 
                             "[2020.02.13]InZones_Cytosine_content_201_100.txt", 
                             sep = ""), header = TRUE)

Ccontent <- Ccontent[,c(1,6,7)]

names(Ccontent) <- c("name","+", "-")

Ccontent_plus <- Ccontent[,c(1,2)]
Ccontent_minus <- Ccontent[,c(1,3)]

Ccontent_plus <- melt(Ccontent_plus, measure.vars = "+")
Ccontent_minus <- melt(Ccontent_minus, measure.vars = "-")

Ccontent <- rbind(Ccontent_plus, Ccontent_minus)

names(Ccontent) <- c("name","strands", "C_content")

mut <- merge(mut, Ccontent, by.x=c("name", "strands"), 
             by.y=c("name", "strands"))

Windows = data.frame(str_split_fixed(mut$name, "_", -1))

mut$window_number = as.numeric(levels(Windows[ , ncol(Windows)]
))[Windows[ , ncol(Windows)]] -101

mut$repdomains <- Windows[,3]

mut$name <- "Initiation_Zones"

mut$norm <- mut$count / mut$C_content

#### plot Initiation Zones 

#### filter the data

mut$repdomains <- factor(mut$repdomains, levels = c("UTZ", "ERD", "DTZ", "LRD"))

#### add plot format

source(paste(sourcePath, "4_plot_format.R", sep = ""))

#### create the plot

p.A.1.1 <- ggplot(mut) + 
  geom_line(aes(x = window_number, y = norm, 
                color = strands)) + 
  facet_grid(~repdomains) +
  xlab("") + ylab("Mutation Count \nper Region per Cytosine") +
  scale_x_continuous(limits = c(-100, 100), 
                     breaks = c(-100, 0, 100), 
                     labels = c("-10 kb", "0", "+10 kb")) + 
  scale_color_manual(values = strand_colors) + 
  ylim(0.002, 0.031) +
  labs(color = "Strands") +
  labs(title = paste("Mutation Profile of Initiation Zones"), 
       subtitle = paste("20 kb sliding windows at 100 bp intervals", sep = "")) 

# add plot format
p.A.1.1 <- p.A.1.1 + p_format + 
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(hjust=c(0.1, 0.5, 0.9))) 

p.A.1.1 # visualize

#### p.A.1.2 ####

mut_plus <- filter(mut, window_number > 0)
mut_plus_agg <- aggregate(x = mut_plus$norm, by = 
                            list(mut_plus$name, mut_plus$repdomains, 
                                 mut_plus$strands), FUN = "mean")
mut_plus_agg$direction <- "Right Replicating"
mut_minus <- filter(mut, window_number < 0)
mut_minus_agg <- aggregate(x = mut_minus$norm, by = 
                            list(mut_minus$name, mut_minus$repdomains, 
                                 mut_minus$strands), FUN = "mean")
mut_minus_agg$direction <- "Left Replicating"
mut_agg <- rbind(mut_plus_agg, mut_minus_agg)

#### create the plot 

p.A.1.2 <- ggplot() + 
  geom_bar(data = mut_agg, aes(x = Group.2, y = x, 
                               fill = Group.2, color = Group.3), 
           stat = "identity", size = 1.5, 
           position=position_dodge2(padding = 0.05)) +
  facet_wrap(~direction) +
  ylim(0, 0.03) +
  xlab("Replication Domains") + 
  ylab("Mutation Count \nper Region \nper Cytosine") +
  scale_fill_manual(values = repdomain_colors, guide = FALSE) + 
  scale_color_manual(values = c("+" = "#0571b0", 
                                "-" = "#ca0020"), 
                     guide = FALSE) +
  labs(color = "Strands", fill = "") 

p.A.1.2 <- p.A.1.2 + p_format + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) # add plot format

p.A.1.2 # visualize

#### p.A.1.3 ####

mut_casted <- dcast(mut, name + repdomains + window_number ~ strands, 
                    value.var = "norm")
mut_casted$plus_min <- mut_casted$"+" / mut_casted$"-" 
mut_plus <- filter(mut_casted, window_number > 0)
mut_plus_agg <- aggregate(x = mut_plus$plus_min, by = 
                            list(mut_plus$name, mut_plus$repdomains), 
                          FUN = "mean")
mut_plus_agg$direction <- "Right Replicating"
mut_minus <- filter(mut_casted, window_number < 0)
mut_minus_agg <- aggregate(x = mut_minus$plus_min, by = 
                             list(mut_minus$name, mut_minus$repdomains), 
                           FUN = "mean")
mut_minus_agg$direction <- "Left Replicating"
mut_agg <- rbind(mut_plus_agg, mut_minus_agg)

#### create the plot

p.A.1.3 <- ggplot() + 
  geom_bar(data = mut_agg, aes(x = Group.2, y = log2(x), fill = Group.2), 
           stat = "identity", position=position_dodge()) +
  facet_wrap(~direction) +
  ylim(-0.25, 0.25) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  xlab("") + ylab("Plus/Minus \n(Log2)") +
  scale_fill_manual(values = repdomain_colors, guide = FALSE) 

# add plot format
p.A.1.3 <- p.A.1.3 + p_format + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())  

p.A.1.3 # visualize


#### p.A.2.1 SNS-seq data ####

sourcePath <- paste("~/Documents/myprojects/replicationRepair/", 
                    "1_code/R/", sep = "") 

date <- "[2020.03.03.18:24]"

fr_name <- "sns_seq_repdomains_mutations_combined_51_100"

#### set variables 

dateout <- Sys.Date()
dateout <- format(dateout, format = "[%Y.%m.%d]")

mut_dir <- paste("~/Documents/myprojects/replicationRepair/", 
                 "4_output/gitignore/1_TextforPlotting/", sep = "")

#### import and rearrange data

mut <- read.table(paste(mut_dir, date, fr_name, ".bed", sep = ""))

names(mut) <- c("chr", "start_pos", "end_pos", "name", "score", "strands", "count")

Ccontent <- read.table(paste(mut_dir, "[2020.03.03]sns_seq_nuc_content_51_100.txt", 
                             sep = ""), header = TRUE)

Ccontent <- Ccontent[,c(1,6,7)]

names(Ccontent) <- c("name","+", "-")

Ccontent_plus <- Ccontent[,c(1,2)]
Ccontent_minus <- Ccontent[,c(1,3)]

Ccontent_plus <- melt(Ccontent_plus, measure.vars = "+")
Ccontent_minus <- melt(Ccontent_minus, measure.vars = "-")

Ccontent <- rbind(Ccontent_plus, Ccontent_minus)

names(Ccontent) <- c("name","strands", "C_content")

mut <- merge(mut, Ccontent, by.x=c("name", "strands"), 
             by.y=c("name", "strands"))

Windows = data.frame(str_split_fixed(mut$name, "_", -1))

mut$window_number = as.numeric(levels(Windows[ , ncol(Windows)]
))[Windows[ , ncol(Windows)]] -26

mut$repdomains <- Windows[,4]

mut$rep <- Windows[,3]

mut$name <- "SNS-seq"

mut$norm <- mut$count / mut$C_content

#### plot SNS-seq

#### filter the data

mut$repdomains <- factor(mut$repdomains, levels = c("UTZ", "ERD", "DTZ", "LRD"))

mut <- filter(mut, rep == "rep1")

#### add plot format

source(paste(sourcePath, "4_plot_format.R", sep = ""))

#### create the plot

p.A.2.1 <- ggplot(mut) + 
  geom_line(aes(x = window_number, y = norm, 
                color = strands)) + 
  facet_grid(~repdomains) +
  xlab("") + ylab("Mutation Count \nper Region per Cytosine") +
  scale_x_continuous(limits = c(-25, 25), 
                     breaks = c(-23, 0, 23), 
                     labels = c("-5 kb", "0", "+5 kb")) + 
  scale_color_manual(values = strand_colors) +
  ylim(0.002, 0.032) +
  labs(color = "Strands") +
  labs(title = paste("Mutation Profile of SNS-seq"), 
       subtitle = paste("10 kb sliding windows at 100 bp intervals", sep = "")) 

p.A.2.1 <- p.A.2.1 + p_format + 
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(hjust=c(0.1, 0.5, 0.9)),
        legend.position = "right") # add plot format

p.A.2.1 # visualize

#### p.A.2.2 ####

mut_plus <- filter(mut, window_number > 0)
mut_plus_agg <- aggregate(x = mut_plus$norm, by = 
                            list(mut_plus$name, mut_plus$repdomains, 
                                 mut_plus$strands), FUN = "mean")
mut_plus_agg$direction <- "Right Replicating"
mut_minus <- filter(mut, window_number < 0)
mut_minus_agg <- aggregate(x = mut_minus$norm, by = 
                             list(mut_minus$name, mut_minus$repdomains, 
                                  mut_minus$strands), FUN = "mean")
mut_minus_agg$direction <- "Left Replicating"
mut_agg <- rbind(mut_plus_agg, mut_minus_agg)

#### create the plot

p.A.2.2 <- ggplot() + 
  geom_bar(data = mut_agg, aes(x = Group.2, y = x, 
                               fill = Group.2, color = Group.3), 
           stat = "identity", size = 1.5, 
           position=position_dodge2(padding = 0.05)) +
  facet_wrap(~direction) +
  ylim(0, 0.03) +
  xlab("") + 
  ylab("Mutation Count \nper Region \nper Cytosine") +
  scale_fill_manual(values = repdomain_colors, guide = FALSE) + 
  scale_color_manual(values = c("+" = "#0571b0", 
                                "-" = "#ca0020"), 
                     guide = FALSE) +
  labs(color = "Strands", fill = "") 

p.A.2.2 <- p.A.2.2 + p_format + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) # add plot format

p.A.2.2 # visualize

#### p.A.2.3 ####

mut_casted <- dcast(mut, name + repdomains + window_number ~ strands, 
                    value.var = "norm")
mut_casted$plus_min <- mut_casted$"+" / mut_casted$"-" 
mut_plus <- filter(mut_casted, window_number > 0)
mut_plus_agg <- aggregate(x = mut_plus$plus_min, by = 
                            list(mut_plus$name, mut_plus$repdomains), 
                          FUN = "mean")
mut_plus_agg$direction <- "Right Replicating"
mut_minus <- filter(mut_casted, window_number < 0)
mut_minus_agg <- aggregate(x = mut_minus$plus_min, by = 
                             list(mut_minus$name, mut_minus$repdomains), 
                           FUN = "mean")
mut_minus_agg$direction <- "Left Replicating"
mut_agg <- rbind(mut_plus_agg, mut_minus_agg)

#### create the plot 

p.A.2.3 <- ggplot() + 
  geom_bar(data = mut_agg, aes(x = Group.2, y = log2(x), fill = Group.2), 
           stat = "identity", position=position_dodge()) +
  facet_wrap(~direction) +
  ylim(-0.25, 0.25) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  xlab("") + ylab("Plus/Minus \n(Log2)") +
  scale_fill_manual(values = repdomain_colors, guide = FALSE) 

p.A.2.3 <- p.A.2.3 + p_format + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) # add plot format

p.A.2.3 # visualize

#### p.B.1.1 Initiation Zone 200kb ####

date <- "[2020.03.03.13:50]"

fr_name <- "inZones_repdomains_mutations_combined_201_1000"

#### set variables 

dateout <- Sys.Date()
dateout <- format(dateout, format = "[%Y.%m.%d]")

mut_dir <- paste("~/Documents/myprojects/replicationRepair/", 
                 "4_output/gitignore/1_TextforPlotting/", sep = "")

#### import and rearrange data 

mut <- read.table(paste(mut_dir, date, fr_name, ".bed", sep = ""))

names(mut) <- c("chr", "start_pos", "end_pos", "name", "score", "strands", "count")

Ccontent <- read.table(paste(mut_dir, 
                             "[2020.03.03]InZones_Cytosine_content_201_1000.txt", 
                             sep = ""), header = TRUE)

Ccontent <- Ccontent[,c(1,6,7)]

names(Ccontent) <- c("name","+", "-")

Ccontent_plus <- Ccontent[,c(1,2)]
Ccontent_minus <- Ccontent[,c(1,3)]

Ccontent_plus <- melt(Ccontent_plus, measure.vars = "+")
Ccontent_minus <- melt(Ccontent_minus, measure.vars = "-")

Ccontent <- rbind(Ccontent_plus, Ccontent_minus)

names(Ccontent) <- c("name","strands", "C_content")

mut <- merge(mut, Ccontent, by.x=c("name", "strands"), 
             by.y=c("name", "strands"))

Windows = data.frame(str_split_fixed(mut$name, "_", -1))

mut$window_number = as.numeric(levels(Windows[ , ncol(Windows)]
))[Windows[ , ncol(Windows)]] -101

mut$repdomains <- Windows[,3]

mut$name <- "Initiation_Zones"

mut$norm <- mut$count / mut$C_content / 10

#### plot Initiation Zones 

#### filter the data 

mut$repdomains <- factor(mut$repdomains, levels = c("UTZ", "ERD", "DTZ", "LRD"))

#### add plot format  

source(paste(sourcePath, "4_plot_format.R", sep = ""))

#### create the plot 

p.B.1.1 <- ggplot(mut) + 
  geom_line(aes(x = window_number, y = norm, 
                color = strands)) + 
  facet_grid(~repdomains) +
  xlab("") + ylab("Mutation Count \nper Region \nper Cytosine") +
  scale_x_continuous(limits = c(-100, 100), 
                     breaks = c(-100, 0, 100), 
                     labels = c("-100 kb", "0", "+100 kb")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands") +
  labs(title = paste("Mutation Profile of Initiation Zones"), 
       subtitle = paste("200 kb sliding windows at 1000 bp intervals", sep = "")) 

p.B.1.1 <- p.B.1.1 + p_format + 
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(hjust=c(0.1, 0.5, 0.9))) # add plot format

p.B.1.1 # visualize

#### p.B.1.2 ####

mut_plus <- filter(mut, window_number > 0)
mut_plus_agg <- aggregate(x = mut_plus$norm, by = 
                            list(mut_plus$name, mut_plus$repdomains, 
                                 mut_plus$strands), FUN = "mean")
mut_plus_agg$direction <- "Right Replicating"
mut_minus <- filter(mut, window_number < 0)
mut_minus_agg <- aggregate(x = mut_minus$norm, by = 
                             list(mut_minus$name, mut_minus$repdomains, 
                                  mut_minus$strands), FUN = "mean")
mut_minus_agg$direction <- "Left Replicating"
mut_agg <- rbind(mut_plus_agg, mut_minus_agg)

#### create the plot 

p.B.1.2 <- ggplot() + 
  geom_bar(data = mut_agg, aes(x = Group.2, y = x, 
                               fill = Group.2, color = Group.3), 
           stat = "identity", size = 1.5, 
           position=position_dodge2(padding = 0.05)) +
  facet_wrap(~direction) +
  xlab("Replication Domains") + 
  ylab("Mutation Count \nper Region \nper Cytosine") +
  scale_fill_manual(values = repdomain_colors, guide = FALSE) + 
  scale_color_manual(values = c("+" = "#0571b0", 
                                "-" = "#ca0020"), 
                     guide = FALSE) +
  labs(color = "Strands", fill = "") 

p.B.1.2 <- p.B.1.2 + p_format + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) + 
  ylim(0, 0.03)  # add plot format

p.B.1.2 # visualize

#### p.B.1.3 ####

mut_casted <- dcast(mut, name + repdomains + window_number ~ strands, 
                    value.var = "norm")
mut_casted$plus_min <- mut_casted$"+" / mut_casted$"-" 
mut_plus <- filter(mut_casted, window_number > 0)
mut_plus_agg <- aggregate(x = mut_plus$plus_min, by = 
                            list(mut_plus$name, mut_plus$repdomains), 
                          FUN = "mean")
mut_plus_agg$direction <- "Right Replicating"
mut_minus <- filter(mut_casted, window_number < 0)
mut_minus_agg <- aggregate(x = mut_minus$plus_min, by = 
                             list(mut_minus$name, mut_minus$repdomains), 
                           FUN = "mean")
mut_minus_agg$direction <- "Left Replicating"
mut_agg <- rbind(mut_plus_agg, mut_minus_agg)

#### create the plot 

p.B.1.3 <- ggplot() + 
  geom_bar(data = mut_agg, aes(x = Group.2, y = log2(x), fill = Group.2), 
           stat = "identity", position=position_dodge()) +
  facet_wrap(~direction) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  xlab("") + ylab("Plus/Minus \n(Log2)") +
  scale_fill_manual(values = repdomain_colors, guide = FALSE) 

p.B.1.3 <- p.B.1.3 + p_format + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) +
  ylim(-0.25, 0.25) # add plot format

p.B.1.3 # visualize


#### p.B.2.1 Initiation Zone 400kb ####

date <- "[2020.03.03.13:30]"

fr_name <- "inZones_repdomains_mutations_combined_201_2000"

#### set variables 

dateout <- Sys.Date()
dateout <- format(dateout, format = "[%Y.%m.%d]")

mut_dir <- paste("~/Documents/myprojects/replicationRepair/", 
                 "4_output/gitignore/1_TextforPlotting/", sep = "")

#### import and rearrange data 

mut <- read.table(paste(mut_dir, date, fr_name, ".bed", sep = ""))

names(mut) <- c("chr", "start_pos", "end_pos", "name", "score", "strands", "count")

Ccontent <- read.table(paste(mut_dir, 
                             "[2020.03.03]inZones_nuc_content_201_2000.txt", 
                             sep = ""), header = TRUE)

Ccontent <- Ccontent[,c(1,6,7)]

names(Ccontent) <- c("name","+", "-")

Ccontent_plus <- Ccontent[,c(1,2)]
Ccontent_minus <- Ccontent[,c(1,3)]

Ccontent_plus <- melt(Ccontent_plus, measure.vars = "+")
Ccontent_minus <- melt(Ccontent_minus, measure.vars = "-")

Ccontent <- rbind(Ccontent_plus, Ccontent_minus)

names(Ccontent) <- c("name","strands", "C_content")

mut <- merge(mut, Ccontent, by.x=c("name", "strands"), 
             by.y=c("name", "strands"))

Windows = data.frame(str_split_fixed(mut$name, "_", -1))

mut$window_number = as.numeric(levels(Windows[ , ncol(Windows)]
))[Windows[ , ncol(Windows)]] -101

mut$repdomains <- Windows[,3]

mut$name <- "Initiation_Zones"

mut$norm <- mut$count / mut$C_content / 20

#### plot Initiation Zones 

#### filter the data 

mut$repdomains <- factor(mut$repdomains, levels = c("UTZ", "ERD", "DTZ", "LRD"))

#### add plot format  

source(paste(sourcePath, "4_plot_format.R", sep = ""))

#### create the plot 

p.B.2.1 <- ggplot(mut) + 
  geom_line(aes(x = window_number, y = norm, 
                color = strands)) + 
  facet_grid(~repdomains) +
  xlab(windows_lab) + ylab("Mutation Count \nper Region \nper Cytosine") +
  scale_x_continuous(limits = c(-100, 100), 
                     breaks = c(-90, 0, 90), 
                     labels = c("-200 kb", "0", "+200 kb")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands") +
  labs(title = paste("Mutation Profile of Initiation Zones"), 
       subtitle = paste("400 kb sliding windows at 2000 bp intervals", sep = "")) 

p.B.2.1 <- p.B.2.1 + p_format + 
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(hjust=c(0.1, 0.5, 0.9))) # add plot format

p.B.2.1 # visualize

#### p.B.2.2 ####

mut_plus <- filter(mut, window_number > 0)
mut_plus_agg <- aggregate(x = mut_plus$norm, by = 
                            list(mut_plus$name, mut_plus$repdomains, 
                                 mut_plus$strands), FUN = "mean")
mut_plus_agg$direction <- "Right Replicating"
mut_minus <- filter(mut, window_number < 0)
mut_minus_agg <- aggregate(x = mut_minus$norm, by = 
                             list(mut_minus$name, mut_minus$repdomains, 
                                  mut_minus$strands), FUN = "mean")
mut_minus_agg$direction <- "Left Replicating"
mut_agg <- rbind(mut_plus_agg, mut_minus_agg)

#### create the plot 

p.B.2.2 <- ggplot() + 
  geom_bar(data = mut_agg, aes(x = Group.2, y = x, 
                               fill = Group.2, color = Group.3), 
           stat = "identity", size = 1.5, 
           position=position_dodge2(padding = 0.05)) +
  facet_wrap(~direction) +
  xlab("Replication Domains") + 
  ylab("Mutation Count \nper Region \nper Cytosine") +
  scale_fill_manual(values = repdomain_colors, guide = FALSE) + 
  scale_color_manual(values = c("+" = "#0571b0", 
                                "-" = "#ca0020"), 
                     guide = FALSE) +
  labs(color = "Strands", fill = "") 

p.B.2.2 <- p.B.2.2 + p_format + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) + 
  ylim(0, 0.03) # add plot format

p.B.2.2 # visualize

#### p.B.2.3 ####

mut_casted <- dcast(mut, name + repdomains + window_number ~ strands, 
                    value.var = "norm")
mut_casted$plus_min <- mut_casted$"+" / mut_casted$"-" 
mut_plus <- filter(mut_casted, window_number > 0)
mut_plus_agg <- aggregate(x = mut_plus$plus_min, by = 
                            list(mut_plus$name, mut_plus$repdomains), 
                          FUN = "mean")
mut_plus_agg$direction <- "Right Replicating"
mut_minus <- filter(mut_casted, window_number < 0)
mut_minus_agg <- aggregate(x = mut_minus$plus_min, by = 
                             list(mut_minus$name, mut_minus$repdomains), 
                           FUN = "mean")
mut_minus_agg$direction <- "Left Replicating"
mut_agg <- rbind(mut_plus_agg, mut_minus_agg)

#### create the plot 

p.B.2.3 <- ggplot() + 
  geom_bar(data = mut_agg, aes(x = Group.2, y = log2(x), fill = Group.2), 
           stat = "identity", position=position_dodge()) +
  facet_wrap(~direction) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  xlab("Replication Domains") + ylab("Plus/Minus \n(Log2)") +
  scale_fill_manual(values = repdomain_colors, guide = FALSE) 

p.B.2.3 <- p.B.2.3 + p_format + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) +
  ylim(-0.25, 0.25) # add plot format

p.B.2.3 # visualize

#### final plot ####


((p.A.2.1 + p.A.2.2 / p.A.2.3) / (p.A.1.1 + p.A.1.2 / p.A.1.3)) / 
  (p.B.1.1 + (p.B.1.2 / p.B.1.3)) / (p.B.2.1 + (p.B.2.2 / p.B.2.3)) +
  plot_layout(guides = "collect", tag_level = 'new') +  
  plot_annotation(tag_levels = c('A'), tag_suffix = ':')

ggsave("~/Desktop/fig2.png", width = 497, height = 410, units = "mm")



