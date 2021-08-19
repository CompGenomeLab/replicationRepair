#### Packages and Libraries ####

library(stringr)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)
library(grid)
library(ggthemes)


#### Variables ####

# name of the sample csv file for plot A
real_name <- paste("/Users/azgarian/Documents/myprojects/replicationRepair/results/",
                   "final/final_reports_hg19_iz_repdomains_hela_windows_201_100_ready.csv",
                   sep = "")

sim_name <- paste("/Users/azgarian/Documents/myprojects/replicationRepair/results/",
                  "final/final_reports_sim_hg19_iz_repdomains_hela_windows_201_100_ready.csv", 
                  sep = "")

name="iz"

#### Default Plot Format ####

source("/Users/azgarian/Documents/myprojects/replicationRepair/workflow/scripts/4plots/4_plot_format.R")

# for naming of simulated samples
method_labs_sim <- c("Simulated \nDamage-seq", "Simulated \nXR-seq")
names(method_labs_sim) <- c("Damage_seq", "XR_seq")

#### Functions ####

source("/Users/azgarian/Documents/myprojects/replicationRepair/workflow/scripts/4plots/4_functions.R")

repair_rate <- function ( df, ds_rep = "A" ){
  ds <- filter(df, method == "Damage_seq" & replicate == ds_rep)
  xr <- filter(df, method == "XR_seq")
  xr_list <- split(xr, xr$replicate)
  df_xr_ds <- xr_list[[1]][0, ]
  
  for ( xr_rep in 1:length(xr_list) ){
    temp <- xr_list[[xr_rep]]
    temp <- rbind(temp, ds) 
    temp <- dcast(temp, chromosomes + start_position + end_position + 
                    dataset + score + dataset_strand + product + windows +
                    cell_line + time_after_exposure + treatment + repdomains +
                    sample_strand + phase ~ method, 
                  value.var = "RPKM")
    temp <- cbind(temp, xr_list[[xr_rep]]["replicate"])
    df_xr_ds <- rbind(df_xr_ds, temp)
  }
  
  df_xr_ds$xr_ds <- df_xr_ds$XR_seq / df_xr_ds$Damage_seq
  df_xr_ds <- select(df_xr_ds, -c("Damage_seq", "XR_seq"))
  
  return(df_xr_ds)
}

real_sim <- function ( df, sim_rep = "A" ){
  sim <- filter(df, type == "sim" & replicate == sim_rep)
  real <- filter(df, type == "real")
  real_list <- split(real, real$replicate)
  df_real_sim <- real_list[[1]][0, ]
  
  for ( real_rep in 1:length(real_list) ){
    temp <- real_list[[real_rep]]
    temp <- rbind(temp, sim) 
    temp <- dcast(temp, chromosomes + start_position + end_position + 
                    dataset + score + dataset_strand + product + windows +
                    cell_line + treatment + method + phase + repdomains +
                    sample_strand + time_after_exposure ~ type, 
                  value.var = "RPKM")
    temp <- cbind(temp, real_list[[real_rep]]["replicate"])
    df_real_sim <- rbind(df_real_sim, temp)
  }
  
  df_real_sim$real_sim <- df_real_sim$real / df_real_sim$sim
  df_real_sim <- select(df_real_sim, -c("real", "sim"))
  
  return(df_real_sim)
}

rr_rs <- function ( df, ds_rep = "A" ){
  ds <- filter(df, method == "Damage_seq" & replicate == ds_rep)
  xr <- filter(df, method == "XR_seq")
  xr_list <- split(xr, xr$replicate)
  df_xr_ds <- xr_list[[1]][0, ]
  
  for ( xr_rep in 1:length(xr_list) ){
    temp <- xr_list[[xr_rep]]
    temp <- rbind(temp, ds) 
    temp <- dcast(temp, chromosomes + start_position + end_position + 
                    dataset + score + dataset_strand + product + windows +
                    cell_line + time_after_exposure + treatment + repdomains +
                    sample_strand + type + phase ~ method, 
                  value.var = "RPKM")
    temp <- cbind(temp, xr_list[[xr_rep]]["replicate"])
    df_xr_ds <- rbind(df_xr_ds, temp)
  }
  
  df_xr_ds$xr_ds <- df_xr_ds$XR_seq / df_xr_ds$Damage_seq
  df_xr_ds <- select(df_xr_ds, -c("Damage_seq", "XR_seq"))
  
  sim <- filter(df_xr_ds, type == "sim" & replicate == ds_rep)
  real <- filter(df_xr_ds, type == "real")
  real_list <- split(real, real$replicate)
  df_real_sim <- real_list[[1]][0, ]
  
  for ( real_rep in 1:length(real_list) ){
    temp <- real_list[[real_rep]]
    temp <- rbind(temp, sim) 
    temp <- dcast(temp, chromosomes + start_position + end_position + 
                    dataset + score + dataset_strand + product + windows +
                    cell_line + treatment + phase + repdomains +
                    sample_strand + time_after_exposure ~ type, 
                  value.var = "xr_ds")
    temp <- cbind(temp, real_list[[real_rep]]["replicate"])
    df_real_sim <- rbind(df_real_sim, temp)
  }
  
  df_real_sim$real_sim <- df_real_sim$real / df_real_sim$sim
  df_real_sim <- select(df_real_sim, -c("real", "sim"))
  
  return(df_real_sim)
}

#### Main ####

# real
real_df <- read.csv( real_name)

window <- cbind(real_df$dataset, real_df$sample_names, real_df$sample_strand, 
                data.frame(str_split_fixed(real_df$dataset, "_", -1)))
window <- window[ , c(1, 6, 2, 3, ncol(window))]
names(window) <- c("dataset","repdomains", "sample_names","sample_strand", "windows")
window$windows <- as.numeric(as.character(window$windows)) - 101
real_df_org <- merge(real_df, window, by=c("dataset","sample_names","sample_strand"))


real_df_org$dataset <- gsub("_.*", "", real_df_org$dataset)

real_df_org$sample_strand <- factor(
  real_df_org$sample_strand, levels = c("+","-"))

real_df_org$repdomains <- factor(
  real_df_org$repdomains, levels = c("UTZ","ERD","LRD","DTZ"))

# for sim
sim_df <- read.csv( sim_name )

window <- cbind(sim_df$dataset, sim_df$sample_names, sim_df$sample_strand, 
                data.frame(str_split_fixed(sim_df$dataset, "_", -1)))
window <- window[ , c(1, 6, 2, 3, ncol(window))]
names(window) <- c("dataset","repdomains", "sample_names","sample_strand", "windows")
window$windows <- as.numeric(as.character(window$windows)) - 101
sim_df_org <- merge(sim_df, window, by=c("dataset","sample_names","sample_strand"))


sim_df_org$dataset <- gsub("_.*", "", sim_df$dataset)

sim_df_org$sample_strand <- factor(
  sim_df_org$sample_strand, levels = c("+","-"))

sim_df_org$repdomains <- factor(
  sim_df_org$repdomains, levels = c("UTZ","ERD","LRD","DTZ"))

#### He-La ####

# filtering for A.1
real_df_org <- filter(real_df_org, product == "CPD", phase== "early", replicate=="A")

sim_df_org <- filter(sim_df_org, product == "CPD", phase== "early", replicate=="A")

hela_real_rr <- repair_rate(real_df_org)

hela_sim_rr <- repair_rate(sim_df_org)

real_df_org$type <- "real"
sim_df_org$type <- "sim"

df_org <- rbind(real_df_org, sim_df_org)

df_rs <- real_sim(df_org)
df_rr_rs <- rr_rs(df_org)

#### Plot hela_real ####

# create the plot 
p_hela_real <- ggplot(real_df_org, aes(x = windows, y = RPKM)) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~product~time_after_exposure~phase~method~repdomains) + 
  xlab("Relative Position to Initiation Zones (kb)") + ylab(fr_lab) +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")

# adding and overriding the default plot format
p_hela_real <- p_hela_real + p_format 

p_hela_real

ggsave(paste("~/Desktop/figs/p_", name, "_real.png", sep=""), width = 22, height = 18, units = "cm")

#### Plot hela_sim ####

# create the plot 
p_hela_sim <- ggplot(sim_df_org, aes(x = windows, y = RPKM)) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~product~time_after_exposure~phase~method~repdomains) + 
  xlab("Relative Position to Initiation Zones (kb)") + ylab(fr_lab) +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")

# adding and overriding the default plot format
p_hela_sim <- p_hela_sim + p_format 

p_hela_sim

ggsave(paste("~/Desktop/figs/p_", name, "_sim.png", sep=""), width = 22, height = 18, units = "cm")

#### Plot hela_real_rr ####

# create the plot 
p_hela_real_rr <- ggplot(hela_real_rr, aes(x = windows, y = log2(xr_ds))) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~product~time_after_exposure~phase~repdomains) + 
  xlab("Relative Position to Edu Peaks (kb)") + ylab("Repair Rate (log2)") +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")

# adding and overriding the default plot format
p_hela_real_rr <- p_hela_real_rr + p_format 

p_hela_real_rr

ggsave(paste("~/Desktop/figs/p_", name, "_real_rr.png", sep=""), width = 22, height = 18, units = "cm")

#### Plot hela_sim_rr ####

# create the plot 
p_hela_sim_rr <- ggplot(hela_sim_rr, aes(x = windows, y = log2(xr_ds))) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~product~time_after_exposure~phase~repdomains) + 
  xlab("Relative Position to Initiation Zones (kb)") + ylab("Repair Rate (log2)") +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")

# adding and overriding the default plot format
p_hela_sim_rr <- p_hela_sim_rr + p_format 

p_hela_sim_rr

ggsave(paste("~/Desktop/figs/p_", name, "_sim_rr.png", sep=""), width = 22, height = 18, units = "cm")

#### Plot real/sim ####

# create the plot 
p_df_rs <- ggplot(df_rs, aes(x = windows, y = log2(real_sim))) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~product~time_after_exposure~phase~method~repdomains) + 
  xlab("Relative Position to Initiation Zones (kb)") + ylab("Real/Simulation (log2)") +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")

# adding and overriding the default plot format
p_df_rs <- p_df_rs + p_format 

p_df_rs

ggsave(paste("~/Desktop/figs/p_", name, "_rs.png", sep=""), width = 22, height = 18, units = "cm")

#### Plot rr real/sim ####

# create the plot 
p_df_rs <- ggplot(df_rr_rs, aes(x = windows, y = log2(real_sim))) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~product~time_after_exposure~phase~repdomains) + 
  xlab("Relative Position to Initiation Zones (kb)") + ylab("Repair Rate Real/Simulation (log2)") +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")

# adding and overriding the default plot format
p_df_rs <- p_df_rs + p_format 

p_df_rs

ggsave(paste("~/Desktop/figs/p_", name, "_rr_rs.png", sep=""), width = 22, height = 18, units = "cm")

