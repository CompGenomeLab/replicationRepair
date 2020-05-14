#### libraries ####

library(ggplot2)
library(ggpubr)

#### label names ####

product_labs <- c("CPD", "(6-4)PP")
names(product_labs) <- c("CPD", "64_PP")

method_labs <- c("XR-seq", "Damage-seq", "DNA-seq")
names(method_labs) <- c("XR_seq", "Damage_seq", "DNA_seq")

taex_labs <- c("12 min.", "120 min.")
names(taex_labs) <- c("12", "120")

phase_labs <- c("Async.", "Early \nPhase", "Late \nPhase")
names(phase_labs) <- c("async", "early", "late")

rep_labs <- c("Rep. A", "Rep. B")
names(rep_labs) <- c("A", "B")

dataset_strand_labs <- c("Leftward Direction", "Rightward Direction")
names(dataset_strand_labs) <- c("-", "+")

#### label colors and shapes ####

strand_colors <- c("#ca0020", "#0571b0")

phase_colors <- c("#7fcdbb", "#2c7fb8", "#edf8b1")
# "#7fcdbb" - slightly desaturated cyan, "#2c7fb8" - strong blue, 
# "#edf8b1" - very soft yellow

state_colors <- c("Active Promoter" = "red", 
                 "Promoter Flanking" = "indianred1", 
                 "Inactive Promoter" = "mediumorchid3",
                 "Candidate Strong Enhancer" = "orange", 
                 "Candidate Weak Enhancer" = "yellow",
                 "Distal CTCF/Candidate Insulator" = "turquoise",
                 "Transcription Associated" = "darkgreen", 
                 "Low Activity Proximal to Active States" = "green3", 
                 "Polycob Repressed" = "gray",
                 "Heterochromatin/Repetitive/ \n  Copy Number Variation" = 
                   "white") 

repdomain_colors <- c("#f1a340", "#998ec3", "pink", "turquoise") 
# "#f1a340" - bright orange, "#998ec3" - slightly desaturated blue

#### xlabs ####

windows_lab <- "Relative Position (kb)"

chrState_lab <- "Chromatin States"

repdo_lab <- "Replication Domains"

phase_lab <- "Phases"

#### ylabs ####

fr_lab <- "RPKM"

fr_xr_ds_lab <- "Repair Rate (log2)"

fr_xr_dna_lab <- "XR/DNA (log2)"

fr_ds_dna_lab <- "DS/DNA (log2)"

fr_ear_la_lab <- "Early/Late Phase Ratio (log2)"

fr_plus_min_lab <- "Plus/Minus Strand Ratio (log2)"

fr_xr_ds_plus_min_lab <- "Repair Rate Plus/Minus Strand Ratio (log2)"

fr_xr_ds_ear_la_lab <- "Repair Rate Early/Late Phase Ratio (log2)"

fr_ear_la_plus_min_lab <- "Early/Late Phase, Plus/Minus Strand Ratio (log2)"

#### plot organization ####

p_format <- theme_pubr() +
            theme(axis.title.x = element_text(size = 14),
                  axis.title.y = element_text(size = 14),
                  axis.text.x = element_text(size = 12, vjust = 0.6), 
                                             #hjust=c(0.1, 0.5, 0.9)),
                  axis.text.y = element_text(size = 12, vjust = 0.1),
                  strip.text.x = element_text(size = 16),
                  strip.text.y = element_text(size = 16), #angle = 360),
                  strip.background = element_blank(),
                  legend.title = element_text(size = 18, face = "bold"),
                  legend.text = element_text(size = 16),
                  legend.position = "bottom")

#### plots ####

# RPKM
p_RPKM <- function( df ){
  p <- ggplot(df, aes(x = windows, y = RPKM)) + 
    geom_line(aes(color = sample_strand)) + 
    facet_grid(~product~time_after_exposure~phase~method, 
               labeller = labeller(product = product_labs, 
                                   method = method_labs, 
                                   time_after_exposure = taex_labs, 
                                   replicate = rep_labs, phase = phase_labs)) + 
    xlab(windows_lab) + ylab(fr_lab) +
    scale_x_continuous(limits = c(-half_window, half_window), 
                       breaks = c(-half_window, 0, half_window), 
                       labels = c(paste("-", rlength,  sep = ""), 
                                  "0", 
                                  paste("+", rlength,  sep = ""))) + 
    scale_color_manual(values = strand_colors) + 
    labs(color = "Strands")
  return(p)
}

# Repair Rate (XR/DS)
p_rr <- function( df ){
  phase_labs <- c("Async.", "Early Phase", "Late Phase")
  names(phase_labs) <- c("async", "early", "late")
  p <- ggplot(df, aes(x = windows, y = log2(xr_ds))) + 
    geom_line(aes(color = sample_strand)) + 
    facet_grid(~product~time_after_exposure~phase, 
               labeller = labeller(product = product_labs, 
                                   time_after_exposure = taex_labs, 
                                   replicate = rep_labs, phase = phase_labs)) + 
    xlab(windows_lab) + ylab(fr_xr_ds_lab) +
    scale_x_continuous(limits = c(-half_window, half_window), 
                       breaks = c(-half_window, 0, half_window), 
                       labels = c(paste("-", rlength,  sep = ""), 
                                  "0", 
                                  paste("+", rlength,  sep = ""))) + 
    scale_color_manual(values = strand_colors) + 
    labs(color = "Strands")
  return(p)
}

# Repair Rate Early/Late Phase
p_rr_el <- function( df ){  
  p <- ggplot(df, aes(x=windows, y=log2(ear_la))) +
    geom_line(aes(color = sample_strand)) +
    facet_grid(~product~time_after_exposure~. ,
               labeller = labeller(product = product_labs, 
                                   time_after_exposure = taex_labs, 
                                   replicate = rep_labs)) + 
    xlab(windows_lab) + ylab(fr_xr_ds_ear_la_lab) +
    scale_x_continuous(limits = c(-half_window, half_window), 
                       breaks = c(-half_window, 0, half_window), 
                       labels = c(paste("-", rlength,  sep = ""), 
                                  "0", 
                                  paste("+", rlength,  sep = ""))) + 
    scale_color_manual(values = strand_colors) + 
    labs(color = "Strands")
  return(p)
}

# Repair Rate Plus/Minus Strand
p_rr_pm <- function( df ){  
  p <- ggplot(d, aes(x = windows, y = log2(plus_min))) + 
    geom_line(aes(color = phase)) + 
    geom_line(y=0, color="red", linetype="dashed") +
    facet_grid(~product~time_after_exposure~. , 
               labeller = labeller(product = product_labs, 
                                   time_after_exposure = taex_labs, 
                                   replicate = rep_labs)) + 
    xlab(windows_lab) + ylab(fr_xr_ds_plus_min_lab) +
    scale_x_continuous(limits = c(-half_window, half_window), 
                       breaks = c(-half_window, 0, half_window), 
                       labels = c(paste("-", rlength,  sep = ""), 
                                  "0", 
                                  paste("+", rlength,  sep = ""))) + 
    scale_color_manual(name = "Phase", 
                       label = c("Early Phase", "Late Phase", "Asyncronized"), 
                       values = phase_colors) + 
    labs(color = "Strands")
  return(p)
}

# Early/Late Phase
p_el <- function( df ){  
  p <- ggplot(df, aes(x=windows, y=log2(ear_la))) +
    geom_line(aes(color = sample_strand)) +
    facet_grid(~product~time_after_exposure~method,
               labeller = labeller(product = product_labs, 
                                   time_after_exposure = taex_labs, 
                                   replicate = rep_labs,
                                   method = method_labs)) + 
    xlab(windows_lab) + ylab(fr_ear_la_lab) +
    scale_x_continuous(limits = c(-half_window, half_window), 
                       breaks = c(-half_window, 0, half_window), 
                       labels = c(paste("-", rlength,  sep = ""), 
                                  "0", 
                                  paste("+", rlength,  sep = ""))) + 
    scale_color_manual(values = strand_colors) + 
    labs(color = "Strands")
  return(p)
}

# Plus/Minus Strand of Early/Late Phase 
p_el_pm <- function( df ){  
  p <- ggplot(df, aes(x = windows, y = log2(plus_min))) + 
    geom_line() + 
    geom_line(y=0, color="red", linetype="dashed") +
    facet_grid(~product~time_after_exposure~. , 
               labeller = labeller(product = product_labs, 
                                   time_after_exposure = taex_labs, 
                                   replicate = rep_labs)) + 
    xlab(windows_lab) + ylab(fr_ear_la_plus_min_lab) +
    scale_x_continuous(limits = c(-half_window, half_window), 
                       breaks = c(-half_window, 0, half_window), 
                       labels = c(paste("-", rlength,  sep = ""), 
                                  "0", 
                                  paste("+", rlength,  sep = ""))) 
  return(p)
}

# Plus/Minus Strand 
p_pm <- function( df ){  
  p <- ggplot(d, aes(x = windows, y = log2(plus_min))) + 
    geom_line(aes(color = phase)) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    facet_grid(~product~time_after_exposure~method, 
               labeller = labeller(product = product_labs, 
                                   method = method_labs, 
                                   time_after_exposure = taex_labs, 
                                   replicate = rep_labs)) + 
    xlab(windows_lab) + ylab(fr_plus_min_lab) +
    scale_x_continuous(limits = c(-half_window, half_window), 
                       breaks = c(-half_window, 0, half_window), 
                       labels = c(paste("-", rlength,  sep = ""), 
                                  "0", 
                                  paste("+", rlength,  sep = ""))) + 
    scale_color_manual(name = "Phase", 
                       label = c("Early Phase", "Late Phase", "Asyncronized"), 
                       values = phase_colors) 
  return(p)
}

# XR/DNA
p_xd <- function( df ){
  phase_labs <- c("Async.", "Early Phase", "Late Phase")
  names(phase_labs) <- c("async", "early", "late")
  p <- ggplot(df, aes(x = windows, y = log2(xr_dna))) + 
    geom_line(aes(color = sample_strand)) + 
    facet_grid(~product~time_after_exposure~phase, 
               labeller = labeller(product = product_labs, 
                                   time_after_exposure = taex_labs, 
                                   replicate = rep_labs, phase = phase_labs)) + 
    xlab(windows_lab) + ylab(fr_xr_dna_lab) +
    scale_x_continuous(limits = c(-half_window, half_window), 
                       breaks = c(-half_window, 0, half_window), 
                       labels = c(paste("-", rlength,  sep = ""), 
                                  "0", 
                                  paste("+", rlength,  sep = ""))) + 
    scale_color_manual(values = strand_colors) + 
    labs(color = "Strands")
  return(p)
}

# DS/DNA
p_dd <- function( df ){
  phase_labs <- c("Async.", "Early Phase", "Late Phase")
  names(phase_labs) <- c("async", "early", "late")
  p <- ggplot(df, aes(x = windows, y = log2(ds_dna))) + 
    geom_line(aes(color = sample_strand)) + 
    facet_grid(~product~time_after_exposure~phase, 
               labeller = labeller(product = product_labs, 
                                   time_after_exposure = taex_labs, 
                                   replicate = rep_labs, phase = phase_labs)) + 
    xlab(windows_lab) + ylab(fr_ds_dna_lab) +
    scale_x_continuous(limits = c(-half_window, half_window), 
                       breaks = c(-half_window, 0, half_window), 
                       labels = c(paste("-", rlength,  sep = ""), 
                                  "0", 
                                  paste("+", rlength,  sep = ""))) + 
    scale_color_manual(values = strand_colors) + 
    labs(color = "Strands")
  return(p)
}
