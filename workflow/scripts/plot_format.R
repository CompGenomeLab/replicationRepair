#### libraries ####

library(ggplot2)
library(ggpubr)
set.seed(1) 

#### label names ####

chr_states <- c("Tss", "TssF", "PromF", "PromP", "Enh", "EnhF", "EnhWF", 
                "EnhW", "DnaseD", "DnaseU", "FaireW", "CtcfO", "Ctcf", "Gen5'", "Elon", "ElonW", "Gen3'", 
                "Pol2", "H4K20", "Low", "ReprD", "Repr", "ReprW", 
                "Quies", "Art")

general_states <- c("Active Promoter\n(1593/41)\n", 
                    "Promoter Flanking\n(415/49)\n", 
                    "Inactive Promoter\n(16/2)\n", 
                    "Candidate Strong\nEnhancer\n(1001/23)\n", 
                    "Candidate Weak\nEnhancer\n(1491/123)\n", 
                    "Distal CTCF/Candidate\nInsulator\n(34/9)\n", 
                    "Transcription Associated\n(8735/313)\n", 
                    "Low Activity Proximal\nto Active States\n(28714/21984)\n", 
                    "Polycomb Repressed\n(702/3487)\n", 
                    "Heterochromatin/\nRepetitive/\nCopy Number Variation\n(2102/9701)\n")

general_states_short <- c("Active Promoter", 
                    "Promoter Flanking", 
                    "Inactive Promoter", 
                    "C. Strong Enh.", 
                    "C. Weak Enh.", 
                    "Distal CTCF/C. Insulator", 
                    "Transcription Associated", 
                    "Low Activity Proximal", 
                    "Polycomb Repressed", 
                    "HET/Rep./C.N.V.")

chrState2generalState <- c("Active Promoter\n(1593/41)\n", 
                           "Active Promoter\n(1593/41)\n", 
                           "Promoter Flanking\n(415/49)\n", 
                           "Inactive Promoter\n(16/2)\n", 
                           "Candidate Strong\nEnhancer\n(1001/23)\n", 
                           "Candidate Strong\nEnhancer\n(1001/23)\n", 
                           "Candidate Weak\nEnhancer\n(1491/123)\n", 
                           "Candidate Weak\nEnhancer\n(1491/123)\n",
                           "Candidate Weak\nEnhancer\n(1491/123)\n",
                           "Candidate Weak\nEnhancer\n(1491/123)\n", 
                           "Candidate Weak\nEnhancer\n(1491/123)\n", 
                           "Distal CTCF/Candidate\nInsulator\n(34/9)\n", 
                           "Distal CTCF/Candidate\nInsulator\n(34/9)\n",
                           "Transcription Associated\n(8735/313)\n", 
                           "Transcription Associated\n(8735/313)\n", 
                           "Transcription Associated\n(8735/313)\n", 
                           "Transcription Associated\n(8735/313)\n", 
                           "Transcription Associated\n(8735/313)\n", 
                           "Transcription Associated\n(8735/313)\n", 
                           "Low Activity Proximal\nto Active States\n(28714/21984)\n", 
                           "Polycomb Repressed\n(702/3487)\n", 
                           "Polycomb Repressed\n(702/3487)\n", 
                           "Polycomb Repressed\n(702/3487)\n", 
                           "Heterochromatin/\nRepetitive/\nCopy Number Variation\n(2102/9701)\n", 
                           "Heterochromatin/\nRepetitive/\nCopy Number Variation\n(2102/9701)\n")

product_labs <- c("CPD", "(6-4)PP")
names(product_labs) <- c("CPD", "64_PP")

method_labs <- c("XR-seq", "Damage-seq", "DNA-seq")
names(method_labs) <- c("XR_seq", "Damage_seq", "DNA_seq")

taex_labs <- c("0 min.", "12 min.", "120 min.", "60 min.")
names(taex_labs) <- c("0", "12", "120", "60")

phase_labs <- c("Async.", "Early S\nPhase", "Late S\nPhase")
names(phase_labs) <- c("async", "early", "late")

rep_labs <- c("Rep. A", "Rep. B")
names(rep_labs) <- c("A", "B")

dataset_strand_labs <- c("Leftward Direction", "Rightward Direction")
names(dataset_strand_labs) <- c("-", "+")

#### label colors and shapes ####

strand_colors <- c("#0571b0", "#ca0020")

phase_colors <- c("#7fcdbb", "#2c7fb8", "#edf8b1")
# "#7fcdbb" - slightly desaturated cyan, "#2c7fb8" - strong blue, 
# "#edf8b1" - very soft yellow

state_colors <- c("Active Promoter\n(1593/41)\n" = "red", 
                 "Promoter Flanking\n(415/49)\n" = "indianred1", 
                 "Inactive Promoter\n(16/2)\n" = "mediumorchid3",
                 "Candidate Strong\nEnhancer\n(1001/23)\n" = "orange", 
                 "Candidate Weak\nEnhancer\n(1491/123)\n" = "yellow",
                 "Distal CTCF/Candidate\nInsulator\n(34/9)\n" = "turquoise",
                 "Transcription Associated\n(8735/313)\n" = "darkgreen", 
                 "Low Activity Proximal\nto Active States\n(28714/21984)\n" = "green3", 
                 "Polycomb Repressed\n(702/3487)\n" = "gray",
                 "Heterochromatin/\nRepetitive/\nCopy Number Variation\n(2102/9701)\n" = 
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

fr_ear_la_lab <- "Early/Late S Phase Ratio (log2)"

fr_plus_min_lab <- "Plus/Minus Strand Ratio (log2)"

fr_xr_ds_plus_min_lab <- "Repair Rate Plus/Minus Strand Ratio (log2)"

fr_xr_ds_ear_la_lab <- 
  "Relative Difference Between Repair Rates of\nEarly and Late S Phases (log2)"

fr_ear_la_plus_min_lab <- "Early/Late S Phase, Plus/Minus Strand Ratio (log2)"

#### plot organization ####

options(scipen = 999)

p_format <- theme_pubr() +
            theme(plot.title = element_text(size = 12),
                  axis.title.x = element_text(size = 12),
                  axis.title.y = element_text(size = 12),
                  axis.text.x = element_text(size = 10, vjust = 0.6), 
                  axis.text.y = element_text(size = 10, vjust = 0.1),
                  strip.text.x = element_text(size = 12),
                  strip.text.y = element_text(size = 12),
                  strip.background = element_blank(),
                  legend.title = element_text(size = 10, face = "bold"),
                  legend.text = element_text(size = 8),
                  legend.position = "bottom")

