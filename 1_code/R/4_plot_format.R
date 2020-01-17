
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

windows_lab <- "Relative Position"

chrState_lab <- "Chromatin States"

repdo_lab <- "Replication Domains"

#### ylabs ####

fr_lab <- "RPKM"

fr_xr_ds_lab <- "Repair Rate (log2)"

fr_min_plus_lab <- "Minus/Plus Strand Ratio (log2)"

fr_xr_ds_min_plus_lab <- "Repair Rate Minus/Plus Strand Ratio (log2)"

fr_xr_ds_ear_la_lab <- "Repair Rate Early/Late Phase Ratio (log2)"

#### plot organization ####

p_format <- theme_light() +
            theme(axis.title.x = element_text(size = 14),
                  axis.title.y = element_text(size = 14),
                  axis.text.x = element_text(size = 12, vjust = 0.6),
                  axis.text.y = element_text(size = 12, vjust = 0.1),
                  strip.text.x = element_text(size = 16),
                  strip.text.y = element_text(size = 16, angle = 360),
                  legend.title = element_text(size = 18, face = "bold"),
                  legend.text = element_text(size = 16),
                  legend.position = "bottom")

