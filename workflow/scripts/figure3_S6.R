#### Packages and Libraries ####

library(ggpubr)
library(ggplot2)
library(stringr)
library(reshape2)
library(dplyr)
library(patchwork)
library(argparser)
set.seed(1) 

######## Arguments ##########
p <- arg_parser("producing the figure 3 and suplementary figure 6")
p <- add_argument(p, "--df", help="region file with read counts")
p <- add_argument(p, "--df_sim", help="region file with simulated read counts")
p <- add_argument(p, "--dtype", help="damage type of the sample to visualize")
p <- add_argument(p, "--data_prefix", help="name prefix of the dataframes that generate the plots")
p <- add_argument(p, "-o", help="output")

# Parse the command line arguments
argv <- parse_args(p)

sample_csv <- argv$df

sample_csv_sim <- argv$df_sim

#### Variables ####

rep <- "_" # replicate
dprod <- argv$dtype # damage product (CPD/64_PP)

#### Default Plot Format ####

source("workflow/scripts/plot_format.R")


#### Fuctions ####

source("workflow/scripts/functions.R")

#### Main ####

# for both plot A and plot B
sample_df <- read.delim( sample_csv, header = F )

colnames(sample_df) <- c("chromosomes", "start_position", "end_position", 
                            "dataset", "score", "dataset_strand", "counts", 
                            "sample_names", "file_names", "layout", "cell_line", 
                            "product", "method", "uv_exposure", "treatment", 
                            "phase", "time_after_exposure", "replicate", 
                            "project", "sample_source", "sample_strand", 
                            "mapped_reads", "RPKM")

sample_df_filt <- filter(sample_df, product == dprod, replicate == rep, 
                         phase != "async", time_after_exposure == "12")

sample_df_rmlow <- rmv_low_counted_regions(sample_df_filt, 2)

df_rr <- repair_rate( sample_df_rmlow )

colnames(df_rr)[12] <- "real"

sample_df_sim <- read.delim( sample_csv_sim, header = F )

colnames(sample_df_sim) <- c("chromosomes", "start_position", "end_position", 
                         "dataset", "score", "dataset_strand", "counts", 
                         "sample_names", "file_names", "layout", "cell_line", 
                         "product", "method", "uv_exposure", "treatment", 
                         "phase", "time_after_exposure", "replicate", 
                         "project", "sample_source", "sample_strand", 
                         "mapped_reads", "RPKM")

sample_df_sim_filt <- filter(sample_df_sim, product == dprod, replicate == rep, 
                         phase != "async", time_after_exposure == "12")

sample_df_sim_rmlow <- rmv_low_counted_regions(sample_df_sim_filt, 2)

df_rr_sim <- repair_rate( sample_df_sim_rmlow )

colnames(df_rr_sim)[12] <- "sim"

df_rr_rs <- merge(df_rr, df_rr_sim, by = c("chromosomes", "start_position", 
                                            "end_position", "dataset", "score",
                                            "dataset_strand", "product", 
                                            "phase", "time_after_exposure", 
                                            "sample_strand", "replicate" ))
df_rr_rs$xr_ds <- df_rr_rs$real / df_rr_rs$sim

# for plot A
pA_df_rr <- chrState_naming( df_rr_rs, chr_states, general_states, 
                             chrState2generalState )

pA_df_rr$state_short <- NA
for (i in 1:length(chr_states)) {
  pA_df_rr <- within(pA_df_rr, state_short[states == general_states[i]] <- 
                       general_states_short[i])
}

# for plot B
df_rr_ear_la <- rrEarly_rrLate( df_rr_rs )
pB_df_rr_ear_la <- chrState_naming( df_rr_ear_la, chr_states, general_states, 
                                    chrState2generalState )

pB_df_rr_ear_la$state_short <- NA
for (i in 1:length(chr_states)) {
  pB_df_rr_ear_la <- within(pB_df_rr_ear_la, state_short[states == general_states[i]] <- 
                       general_states_short[i])
}

#### Filtering Samples ####

# for plot A
pA_filt <- filter(pA_df_rr, xr_ds != "NaN", xr_ds != "Inf", 
                  dataset_strand != "DTZ", dataset_strand != "UTZ", 
                  product == dprod, replicate == rep, phase != "async", 
                  time_after_exposure == "12")

pA_filt$state_short <- factor(pA_filt$state_short, levels=general_states_short)

# for plot B
pB_filt <- filter(pB_df_rr_ear_la, ear_la != "NaN", ear_la != "Inf", 
                  dataset_strand != "DTZ", dataset_strand != "UTZ", 
                  product == dprod, replicate == rep,
                  time_after_exposure == "12")

pB_filt$state_short <- factor(pB_filt$state_short, levels=general_states_short)

pB_filt$ear_la <- log2(pB_filt$ear_la)

stat.test = compare_means(ear_la ~ 1, paired = FALSE, data = pB_filt, method = "wilcox.test", 
                          group.by = "state_short", mu = 0) %>%
  mutate(y.position = 4.5)

#### Plot A ####

pA_filt$states <- NA

write.table(pA_filt, file = paste0(argv$data_prefix, "A.csv"), quote = FALSE, 
            row.names = FALSE, sep = ",")

pA_filt <- chrState_naming( pA_filt, chr_states, general_states, 
                             chrState2generalState )

pA_filt$phase[pA_filt$phase == "early"] <- "Early \nS Phase"
pA_filt$phase[pA_filt$phase == "late"] <- "Late \nS Phase"

# create the plot
p.A <- ggplot(pA_filt, aes(x = state_short, y = log2(xr_ds))) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_boxplot(aes(fill=factor(states)), 
               outlier.shape = NA, lwd=0.5) +
  facet_grid(~phase~dataset_strand) +
  xlab("") + ylab("n. Repair Rate (RR) (log2)") +
  scale_fill_manual(name = "Chromatin States\n(#s in ERD/LRD)", values = state_colors) +
  #ylim(-4,5) +
  guides(size = "none") 

# adding and overriding the default plot format
p.A <- p.A + p_format +
  theme(panel.border=element_rect(size=1, fill = NA),
        strip.text = element_blank(),
        legend.position = "right", 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

#### Plot B ####

pB_filt$states <- NA

write.table(pB_filt, file = paste0(argv$data_prefix, "B.csv"), quote = FALSE, 
            row.names = FALSE, sep = ",")

pB_filt <- chrState_naming( pB_filt, chr_states, general_states, 
                             chrState2generalState )

# create the plot
p.B <- ggplot(pB_filt, aes(x = state_short, y = ear_la)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_boxplot(aes(fill = factor(states)), outlier.shape = NA, lwd=0.5) +
  facet_wrap(~dataset_strand) +
  xlab(chrState_lab) + ylab(expression(RR[E] / RR[L] (log2))) +
  scale_fill_manual(name = "Chromatin States\n(#s in ERD/LRD)", values = state_colors) +
  #ylim(-1,1) +
  guides(fill="none") 

# adding and overriding the default plot format
p.B <- p.B + p_format +
  theme(panel.border=element_rect(size=1, fill = NA),
        strip.text = element_blank(),
        legend.position = "right", 
        axis.text.x = element_text(angle = 60, vjust = 0.95, hjust = 0.95)) +
  geom_text(data = stat.test, aes(y=5, label = p.signif))

p.A / p.B + 
  plot_layout(guides = "collect") & 
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 12, face="bold"),
        legend.margin=margin())

ggsave(argv$o, width = 22, height = 18, units = "cm")
