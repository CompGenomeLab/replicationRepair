#### Packages and Libraries ####

library(argparser)
library(dplyr)
library(tidyr)

######## Arguments ##########
p <- arg_parser("producing non-overlapping/overlapping hela initiation zones")
p <- add_argument(p, "--hela", help="HeLa initiation zones intersected to GM06990 and IMR90")
p <- add_argument(p, "--noverlap", help="non-overlapping HeLa initiation zones")
p <- add_argument(p, "--overlap", help="overlapping HeLa initiation zones")

# Parse the command line arguments
argv <- parse_args(p)

hela <- argv$hela

#hela <- paste0("/home/azgarian/Desktop/",
#               "SRR2913039_intersect2_SRR2913063_ERR2760855.txt")

#gm <- paste0("/home/azgarian/Desktop/",
#             "SRR2913063_intersect2_SRR2913039_ERR2760855.txt")

#imr90 <- paste0("/home/azgarian/",
#                "ERR2760855_intersect2_SRR2913039_SRR2913063.txt")

hela_df <- read.table(hela)

hela_no_overlap <- filter(hela_df, V7 == -1)
hela_no_overlap$name <- "HeLa" 
hela_no_overlap <- hela_no_overlap[,c(1,2,3,11,4)]
write.table(hela_no_overlap, file = argv$noverlap, sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

hela_df$V10 <- paste("HeLa",hela_df$V1,hela_df$V2,hela_df$V3,hela_df$V4,sep=",")
hela_df$V11 <- paste(hela_df$V5,hela_df$V6,hela_df$V7,hela_df$V8,hela_df$V9,sep=",")
hela_df <- hela_df[,c(10,11)]
hela_agg <- aggregate(V11 ~ V10, data = hela_df, paste, collapse = ",")
hela_overlap <- filter(hela_agg, grepl("GM06990.*IMR90|IMR90.*GM06990", V11))
hela_final <- separate(hela_overlap, V10, c("name","chr","start","end","score"), sep=",")
hela_final <- hela_final[,c(2,3,4,1,5)]
write.table(hela_final, file = argv$overlap, sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)