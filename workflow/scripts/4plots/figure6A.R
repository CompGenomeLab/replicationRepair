#### Packages and Libraries ####

library(argparser)
library(dplyr)
library(tidyr)
library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

######## Arguments ##########
p <- arg_parser("producing the figure 6A")
p <- add_argument(p, "--hela", help="HeLa initiation zones intersected to GM06990 and IMR90")
p <- add_argument(p, "--gm06990", help="GM06990 initiation zones intersected to HeLa and IMR90")
p <- add_argument(p, "--imr90", help="IMR90 initiation zones intersected to GM06990 and HeLa")
p <- add_argument(p, "--noverlap", help="non-overlapping HeLa initiation zones")
p <- add_argument(p, "--overlap", help="overlapping HeLa initiation zones")
p <- add_argument(p, "--fig6A", help="figure output")

# Parse the command line arguments
argv <- parse_args(p)

hela <- argv$hela

gm06990 <- argv$gm06990

imr90 <- argv$imr90

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
hela_all <- paste(hela_agg$V10,hela_agg$V11,sep=",")


gm_df <- read.table(gm)
gm_df$V10 <- paste("GM06990",gm_df$V1,gm_df$V2,gm_df$V3,gm_df$V4,sep=",")
gm_df$V11 <- paste(gm_df$V5,gm_df$V6,gm_df$V7,gm_df$V8,gm_df$V9,sep=",")
gm_df <- gm_df[,c(10,11)]
gm_agg <- aggregate(V11 ~ V10, data = gm_df, paste, collapse = ",")
gm_sep <- gm_agg %>% separate(V11, c("A","B"), sep=",IMR90,")

gm_others <- filter(gm_sep, grepl(".,.,-1,-1,.|IMR90", A))
gm_others$all <- paste(gm_others$V10,gm_others$A,sep=",")
gm_hela <- filter(gm_sep, !grepl(".,.,-1,-1,.|IMR90", A))
gm_hela_na <- filter(gm_hela, is.na(B))
gm_hela_na$all <- paste(gm_hela_na$A,gm_hela_na$V10,sep=",")
gm_hela_imr <- filter(gm_hela, !is.na(B))
gm_hela_imr$all <- paste(gm_hela_imr$A,gm_hela_imr$V10,"IMR90",gm_hela_imr$B,sep=",")
gm_bind <- rbind(gm_others,gm_hela_na,gm_hela_imr)
gm_all <- as.character(gm_bind$all)


imr90_df <- read.table(imr90)
imr90_df$V10 <- paste("IMR90",imr90_df$V1,imr90_df$V2,imr90_df$V3,imr90_df$V4,sep=",")
imr90_df$V11 <- paste(imr90_df$V5,imr90_df$V6,imr90_df$V7,imr90_df$V8,imr90_df$V9,sep=",")
imr90_df <- imr90_df[,c(10,11)]
imr90_agg <- aggregate(V11 ~ V10, data = imr90_df, paste, collapse = ",")
imr90_all <- paste(imr90_agg$V11,imr90_agg$V10,sep=",")


hela_no_gm <- hela_all[!(hela_all %in% gm_all)]
hela_no_gm_imr <- data.frame(name=hela_no_gm[!(hela_no_gm %in% imr90_all)])
hela_discard <- filter(hela_no_gm_imr, !grepl(",-1,", name))
hela_all2 <- data.frame(name=hela_all)
hela_all <- hela_all2$name[!(hela_all2$name %in% hela_discard$name)]

gm_no_hela <- gm_all[!(gm_all %in% hela_all)]
gm_no_hela_imr <- data.frame(name=gm_no_hela[!(gm_no_hela %in% imr90_all)])
gm_discard <- filter(gm_no_hela_imr, !grepl(",-1,", name))
gm_all2 <- data.frame(name=gm_all)
gm_all <- gm_all2$name[!(gm_all2$name %in% gm_discard$name)]

imr90_no_hela <- imr90_all[!(imr90_all %in% hela_all)]
imr90_no_hela_gm <- data.frame(name=imr90_no_hela[!(imr90_no_hela %in% gm_all)])
imr90_discard <- filter(imr90_no_hela_gm, !grepl(",-1,", name))
imr90_all2 <- data.frame(name=imr90_all)
imr90_all <- imr90_all2$name[!(imr90_all2$name %in% imr90_discard$name)]

hela_gm <- hela_all[hela_all %in% gm_all]
hela_combined <- data.frame(combined=hela_gm[hela_gm %in% imr90_all])
hela_combined <- hela_combined %>% separate(combined, c("name","chr","start","end","score",
                                                        ".1",".2",".3",".4",".5",
                                                        ".6",".7",".8",".9",".10"), sep=",")

hela_only <- hela_combined[,c(2,3,4,1,5)]
write.table(hela_only, file = argv$overlap, sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

temp <- venn.diagram(
  x = list(hela_all, gm_all, imr90_all),
  category.names = c("HeLa-S3" , "GM06990" , "IMR90"),
  #filename = '~/Desktop/intersections_venn_diagramm.pdf',
  #output=TRUE,
  
  # Output features
  #imagetype="pdf" ,
  #height = 1200 , 
  #width = 1200 , 
  #resolution = 800,
  #compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = 2,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 2,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1,
  filename = NULL
)
grid.draw(temp)
pdf(file=argv$fig6A)
  grid.draw(temp)
dev.off()
