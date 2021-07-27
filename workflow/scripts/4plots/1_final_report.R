#### library #####


#### import data ####

dir_fr <- "~/Documents/myprojects/replicationRepair/final/"

temp <- list.files(path = dir_fr, pattern = ".txt")

temp <- temp[grepl("*final_report.*\\.txt", temp)]

setwd(dir_fr)

for (report_name in temp) {

fr <- read.delim(report_name, header = FALSE, sep = "\t")

colnames(fr) <- c("chromosomes", "start_position", "end_position", "dataset", 
                  "score", "dataset_strand", "counts", "sample_names", 
                  "file_names", "layout", "cell_line", "product", "method", 
                  "uv_exposure", "treatment", "phase", "time_after_exposure", 
                  "replicate", "project", "sample_source", "sample_strand", 
                  "mapped_reads", "RPKM")

#### data information ####

new_name <- sub(".txt", "", report_name)

length <- as.numeric(fr$end_position) - as.numeric(fr$start_position)

write("", file = paste(new_name, "_info.TXT", sep = ""))

write(paste("project:", list(levels(fr$project)), "\n\nreport file name:", 
            report_name, "\n\nlength min:", min(length), "\n\nlength max:", 
            max(length), "\n\nnumber of datasets:", nlevels(fr$dataset), 
            "\n\ndataset:", list(levels(fr$dataset)), "\n\nnumber of samples:", 
            nlevels(fr$sample_names), "\n\nsamples:", 
            list(levels(fr$sample_names)), "\n\nmethods:", 
            list(levels(fr$method)), "\n\n\n\n\n"), 
      file = paste(new_name, "_info.TXT", sep = ""), append = TRUE)

write.csv(fr, file = paste(new_name, "_ready.csv", sep = ""), 
          row.names = FALSE)

}

rm(fr, length, new_name, report_name, temp)
