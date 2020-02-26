
#### run data ####

sourcePath <- "~/Documents/myprojects/replicationRepair/1_code/R/"

date <- "[2020.02.10]"

fr_name_list <- c(
                  # "inZones_windows_201_100", 
                  # "inZones_windows_201_1000", 
                  # "inZones_windows_201_10000",
                  # "hiRFD_windows_201_100",
                  # "hiRFD_windows_201_1000",
                  # "hiRFD_windows_201_10000",
                  # "sns_seq_windows_201_100",
                  # "sns_seq_windows_201_1000",
                  # "S1toS2_windows_201_100",
                  # "S1toS2_windows_201_1000",
                  "repdomains_windows_201_100",
                  "repdomains_windows_201_1000",
                  "repdomains_windows_201_10000"
                  # "repdomains",
                  # "genome_1Mb",
                  # "genome_100kb",
                  # "genome_20kb",
                  # "_repdomains_windows",
                  # "chromhmm_windows",
                  # "chromhmm_windows_chr"
                  )

for (fr_name in fr_name_list) {
  
  if (grepl("inZones_windows", fr_name)) { 
    
    source(paste(sourcePath, "3_initiation_zones_windows.R", sep = "")) 
    
  } else if (grepl("hiRFD_windows", fr_name)) { 
      
    source(paste(sourcePath, "3_highRFD_windows.R", sep = "")) 
    
  } else if (grepl("sns_seq_windows", fr_name)) { 

    source(paste(sourcePath, "3_sns_seq_rep1_rep2_windows.R", sep = ""))
  
  } else if (grepl("S1toS2_windows", fr_name)) {

    source(paste(sourcePath, "3_S1toS2_peak_windows.R", sep = ""))
      
  } else if (fr_name == "repdomains") { 
    
    source(paste(sourcePath, "3_repdomains.R", sep = ""))
    
  } else if (grepl("^repdomains_windows", fr_name)) { 
    
    source(paste(sourcePath, "3_repdomains_windows.R", sep = ""))

  } else if (grepl("_repdomains_windows", fr_name)) { 
    
    source(paste(sourcePath, "3_repdomains_intersect.R", sep = ""))
    
  } else if (grepl("chromhmm", fr_name)) { 
    
    source(paste(sourcePath, "3_chromhmm_chromatin_states.R", sep = ""))

  } else if (grepl("genome", fr_name)) {

    source(paste(sourcePath, "3_genome.R", sep = ""))
    }
}

