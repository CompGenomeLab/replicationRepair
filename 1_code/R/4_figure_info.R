
dataInfoPath <- paste("~/Documents/My_Projects/Project_Repair_Replication/", 
                      "Results/1_TextforPlotting", sep = "")

figureInfoPath <- paste("~/Documents/My_Projects/Project_Repair_Replication/", 
                        "Scripts/R", sep = "") 

# dataInfo, figureInfo, figureName, and figurePath will be delivered by the 
# main plot scripts (3_*).

figureFile <- sub(".pdf", ".txt", figureName)

file.copy(from = file.path(dataInfoPath, dataInfo), 
          to = file.path(figurePath, figureFile), 
          overwrite = TRUE)

file.append(file.path(figurePath, figureFile), 
            file.path(figureInfoPath, figureInfo))


