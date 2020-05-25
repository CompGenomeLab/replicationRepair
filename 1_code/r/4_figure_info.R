
figureInfoPath <- sourcePath

# dataInfo, figureInfo, figureName, and figurePath will be delivered by the 
# main plot scripts (3_*). dataInfoPath will be delivered by the 
# 2_report_sub_dfs.R.

figureFile <- sub(".pdf", ".txt", figureName)

if (file.exists(paste(dataInfoPath, dataInfo, sep = ""))){
  
  print(paste(dataInfo, "exist!"))
  
  file.copy(from = file.path(dataInfoPath, dataInfo), 
            to = file.path(figurePath, figureFile), 
            overwrite = TRUE)
  
  file.append(file.path(figurePath, figureFile), 
              file.path(figureInfoPath, figureInfo))

} else {
  
  print(paste(dataInfo, "don't exist"))
  
  file.copy(from = file.path(figureInfoPath, figureInfo),
            to = file.path(figurePath, figureFile), 
            overwrite = TRUE)
}
