#### libraries ####

library(stringr)


#### functions ####

output_date <- function (){
  dateout <- Sys.Date()
  dateout <- format(dateout, format = "[%Y.%m.%d]")
  return(dateout)
}


get_sub_str <- function ( whole_str, substr1, substr2 ){
  target_str <- str_match(whole_str, paste(substr1, "(.*?)", substr2, sep = ""))
  return(target_str[2])
}


middle <- function( window_number ){
  if (window_number %% 2 == 0) {
    middle <- window_number / 2
  } else {
    middle <- (window_number - 1) / 2 + 1
  }
  return(middle)
}

region_length <- function( window_number, window_length ){
  rlength <- (window_number - 1) * window_length / 2000 # in kb and half of it. 
  return(rlength)
}

window_numbering <- function( mydata, name_column, middle ){
  window <- cbind(mydata[ , name_column], 
                  data.frame(str_split_fixed(mydata[ , name_column], "_", -1)))
  window <- window[ , c(1, ncol(window))]
  names(window) <- c("dataset", "windows")
  window$windows <- as.numeric(as.character(window$windows)) - middle
  newdata <- merge(mydata, window, by.x="dataset", by.y="dataset")
  return(newdata)
}

domain_name <- function( mydata, name_column ){
  window <- cbind(mydata[ , name_column], 
                  data.frame(str_split_fixed(mydata[ , name_column], "_", -1)))
  window <- window[ , c(1, ncol(window)-1)]
  names(window) <- c("dataset", "repdomains")
  mydata$repdomains <- window[,2]
  mydata$repdomains <- factor(mydata$repdomains, 
                              levels = c("UTZ", "ERD", "DTZ", "LRD"))
  return(mydata)
}

fig_save <- function( figpath, figname ) {
  ggsave(path = figurePath, filename = figureName, 
         width = 297, height = 210, units = "mm")
  figurePNG <- sub(".pdf", ".png", figureName)
  ggsave(path = figurePath, filename = figurePNG, 
         width = 297, height = 210, units = "mm")
}

