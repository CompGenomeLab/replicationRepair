function (x) 
{
  if (is.list(x)) {
    if (length(x[[1]]) == 1 && x[[1]] %in% c("a", "A", "1", 
                                             "i", "I", "fig2", "fig4", 
                                             "supp_rr", "supp_muts",
                                             "fig1", "supp2", "supfig11")) {
      x <- x[[1]]
    }
    else {
      return(x[[1]])
    }
  }
  switch(as.character(x), a = letters, A = LETTERS, `1` = 1:100, 
         fig2 = c("A", "", "B"), i = tolower(as.roman(1:100)), 
         fig4 = c("A", "B", "", "", "", "C", "", "", ""),
         supfig11 = c("A", "", "B", "", "", "", "C", "", "", ""),
         supp_rr = c("A", "", "", "", "B", "", "", ""),
         supp_muts = c("A", "", "", "B", "", "","C", "", ""),
         fig1 = c("A", "B", "", "C", "", "","", "D"),
         supp2 = c("A", "B", "C", "", "D"),
         I = as.roman(1:100), stop("Unknown tag type: ", x, call. = FALSE))
}

library(patchwork)
trace("get_level", edit = TRUE, where = asNamespace("patchwork"))

