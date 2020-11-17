extract_hugo_symbol_colnames <- function(x){
  
  stopifnot(is.data.frame(x) | is.matrix(x))
  
  if (is.data.frame(x) | is.matrix(x)) {

    
    if (any(!is.na(stringr::str_match(colnames(x), " \\([0-9\\&]+\\)$")[, 1]))) {
      colnames(x) <- str_match(colnames(x), "^(.+) \\([0-9\\&]+\\)$")[, 2]
    }
    else if (any(!is.na(stringr::str_match(colnames(x), " \\(ENSG[0-9\\.]+\\)$")[, 1]))) {
      colnames(x) <- str_match(colnames(x), "^(.+) \\(ENSG[0-9\\.]+\\)$")[, 2]
    }
    x <- x[, !is.na(colnames(x)) & colnames(x) != ""]
    
  }
  return(x)
}