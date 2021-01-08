Beckettnested <- function(web, result){
  # a function to be called when a recursive modularisation is desired
  # result@modules # 1st and 2nd column can be ignored; other values give the module of which this species is a member
  # species are in columns, with row-species first (e.g. plants)
  
  # 1st ROW is the original web
  # 2nd row is the first module etc
  
  nr <- NROW(web)
  nc <- NCOL(web)
  
  nr.modules <- nrow(result@modules) - 1
  result.nested <- list()
  
  for (i in seq_len(nr.modules)){
    subm <- which(result@modules[i+1, -c(1:2)] != 0) # identify module members
    subweb.rows <- subm[subm <= nr]
    subweb.cols <- subm[subm > nr & subm <= (nr +nc)]
    subweb <- web[subweb.rows, subweb.cols - nr]
    
    # re-apply LPA:
    result.nested[[i]] <- computeModules(subweb)
  }
  
}