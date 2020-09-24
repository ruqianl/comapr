#' getDistortedMarkers
#'
#' Marker disortation detection using chisq-test
#'
#' @details
#' We expect the genotypes to appear with the frequenceis of 1:1 homo_alt:hets.
#' We use chisq-test for finding markers that have genotypes among samples that
#' are significantly diffferent from the 1:1 ratio and report them
#'
#' @param p the expected geneotype ratio in a numeric vector, defauls to c(0.5,0.5)
#' @importFrom stats chisq.test p.adjust
#' @inheritParams countGT
#' @author Ruqian Lyu
#'
#' @export

getDistortedMarkers <- function(geno, p = c(0.5,0.5)){
  geno_table <- sapply(rownames(geno), function(marker){
    list(Het = ifelse(is.na(table(geno[marker,],useNA = "no")["Het"]),
                      0,
                      table(geno[marker,],useNA = "no")["Het"]),
         Homo_alt = ifelse(is.na(table(geno[marker,],useNA = "no")["Homo_alt"]),
                           0,
                           table(geno[marker,],useNA = "no")["Homo_alt"]))
  })
  geno_table <- data.frame(Markers = colnames(geno_table),
                           No.Hets =  as.character(unlist(geno_table["Het",])),
                           No.Homo_alt = as.character(unlist(geno_table["Homo_alt",])),
                           stringsAsFactors = FALSE)
  
  pvals <- sapply(as.character(geno_table$Markers), function(marker){
    ctest <- chisq.test(as.numeric(geno_table[geno_table$Markers==marker,2:3]),
                        p = p)
    ctest$p.value
  })
  
  #names(pvals) == geno_table$Markers
  geno_table$Pvals <- pvals
  geno_table$Adj.pvals <- p.adjust(pvals, method = "BH")
  
  return(geno_table[order(geno_table$Adj.pvals),])
}
