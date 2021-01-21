#' combineHapState
#' 
#' combine two `RangedSummarizedExperiment` objects, each contains the haplotype
#' state for a list of SNPs across a set of cells. The combined result will have
#' cells from two individuals and merged list of SNPs from the two.
#' 
#' @importFrom  Matrix Matrix
#' @importFrom IRanges findOverlaps
#' @importFrom SummarizedExperiment rowRanges assay assay<- rowRanges<- assays assays<- 
#' @importFrom S4Vectors merge
#' @param rse1, the first `RangedSummarizedExperiment`
#' @param rse2, the second `RangedSummarizedExperiment`
#' @param groupName, a character vector of length 2 that contains the first and 
#' the second group's names
#' @examples
#' demo_path <-paste0(system.file("extdata",package = "comapr"),"/")
#' s1_rse_state <- readHapState("s1",chroms=c("chr1"),
#'                              path=demo_path,barcodeFile=NULL,minSNP = 0,
#'                              minlogllRatio = 50,
#'                              bpDist = 100,maxRawCO=10)
#' 
#' s2_rse_state <- readHapState("s2",chroms=c("chr1"),
#'                              path=demo_path,
#'                              barcodeFile=paste0(demo_path,"s2_barcodes.txt"),minSNP = 0, 
#'                              minlogllRatio = 50,
#'                              bpDist = 100,maxRawCO=10)
#' sb <- combineHapState(s1_rse_state,s2_rse_state)
#' @return A `RangedSummarizedExperiment` that contains the cells and SNPs in 
#' both `rse`
#' 
#' @export
#' 
#' @author Ruqian Lyu
#' 
combineHapState <- function(rse1,rse2,
                             groupName=c("Sample1","Sample2")){

  stopifnot(length(groupName)==2)
  
  merged_anno <- S4Vectors::merge(rowRanges(rse1),rowRanges(rse2),all=T)
  
  stopifnot(!is.data.frame(merged_anno))
  t <- Matrix(data=0,nrow=length(merged_anno),ncol=(ncol(rse1)+ncol(rse2)),
              sparse = T)
  se1_rows <- IRanges::findOverlaps(merged_anno,rowRanges(rse1),select = "all")
  se1_rows <- se1_rows@from
  t[se1_rows,1:ncol(rse1)] <- SummarizedExperiment::assays(rse1)[["vi_state"]]
  se2_rows <- findOverlaps(merged_anno,rowRanges(rse2),select = "all")
  se2_rows <- se2_rows@from
  t[se2_rows,(ncol(rse1)+1):ncol(t)] <- assay(rse2)
  combinedCols <- rbind(colData(rse1),colData(rse2))
  
  ## should combine meta data as well
  
  SummarizedExperiment(assay=list(vi_state=t),
                         colData = combinedCols,
                         rowRanges = merged_anno) 
  
 

}