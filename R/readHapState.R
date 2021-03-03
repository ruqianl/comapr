#' readHapState
#' 
#' A helper function that parses the viterbi state matrix (in .mtx format), 
#' barcode.txt and snpAnno.txt files for each individual.
#' 
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment colData assay assay<- SummarizedExperiment
#' @importFrom SummarizedExperiment rowRanges metadata<- assays assays<-
#' @importFrom Matrix readMM
#' @importFrom utils read.table

#' @param path, the path to the files, with name patterns *{chrom}_vi.mtx, 
#' *{chrom}_viSegInfo.txt, end with slash
#' @param barcodeFile, if NULL, it is assumed to be in the same directory as the 
#' other files and with name sampleName_barcodes.txt
#' @param chroms, the character vectors of chromosomes to parse. Multiple chromosomes'
#' results will be concated together.
#' @param sampleName, the name of the sample to parse which is used as prefix for
#' finding relevant files for the underlying sample
#' @inheritParams .filterCOsExtra
#' @return a RangedSummarizedExperiment with rowRanges as SNP positions that contribute
#' to crossovers in any cells. colData contains cells annotation including barcodes and
#' sampleName.
#' 
#' @examples
#' demo_path <-system.file("extdata",package = "comapr")
#' s1_rse_state <- readHapState(sampleName="s1",chroms=c("chr1"),
#' path=paste0(demo_path,"/"),
#' barcodeFile=NULL,minSNP = 0, minlogllRatio = 50,
#' bpDist = 100,maxRawCO=10,minCellSNP=10)
#' s1_rse_state
#' @export
#' @author Ruqian Lyu
#' 

readHapState <- function(sampleName,chroms=c("chr1"),path,
                         barcodeFile=NULL,minSNP = 30, minlogllRatio = 200,
                         bpDist = 100,maxRawCO=10,nmad=1.5,minCellSNP = 200,
                         biasTol=0.45){
  if(is.null(barcodeFile)){
    barcodeFile <- paste0(path,sampleName,"_barcodes.txt")
  }
  barcodes <- read.table(file=barcodeFile,stringsAsFactors = F,
                           col.names = "barcodes")
  se_list <- bplapply(chroms,function(chr){

    snpAnno <- read.table(file=paste0(path,sampleName,"_",chr,"_snpAnnot.txt"),
                          stringsAsFactors = F,
                          header=T)
    segInfo <- read.table(file = paste0(path,sampleName,"_",chr,"_viSegInfo.txt"),
                          stringsAsFactors = F,
                          col.names = c("ithSperm","Seg_start","Seg_end",
                                        "logllRatio","nSNP","State"))
    vi_mtx <- readMM(file = paste0(path,sampleName,"_",chr,"_vi.mtx"))
    
    se <- SummarizedExperiment::SummarizedExperiment(assays = 
                                                       list(vi_state = vi_mtx),
                               colData = barcodes,
                               rowRanges = GRanges(
                                 seqnames = chr,
                                 ranges = IRanges::IRanges(start = snpAnno$POS,
                                                           width = 1)),
                               metadata = data.frame(segInfo,chr=chr,
                                                     sampleGroup=sampleName))

    se <- .filterCOsExtra(se ,minSNP = minSNP, 
                         minlogllRatio = minlogllRatio,
                         minCellSNP = minCellSNP,
                         bpDist = bpDist,
                         maxRawCO=maxRawCO,
                         biasTol=biasTol,
                         nmad=nmad)
    se
  })
  
  rbind_se <- function(se1,se2){
    cc <- intersect(colnames(se1),colnames(se2))
    ## rbind does not combined the metadata, so
    rbind(se1[,cc],se2[,cc])
  }
  suppressWarnings(Reduce(rbind_se,se_list))
  #do.call(rbind,se_list)
}

#'Filter out doublet cells and uninformative SNPs
#'
#'This function filter out cells that have been called too many crossovers due
#'to diploid cell contamination or doublets. It also only keeps SNPs (rows) that
#'ever contribute to a crossover interval. This function should be run for 
#'individual chromosomes and is called internaly by `readHapState`

#'@param se, the SummarizedExperiment object that contains the called haplotype
#'state matrix in the assay field and haplotype segment information in the metadata
#'field.
#'@param maxRawCO, if a cell has more than `maxRawCO` number of raw crossovers
#'called across one chromosome, the cell is filtered out
#'@param minSNP, the crossover(s) will be filtered out if introduced by a segment
#'that has fewer than `minSNP` SNPs to support.
#'@param bpDist, the crossover(s) will be filtered out if introduced by a segment
#'that is shorter than `bpDist` basepairs.
#'@param minlogllRatio, the crossover(s) will be filtered out if introduced by a
#'segment that has lower than `minlogllRatio` to its reversed state.
#'@param minCellSNP, the minimum number of SNPs detected for a cell to be kept, 
#'used with `nmads`
#'@param biasTol, the SNP's haplotype ratio across all cells is assumed 
#'to be 1:1. This argument can be used for removing SNPs that have a biased
#'haplotype. i.e. almost always inferred to be haplotype state 1. It specifies
#'a bias tolerance value, SNPs with bias from 0.5 smaller than this value are 
#'kept.
#'@param nmad, how many mean absolute deviations lower than the median number 
#'of SNPs per cellfor a cell to be considered as low coverage cell and filtered 
#' This or `minCellSNP`, whichever is larger, is applied
#' 
#'@importFrom Matrix Matrix
#'@importFrom SummarizedExperiment metadata<- assay<- colData rowRanges rowRanges<-
#'@importFrom S4Vectors metadata
#'@importFrom IRanges ranges
#'@importFrom stats mad median
#'@return A `RangedSummarizedExperment` object that have different dims with input.
#'the colnames are the cell barcodes, rowRanges specify the location of SNPs that
#'contribute to crossovers.
#'@details The `logllRatio` value is returned by `sscocaller` for each haplotype
#'segment formed by consecutive SNPs that are called to have a same state. It is
#'calculated by taking log of ratio (likelihood of SNPs with inferred states) 
#'and (likelihood of SNPs with reversed states)
#'
#'@author Ruqian Lyu
#'@keywords internal

.filterCOsExtra <- function(se,minSNP = 30, minlogllRatio = 200,
                            minCellSNP = 200,
                            bpDist = 100,
                            maxRawCO = 10,
                            biasTol = 0.45,
                            nmad = 1.5){
  total_SNP <- total_co <- barcode <- nSNP <- NULL
  segInfo <- data.frame(S4Vectors::metadata(se),stringsAsFactors = F)
  
  segInfo$bp_dist <- segInfo$Seg_end - segInfo$Seg_start
  # ithSperm counts from 0. thus +1
  segInfo$barcode <- SummarizedExperiment::colData(se)$barcodes[as.numeric(gsub("ithSperm",
                                                                                "",segInfo$ithSperm))+1]
  
  #keepCells <- names(table(segInfo$barcode))[table(segInfo$barcode)-1 <= maxRawCO]
  
  suppressMessages(
    keepCells <- segInfo %>% dplyr::group_by(barcode) %>% 
                     dplyr::summarise(total_SNP = sum(nSNP),
                                      total_co = length(nSNP)-1))
  ## provided minCellSNP or nmads smaller from median, whichever is larger
  minCellSNP <- ifelse((median(keepCells$total_SNP)-nmad*mad(keepCells$total_SNP))<minCellSNP,
                       minCellSNP,(median(keepCells$total_SNP)-nmad*mad(keepCells$total_SNP)))
  
  keepCells <- keepCells %>%
    dplyr::filter(total_co <= maxRawCO & 
                    total_SNP > minCellSNP )
  
  keepCells <-  keepCells$barcode
  
  segInfo_f <- segInfo[segInfo$barcode %in% keepCells
                       & segInfo$nSNP > minSNP 
                       & segInfo$logllRatio > minlogllRatio 
                       & segInfo$bp_dist > bpDist,]
  
  ## get Crossover contributing SNPs only
  co_contr_snp <- unique(c(segInfo_f$Seg_end,segInfo_f$Seg_start))
  co_contr_snp <- co_contr_snp[order(co_contr_snp)]
  
  ## remove some SNPs if their inferred states are biased
  markers <- co_contr_snp 
  queryR <-  IRanges::IRanges(start = markers,
                     width = 1)
  
  markers_states <- lapply(unique(segInfo_f$ithSperm),function(ithS){
    #segInfo_f[segInfo_f$ithSperm==ithS,]
    subjectR <- IRanges(start = segInfo_f[segInfo_f$ithSperm==ithS,]$Seg_start,
                        end = segInfo_f[segInfo_f$ithSperm==ithS,]$Seg_end,
                        state = segInfo_f[segInfo_f$ithSperm==ithS,]$State)
    
    temp <- rep(0,length(markers))
    myHits <- IRanges::findOverlaps(queryR,subjectR)
    temp[myHits@from] <- subjectR@elementMetadata@listData$state[myHits@to]  
    temp
  })
  
  markers_states <- do.call(cbind,markers_states)
  
  ratio <- rowSums(markers_states==1)/(rowSums(markers_states==2) + 
                                         rowSums(markers_states==1))
  
  #hist(ratio)
  
  co_contr_snp <- markers[abs(ratio-0.5)<biasTol]
  
  #  markers[which(ratio==0)]
  ## row numbers of these SNPs in the assay
  ithSNP <- match(co_contr_snp,IRanges::ranges(rowRanges(se ))@start)
  # 
  ## subset the columns to only keep the Good cells
  colnames(se) <- colData(se)$barcodes
  
  se <- se[ithSNP,keepCells]
  
  vi_m <- Matrix(data=0,nrow=nrow(se),ncol=ncol(se))
  queryR <- IRanges::IRanges(start=co_contr_snp,
                             width=1)
  
  ##<<<<<----UPDATE,TEST--------->>>>
  vi_m <- sapply(colnames(se),function(bc){
    nthCol <- which(colnames(se)==bc)
    #  message(nthCol)
    segs <- segInfo_f[segInfo_f$barcode==bc,]
    
    # subjectR <- IRanges(start =segs[,"Seg_start"],
    #                     end = segs[,"Seg_end"],
    #                     state = segs[,"State"])
    # 
    # temp <- rep(0,length(co_contr_snp))
    # myHits <- IRanges::findOverlaps(queryR,subjectR)
    # temp[myHits@from] <- subjectR@elementMetadata@listData$state[myHits@to]  
    # temp
    
    # 
    # 
    # posS <- as.numeric(segs[,"Seg_start"])
    # posE <- as.numeric(segs[,"Seg_end"])
    # # 
    # mthrows <- match(c(posS,posE),co_contr_snp)
    temp <- rep(0,length(co_contr_snp))
    
    mthrows <- match(co_contr_snp,as.numeric(segs[,"Seg_start"]))
    
    temp[!is.na(mthrows)] <- segs[,"State"][mthrows[!is.na(mthrows)]]
    
    mthrows <- match(co_contr_snp,as.numeric(segs[,"Seg_end"]))
    
    temp[!is.na(mthrows)] <- segs[,"State"][mthrows[!is.na(mthrows)]]
    
    temp
  })
  vi_m <- Matrix(vi_m,sparse = T)
  SummarizedExperiment::assays(se)[["vi_state"]] <- vi_m 
  SummarizedExperiment::metadata(se) <- segInfo_f
  se
}
