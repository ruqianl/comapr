#' cal_marker_dist
#' 
#' Covert the recombination rate to cM by mapping functions
#' @importFrom GenomicRanges mcols
#' @noRd
cal_marker_dist <- function(co_count,mapping_fun="k",by_group = NULL){
  

  ## just covert the recombination rate to cM by mapping functions

  new_gr <- co_count[,0]
  
  if(is.null(by_group)){ 
    co_rate <- rowMeans(as.matrix(GenomicRanges::mcols(co_count)))
    stopifnot(sum(co_rate>=0.5)==0)
    if(mapping_fun=="k"){
    kosambi_cm <- 25*log((1+2*co_rate)/(1-2*co_rate))
    GenomicRanges::mcols(new_gr)$kosambi_cm <- kosambi_cm
    }
    if(mapping_fun=="h"){
      haldane_cm <- -50 * log(1 - 2 * co_rate)
      GenomicRanges::mcols(new_gr)$haldane_cm <- haldane_cm
    }
  } else {
    ## by_group has the character of prefix of groups
     groups_rb <- bplapply(by_group,function(group_prefix,mapping_fun){
       sids <- grep(group_prefix,colnames(GenomicRanges::mcols(co_count)))
       co_rate <- rowMeans(as.matrix(GenomicRanges::mcols(co_count)[,sids]))
       
      if(mapping_fun=="k"){
        x <- 25*log((1+2*co_rate)/(1-2*co_rate))
      }
      if(mapping_fun=="h"){
       x <-  -50 * log(1 - 2 * co_rate)
      }
       x
    },mapping_fun=mapping_fun)
     
     GenomicRanges::mcols(new_gr) <- do.call(cbind,groups_rb)
      colnames(GenomicRanges::mcols(new_gr)) <- paste0(by_group,"_",mapping_fun)
      
  }
  new_gr
}

#' calGeneticDist
#' 
#' Calculate genetic distances of marker intervals or binned-chromosome 
#' Given whether crossover happens in each marker interval, calculate the
#' recombination fraction in samples and then derive the Haldane or Kosambi
#' genetic distances via mapping functions
#'
#' @param co_count
#' GRange object, returned by \code{countCO}
#' @param bin_size
#' The binning size for grouping marker intervals into bins. If not supplied,the
#' orginial marker intervals are returned with converted genetic distancens
#' based on recombination rate
#' @param mapping_fun
#' The mapping function to use, can be one of "k" or "h" (kosambi or haldane)
#' @param ref_genome
#' The reference genome name. It is used to fetch the chromosome Information
#' from UCSC database.
#' @param by_group, character vector contains the unique prefix of sample names 
#' that are used for defining different sample groups. Or the column name in 
#' colData(co_count) that specify the group factor. If missing all samples are 
#' assumed to be from one group
#' 
#' @importFrom rlang .data
#' @importFrom GenomeInfoDb getChromInfoFromUCSC genome genome<-   
#' @importFrom GenomicRanges GRanges mcolAsRleList sort width 
#' @importFrom GenomicRanges tileGenome binnedAverage mcols
#' @importFrom GenomeInfoDb sortSeqlevels seqlevels seqinfo
#' @importFrom S4Vectors Rle
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment colData rowRanges
#' @return GRanges object
#' GRanges for marker intervals or binned intervals with Haldane or Kosambi 
#' centiMorgans
#' @export
#' 

setGeneric("calGeneticDist",
           function(co_count,
                    bin_size=NULL,
                    mapping_fun="k",
                    ref_genome="mm10",
                    by_group = NULL) 
             standardGeneric("calGeneticDist"))

#'@noRd
calGenetic_dist <- function(co_count,
                            bin_size=NULL,
                            mapping_fun="k",
                            ref_genome="mm10",
                            by_group = NULL){
  stopifnot(is.null(bin_size) | is.numeric(bin_size))
  stopifnot(mapping_fun %in% c("k","h"))
  new_gr <- cal_marker_dist(co_count = co_count,
                          mapping_fun = mapping_fun,
                          by_group = by_group)
  if(is.null(bin_size)){
    return(new_gr)
  } else {
    binned_dna_mm10_gr <- cal_bin_dist(new_gr=new_gr,bin_size=bin_size,
                                       ref_genome="mm10")
    return(binned_dna_mm10_gr)
  }
  
}

#'@noRd
cal_bin_dist <- function(new_gr,bin_size,
                         ref_genome="mm10"){
  ## bin_size supplied then.
  ## fetch the chromoInfo from GenomeInfoDb. 
  ## This is only for getting the basepair lengths of the genome
  
  
  chrom_info <- GenomeInfoDb::getChromInfoFromUCSC(ref_genome)
  ## only for chr1-M
  chrom_info <- chrom_info[grep("_",chrom_info$chrom,invert = TRUE),]
  
  ## Check what seqnames is in new_gr and make it consistent
  if(!grepl("chr",as.character(seqnames(new_gr)[1]))){
    chrom_info$chrom <- gsub("chr","",chrom_info$chrom)
  }
  chrom_info <- chrom_info[chrom_info$chrom %in% GenomeInfoDb::seqlevels(new_gr),]
  ## create Granges object for chromosomes
  seq_length <- chrom_info$size
  names(seq_length) <- chrom_info$chrom
  
  dna_mm10_gr <- GenomicRanges::GRanges(
    seqnames = Rle(names(seq_length)),
    ranges = IRanges(1, end = seq_length, names = names(seq_length)),
    seqlengths = seq_length)
  GenomeInfoDb::genome(dna_mm10_gr) <- ref_genome
  #dna_mm10_gr
  
  
  ## per bp distances
  GenomicRanges::mcols(new_gr) <- apply(GenomicRanges::mcols(new_gr),2,
                                        function(x) x/GenomicRanges::width(new_gr))
  
  tilewidth <- bin_size
  tiles <- GenomicRanges::tileGenome(seqinfo(dna_mm10_gr),
                                     tilewidth = tilewidth)
  binned_dna_mm10_gr <- unlist(tiles)
  # binned_dna_mm10_gr
  new_gr <- GenomicRanges::sort(GenomeInfoDb::sortSeqlevels(new_gr))
  
  bin_dist <-  bplapply(colnames(mcols(new_gr)), function(group_col){
    # dist_rle <- GenomicRanges::mcolAsRleList(new_gr,group_col)
    # runValue(dist_rle)[is.na(runValue(dist_rle))] <- 0
    dist_rle <- GenomicRanges::coverage(new_gr,weight = mcols(new_gr)[,group_col])
    dist_bined <- binnedAverage(binned_dna_mm10_gr,dist_rle,
                                "dist_bin_ave")
    
    return(dist_bined$dist_bin_ave*width(dist_bined))
    
  })
  
  mcols(binned_dna_mm10_gr) <- do.call(cbind,bin_dist)
  colnames(mcols(binned_dna_mm10_gr)) <- colnames(mcols(new_gr))
  binned_dna_mm10_gr
}

#'@rdname calGeneticDist
setMethod("calGeneticDist",signature = c(co_count = 'GRanges',
                                         bin_size='missing',
                                         by_group='missing'),
          
          function(co_count,
                   bin_size=NULL,
                   mapping_fun="k",
                   ref_genome="mm10",
                   by_group = NULL){
            stopifnot(mapping_fun %in% c("k","h"))
            
            new_gr <- cal_marker_dist(co_count = co_count,
                                      mapping_fun = mapping_fun,
                                      by_group = NULL)
            new_gr
})
#'@rdname calGeneticDist
setMethod("calGeneticDist",signature = c(co_count = 'GRanges',
                                         bin_size='numeric',
                                         by_group='missing'),
          
          function(co_count,
                   bin_size=NULL,
                   mapping_fun="k",
                   ref_genome="mm10",
                   by_group = NULL){
            stopifnot(mapping_fun %in% c("k","h"))
            
            new_gr <- calGenetic_dist(co_count = co_count,
                                      mapping_fun = mapping_fun,
                                      bin_size = bin_size,
                                      by_group = NULL)
            new_gr
          })

#'@rdname calGeneticDist
setMethod("calGeneticDist",signature = c(co_count = 'GRanges',
                                         bin_size='missing',
                                         by_group='character'),
          
          function(co_count,
                   bin_size=NULL,
                   mapping_fun="k",
                   ref_genome="mm10",
                   by_group  ){
            stopifnot(mapping_fun %in% c("k","h"))
            
            new_gr <- cal_marker_dist(co_count = co_count,
                                      mapping_fun = mapping_fun,
                                      by_group = by_group)
            new_gr
})

#'@rdname calGeneticDist
setMethod("calGeneticDist",signature = c(co_count = 'GRanges',
                                         bin_size='numeric',
                                         by_group='character'),
          
          function(co_count,
                   bin_size,
                   mapping_fun="k",
                   ref_genome="mm10",
                   by_group){
            stopifnot(mapping_fun %in% c("k","h"))
            
            new_gr <- calGenetic_dist(co_count = co_count,
                                      mapping_fun = mapping_fun,
                                      bin_size = bin_size,
                                      by_group = by_group)
            new_gr
          })

#'@rdname calGeneticDist
#'@importFrom SummarizedExperiment rowRanges rowRanges<-
#'
setMethod("calGeneticDist",signature = c(co_count = 'RangedSummarizedExperiment',
                                         bin_size='missing',
                                         by_group='missing'),
          
          function(co_count,
                   bin_size=NULL,
                   mapping_fun="k",
                   ref_genome="mm10",
                   by_group = NULL){
            stopifnot(mapping_fun %in% c("k","h"))
            co_rate <- rowMeans(as.matrix(assay(co_count)))
            rowRanges(co_count)$raw_rate <- co_rate
            if(mapping_fun == "k")
              rowRanges(co_count)$kosambi <- 25*log((1+2*co_rate)/(1-2*co_rate))
            else
              rowRanges(co_count)$haldane <-  -50 * log(1 - 2 * co_rate)
            
            co_count
          })

#'@rdname calGeneticDist
setMethod("calGeneticDist",
          signature = c(co_count = 'RangedSummarizedExperiment',
                        bin_size='missing',
                        by_group='character'),
          
          function(co_count,
                   bin_size=NULL,
                   mapping_fun="k",
                   ref_genome="mm10",
                   by_group ){
            stopifnot(mapping_fun %in% c("k","h"))
            stopifnot(length(by_group)==1)
            ##by_group should refer to a column in colData(co_count)
            stopifnot(by_group %in% colnames(colData(co_count)))
            fcts <- colData(co_count)[,by_group]
            
            co_rate <- sapply(unique(as.character(fcts)),
                              function(fct){
                                rowMeans(as.matrix(assay(co_count)[,fcts==fct]))})
            
            rowRanges(co_count)$raw_rate <- co_rate
            if(mapping_fun == "k")
              rowRanges(co_count)$kosambi <- 25*log((1+2*co_rate)/(1-2*co_rate))
            else
              rowRanges(co_count)$haldane <-  -50 * log(1 - 2 * co_rate)
            
            co_count
          })


#'@rdname calGeneticDist
setMethod("calGeneticDist",
          signature = c(co_count = 'RangedSummarizedExperiment',
                        bin_size='numeric',
                        by_group='character'),
          
          function(co_count,
                   bin_size,
                   mapping_fun="k",
                   ref_genome="mm10",
                   by_group ){
            
            stopifnot(mapping_fun %in% c("k","h"))
            stopifnot(length(by_group)==1)
            ##by_group should refer to a column in colData(co_count)
            stopifnot(by_group %in% colnames(colData(co_count)))
            fcts <- colData(co_count)[,by_group]
            
            #####
            # divide the cos into bins
            ####
            co_rate <- sapply(unique(as.character(fcts)),
                              function(fct){
                                rowMeans(as.matrix(assay(co_count)[,fcts==fct]))})
            
            rowRanges(co_count)$raw_rate <- co_rate 
            ## divide into bins

            if(mapping_fun == "k"){
              new_gr <- granges(co_count)
              kosambi <- 25*log((1+2*co_rate)/(1-2*co_rate))
              mcols(new_gr) <- kosambi
              dist_gr <- cal_bin_dist(new_gr,bin_size = bin_size)
              
              } else {
              
              new_gr <- granges(co_count)
              haldane <-  -50 * log(1 - 2 * co_rate)
              mcols(new_gr) <- haldane
              dist_gr <- cal_bin_dist(new_gr,bin_size = bin_size)
            }
            dist_gr
            
          })



#'@rdname calGeneticDist
setMethod("calGeneticDist",
          signature = c(co_count = 'RangedSummarizedExperiment',
                        bin_size='numeric',
                        by_group='missing'),
          
          function(co_count,
                   bin_size,
                   mapping_fun="k",
                   ref_genome="mm10",
                   by_group =NULL){
            
            stopifnot(mapping_fun %in% c("k","h"))

            #####
            # divide the cos into bins
            ####
            co_rate <- rowMeans(as.matrix(assay(co_count)))
            
            rowRanges(co_count)$raw_rate <- co_rate 
            ## divide into bins
            
            if(mapping_fun == "k"){
              new_gr <- granges(co_count)
              kosambi <- 25*log((1+2*co_rate)/(1-2*co_rate))
              mcols(new_gr) <- kosambi
              dist_gr <- cal_bin_dist(new_gr,bin_size = bin_size)
              
            } else {
              new_gr <- granges(co_count)
              haldane <-  -50 * log(1 - 2 * co_rate)
              mcols(new_gr) <- haldane
              dist_gr <- cal_bin_dist(new_gr,bin_size = bin_size)
            }
          dist_gr
          
        })

