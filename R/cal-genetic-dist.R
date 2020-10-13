#' cal_marker_dist
#' @importFrom GenomicRanges mcols
#' @noRd
cal_marker_dist <- function(co_geno_gr,mapping_fun="k",by_group = NULL){
  

  ## just covert the recombination rate to cM by mapping functions

  new_gr <- co_geno_gr[,0]
  
  if(is.null(by_group)){ 
    co_rate <- rowMeans(as.matrix(GenomicRanges::mcols(co_geno_gr)))
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
       sids <- grep(group_prefix,colnames(GenomicRanges::mcols(co_geno_gr)))
       co_rate <- rowMeans(as.matrix(GenomicRanges::mcols(co_geno_gr)[,sids]))
       
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
#' Calculate genetic distances of marker intervals or bins 
#'
#' Given whether crossover happens in each marker interval, calculate the
#' recombination fraction in samples and then derive the Haldane or Kosambi
#' genetic distances via mapping functions
#'
#' @param co_geno_gr
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
#' @param by_group, the unique prefix of sample names that are used for defining
#' different sample groups. If NULL all samples are assumed to be from one group
#' @param BPPARAM, the back-end of type bpparamClass
#' @importFrom rlang .data
#' @importFrom GenomeInfoDb fetchExtendedChromInfoFromUCSC genome genome<-   
#' @importFrom GenomicRanges GRanges mcolAsRleList sort width 
#' @importFrom GenomicRanges tileGenome binnedAverage mcols
#' @importFrom GenomeInfoDb sortSeqlevels seqlevels seqinfo
#' @importFrom S4Vectors Rle
#' @return GRanges object
#' GRanges for marker intervals or binned intervals with Haldane or Kosambi 
#' centiMorgans
#' @export
calGeneticDist <- function(co_geno_gr,bin_size=NULL,mapping_fun="k",
                           ref_genome="mm10",by_group = NULL,
                           BPPARAM=BiocParallel::bpparam()){
  stopifnot(is.null(bin_size) | is.numeric(bin_size))
  stopifnot(mapping_fun %in% c("k","h"))
  new_gr <- cal_marker_dist(co_geno_gr = co_geno_gr,
                          mapping_fun = mapping_fun,
                          by_group = by_group)
  if(is.null(bin_size)){
    return(new_gr)
  } else {
    ## bin_size supplied then.
    ## fetch the chromoInfo from GenomeInfoDb. 
    ## This is only for getting the basepair lengths of the genome
    
    chrom_info <- GenomeInfoDb::fetchExtendedChromInfoFromUCSC(ref_genome)
    ## only for chr1-M
    chrom_info <- chrom_info[grep("_",chrom_info$UCSC_seqlevel,invert = TRUE),]

    ## Check what seqnames is in new_gr and make it consistent
    if(!grepl("chr",as.character(seqnames(new_gr)[1]))){
      chrom_info$UCSC_seqlevel <- gsub("chr","",chrom_info$UCSC_seqlevel)
    }
    chrom_info <- chrom_info[chrom_info$UCSC_seqlevel %in% GenomeInfoDb::seqlevels(new_gr),]
    ## create Granges object for chromosomes
    seq_length <- chrom_info$UCSC_seqlength
    names(seq_length) <- chrom_info$UCSC_seqlevel

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
    tiles <- GenomicRanges::tileGenome(seqinfo(dna_mm10_gr),tilewidth = tilewidth)
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
    return(binned_dna_mm10_gr)
  }
  
}

