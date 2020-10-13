#' calculate observed double crossover rate for paired intervals
#' @noRd
#' @importFrom GenomeInfoDb fetchExtendedChromInfoFromUCSC genome seqinfo  
#' @importFrom GenomicRanges GRanges mcolAsRleList sort width 
#' @importFrom GenomicRanges tileGenome binnedAverage mcols mcols<-
#' @importFrom GenomeInfoDb sortSeqlevels seqlevels
#' @importFrom S4Vectors Rle
#' @importFrom utils combn
#' @keywords internal
#' @author Ruqian Lyu

cal_coc_by_group <- function(binned_dna_mm10_gr, group_prefix = NULL){
  chrs <- unique(as.character(seqnames(binned_dna_mm10_gr)))
  all_chr_id <- bplapply(chrs,function(chr){
    if(is.null(group_prefix)){
      binned_dna_mm10_gr_b <- binned_dna_mm10_gr[as.character(seqnames(binned_dna_mm10_gr))==chr,]
    } else {
      binned_dna_mm10_gr_b <- 
        binned_dna_mm10_gr[as.character(seqnames(binned_dna_mm10_gr))==chr,
                           grep(group_prefix,colnames(mcols(binned_dna_mm10_gr)))]
    }
    
    
    rbs <- rowMeans(as.matrix(mcols(binned_dna_mm10_gr_b)))
    mcols(binned_dna_mm10_gr_b) <- as.matrix(mcols(binned_dna_mm10_gr_b)) > 0.5
    #mcols(binned_dna_mm10_gr_b)$rbs <- rbs
    df_cos <- data.frame(mcols(binned_dna_mm10_gr_b),check.names=F)
    
    n <- seq_len( nrow(df_cos))
    
    #  Make combinations
    id <- data.frame(t(combn(nrow(df_cos),2)))
    id$inter_dist <- abs(id$X1 - id$X2)/(nrow(df_cos)-1)
    #table(id$inter_dist)
    #  Get result
    
    dbl_co_two_interval <- apply(id,1,function(row_id){
      colSums(df_cos[c(row_id[1],row_id[2]),])==2
    })
    
    id$observe_co_rate <- colMeans(dbl_co_two_interval)
    id$j_dbl_co_rate <- rbs[id$X1]*rbs[id$X2]
    id$coc <- id$observe_co_rate/id$j_dbl_co_rate
    id$chr <- chr
    id$group_prefix <- group_prefix
    id
  })
  all_chr_id <- do.call(rbind,all_chr_id)
  all_chr_id
}

#' calculate the  coefficient of coincidence for grouped samples
#' 
#' It takes the crossover count genomic range object and distribute the crossovers
#' to `n_intervals` per chromosome. The intervals along each chromosome are then
#' paired and the observed double crossover rate is calculated. The coc is then
#' derived by taking ratio of observed_db_co_rate over expected_db_co_rate.
#' 
#' @param count_co_gr, the crossover count genomic ranges object
#' @param n_intervals, int, the number of intervals to subset each chromosome 
#' @param ref_genome, character, the reference genome
#' @param group_prefix, character vector, the prefix for grouping samples in 
#' sample names.
#' 
#' @importFrom GenomeInfoDb fetchExtendedChromInfoFromUCSC genome seqinfo  
#' @importFrom GenomicRanges GRanges mcolAsRleList sort width 
#' @importFrom GenomicRanges tileGenome binnedAverage mcols
#' @importFrom GenomeInfoDb sortSeqlevels seqlevels
#' @importFrom S4Vectors Rle
#' 
#' @return a data.frame with CoC calculated for interval pairs in each
#' chromosome
#' @export
#' @author Ruqian Lyu
#' 
calCoC <- function(count_co_gr,n_intervals=11,ref_genome="mm10",
                   group_prefix = NULL ){
  ## calculate C.o.C for all chromosomes
  ## n_intervals supplied
  ## fetch the chromoInfo from GenomeInfoDb. 
  ## This is only for getting the basepair lengths of the genome
  chrom_info <- GenomeInfoDb::fetchExtendedChromInfoFromUCSC(ref_genome)
  ## only for chr1-M
  chrom_info <- chrom_info[grep("_",chrom_info$UCSC_seqlevel,invert = TRUE),]
  
  ## Check what seqnames is in new_gr and make it consistent
  if(!grepl("chr",as.character(seqnames(count_co_gr)[1]))){
    chrom_info$UCSC_seqlevel <- gsub("chr","",chrom_info$UCSC_seqlevel)
  }
  chrom_info <- chrom_info[chrom_info$UCSC_seqlevel %in% 
                             GenomeInfoDb::seqlevels(count_co_gr),]
  ## create Granges object for chromosomes
  seq_length <- chrom_info$UCSC_seqlength
  names(seq_length) <- chrom_info$UCSC_seqlevel
  
  dna_mm10_gr <- GenomicRanges::GRanges(
    seqnames = Rle(names(seq_length)),
    ranges = IRanges(1, end = seq_length, names = names(seq_length)),
    seqlengths = seq_length)
  genome(dna_mm10_gr) <- ref_genome
  #dna_mm10_gr
  
  
  ## per bp distances
  GenomicRanges::mcols(count_co_gr) <- apply(GenomicRanges::mcols(count_co_gr),2,
                                             function(x) x/GenomicRanges::width(count_co_gr))
  
  ntile <- n_intervals
  tiles <- lapply(seqlevels(dna_mm10_gr),function(seqn){
    unlist(GenomicRanges::tileGenome(seqinfo(dna_mm10_gr)[seqn,],ntile = ntile))
  })
  
  
  binned_dna_mm10_gr <- suppressWarnings(do.call(c,tiles))
  # binned_dna_mm10_gr
  count_co_gr <- GenomicRanges::sort(GenomeInfoDb::sortSeqlevels(count_co_gr))
  
  bin_cos <- 
    BiocParallel::bplapply(mcols(count_co_gr),
                           function(sid){
                             dist_rle <- GenomicRanges::coverage(count_co_gr,
                                                                 weight = sid)
                             # runValue(dist_rle)[is.na(runValue(dist_rle))] <- 0
                             dist_bined <- GenomicRanges::binnedAverage(binned_dna_mm10_gr,dist_rle,
                                                                        "dist_bin_ave")
                             return(dist_bined$dist_bin_ave*width(dist_bined))})
  
  mcols(binned_dna_mm10_gr) <- do.call(cbind,bin_cos)
  ## if all samples in the same group /  no group prefix is provided 
  if(is.null(group_prefix)){
    coc <- cal_coc_by_group(binned_dna_mm10_gr)
  } else{
    coc <- bplapply(group_prefix,function(pref){
    cal_coc_by_group(binned_dna_mm10_gr,group_prefix = pref)})
    coc <- do.call(rbind,coc)
  }
  
  
  return(coc)

}