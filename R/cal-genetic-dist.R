#' check supplied mapping function is available
#' @param mapping_fun, the mapping function to use, one of c("k","h")
#' @noRd
.check_mapping_fun <- function(mapping_fun){
  if(!mapping_fun %in% c("k","h")){
    message("only Kosambi and Haldane mapping function are supported at current
            stage")
    quit()
  }

}
#' cal_marker_dist
#'
#' Covert the recombination rate to cM by mapping functions
#' @importFrom GenomicRanges mcols
#' @noRd
cal_marker_dist <- function(co_count,mapping_fun="k",group_by = NULL){
  ## just covert the recombination rate to cM by mapping functions

  new_gr <- co_count[,0]

  if(is.null(group_by)){
    co_rate <- rowMeans(as.matrix(mcols(co_count)))
    stopifnot(sum(co_rate>=0.5)==0)
    cm_dist <- .rb_to_dist(co_rate, mapping_fun=mapping_fun)

    switch(mapping_fun,
      k =     mcols(new_gr)$kosambi_cm <- cm_dist,
      h =     mcols(new_gr)$haldane_cm <- cm_dist)

  } else {
    ## group_by has the character of prefix of groups
     groups_rb <- bplapply(group_by,function(group_prefix,mapping_fun){
       sids <- grep(group_prefix,colnames(mcols(co_count)))
       co_rate <- rowMeans(as.matrix(mcols(co_count)[,sids]))
       cm_dist <- .rb_to_dist(co_rate, mapping_fun=mapping_fun)

      # if(mapping_fun=="k"){
      #   x <- 25*log((1+2*co_rate)/(1-2*co_rate))
      # }
      # if(mapping_fun=="h"){
      #  x <-  -50 * log(1 - 2 * co_rate)
      # }
      # switch(mapping_fun,
      #         k =  x <- mcols(new_gr)$kosambi_cm <- cm_dist,
      #         h =   x <- mcols(new_gr)$haldane_cm <- cm_dist)

      cm_dist
    },mapping_fun=mapping_fun)
     mcols(new_gr) <- do.call(cbind,groups_rb)
     colnames(mcols(new_gr)) <- paste0(group_by,"_",mapping_fun)

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
#' GRange or RangedSummarizedExperiment object, returned by \code{countCO}
#' @param bin_size
#' The binning size for grouping marker intervals into bins. If not supplied,the
#' orginial marker intervals are returned with converted genetic distancens
#' based on recombination rate
#' @param mapping_fun
#' The mapping function to use, can be one of "k" or "h" (kosambi or haldane)
#' @param ref_genome
#' The reference genome name. It is used to fetch the chromosome size
#' information from UCSC database.
#' @param chrom_info
#' A user supplied data.frame containing two columns with column names
#' chrom and size, describing the chromosome names and lengths if not using
#' ref_genome from UCSC. If supplied, the `ref_genome` is ignored.
#'
#' @param group_by, character vector contains the unique prefix of sample names
#' that are used for defining different sample groups. Or the column name in
#' colData(co_count) that specify the group factor. If missing all samples are
#' assumed to be from one group
#'
#' @importFrom rlang .data
#' @importFrom GenomeInfoDb getChromInfoFromUCSC genome genome<-
#' @importFrom GenomicRanges GRanges mcolAsRleList sort width
#' @importFrom GenomicRanges tileGenome binnedAverage mcols coverage
#' @importFrom GenomeInfoDb sortSeqlevels seqlevels seqinfo
#' @importFrom S4Vectors Rle
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment colData rowRanges
#' @return GRanges object
#' GRanges for marker intervals or binned intervals with Haldane or Kosambi
#' centiMorgans
#' @examples
#' data(coCount)
#' dist_se <- calGeneticDist(coCount)
#' # dist_se <- calGeneticDist(coCount,group_by="sampleGroup")
#'
#' @export
#'

setGeneric("calGeneticDist",
           function(co_count,
                    bin_size=NULL,
                    mapping_fun="k",
                    ref_genome="mm10",
                    group_by = NULL,
                    chrom_info = NULL)
             standardGeneric("calGeneticDist"))

#'@noRd
calGenetic_dist <- function(co_count,
                            bin_size=NULL,
                            mapping_fun="k",
                            ref_genome="mm10",
                            group_by = NULL,
                            chrom_info = NULL){
  stopifnot(is.null(bin_size) | is.numeric(bin_size))
  .check_mapping_fun(mapping_fun)
  new_gr <- cal_marker_dist(co_count = co_count,
                          mapping_fun = mapping_fun,
                          group_by = group_by)
  if(is.null(bin_size)){
    return(new_gr)
  } else {
    binned_dna_mm10_gr <- cal_bin_dist(new_gr=new_gr,bin_size=bin_size,
                                       ref_genome="mm10",
                                       chrom_info = chrom_info)
    return(binned_dna_mm10_gr)
  }

}

#'@noRd
cal_bin_dist <- function(new_gr,bin_size,
                         ref_genome="mm10",chrom_info = NULL){
  ## bin_size supplied then.
  ## fetch the chromoInfo from GenomeInfoDb.
  ## This is only for getting the basepair lengths of the genome

  if(is.null(chrom_info)){
      out <- tryCatch({
      chrom_info <- getChromInfoFromUCSC(ref_genome)
      chrom_info <- chrom_info[grep("_",chrom_info$chrom,invert = TRUE),]
  },
  error=function(cond) {
    message("The required genome is not available from UCSC.\n
            You can supply a chrom_info in data.frame that have\n
            the chrom names and lengths columns with 'chrom','size'\n
            as column names")
    stop()
  })
}
  stopifnot(c("chrom","size") %in% colnames(chrom_info))
  ## Check what seqnames is in new_gr and make it consistent
  if(!grepl("chr",as.character(seqnames(new_gr)[1]))){
    chrom_info$chrom <- gsub("chr","",chrom_info$chrom)
  }
  chrom_info <- chrom_info[chrom_info$chrom %in%
                             seqlevels(new_gr),]
  ## create Granges object for chromosomes
  seq_length <- chrom_info$size
  names(seq_length) <- chrom_info$chrom

  dna_mm10_gr <- GRanges(
    seqnames = Rle(names(seq_length)),
    ranges = IRanges(1, end = seq_length, names = names(seq_length)),
    seqlengths = seq_length)
  genome(dna_mm10_gr) <- ref_genome
  #dna_mm10_gr


  ## per bp distances
  mcols(new_gr) <- apply(mcols(new_gr),2,
                                        function(x) {
                                          x/width(new_gr)})

  tilewidth <- bin_size
  tiles <- tileGenome(seqinfo(dna_mm10_gr),
                                     tilewidth = tilewidth)
  binned_dna_mm10_gr <- unlist(tiles)
  # binned_dna_mm10_gr
  new_gr <- sort(sortSeqlevels(new_gr))

  bin_dist <-  bplapply(colnames(mcols(new_gr)), function(group_col){
    # dist_rle <- GenomicRanges::mcolAsRleList(new_gr,group_col)
    # runValue(dist_rle)[is.na(runValue(dist_rle))] <- 0
    dist_rle <- coverage(new_gr,weight = mcols(new_gr)[,group_col])
    dist_bined <- binnedAverage(binned_dna_mm10_gr,dist_rle,
                                "dist_bin_ave")

    return(dist_bined$dist_bin_ave*width(dist_bined))

  })

  mcols(binned_dna_mm10_gr) <- do.call(cbind,bin_dist)
  colnames(mcols(binned_dna_mm10_gr)) <- colnames(mcols(new_gr))
  binned_dna_mm10_gr
}

#'@rdname calGeneticDist
setMethod("calGeneticDist",signature = c(co_count= 'GRanges',
                                         bin_size='missing',
                                         group_by='missing'),

          function(co_count,
                   bin_size=NULL,
                   mapping_fun="k",
                   ref_genome,
                   group_by = NULL,
                   chrom_info){
            .check_mapping_fun(mapping_fun)

            new_gr <- cal_marker_dist(co_count = co_count,
                                      mapping_fun = mapping_fun,
                                      group_by = NULL)
            new_gr
})
#'@rdname calGeneticDist
setMethod("calGeneticDist",signature = c(co_count = 'GRanges',
                                         bin_size='numeric',
                                         group_by='missing'),

          function(co_count,
                   bin_size=NULL,
                   mapping_fun="k",
                   ref_genome,
                   group_by = NULL,
                   chrom_info){
            .check_mapping_fun(mapping_fun)

            new_gr <- calGenetic_dist(co_count = co_count,
                                      mapping_fun = mapping_fun,
                                      bin_size = bin_size,
                                      ref_genome = ref_genome,
                                      group_by = NULL,
                                      chrom_info = chrom_info)
            new_gr
          })

#'@rdname calGeneticDist
setMethod("calGeneticDist",signature = c(co_count = 'GRanges',
                                         bin_size='missing',
                                         group_by='character'),

          function(co_count,
                   bin_size=NULL,
                   mapping_fun="k",
                   ref_genome,
                   group_by = NULL,
                   chrom_info){
            .check_mapping_fun(mapping_fun)

            new_gr <- cal_marker_dist(co_count = co_count,
                                      mapping_fun = mapping_fun,
                                      group_by = group_by)
            new_gr
})

#'@rdname calGeneticDist
setMethod("calGeneticDist",signature = c(co_count = 'GRanges',
                                         bin_size='numeric',
                                         group_by='character'),

          function(co_count,
                   bin_size,
                   mapping_fun="k",
                   ref_genome,
                   group_by,
                   chrom_info ){
            .check_mapping_fun(mapping_fun)

            new_gr <- calGenetic_dist(co_count = co_count,
                                      mapping_fun = mapping_fun,
                                      bin_size = bin_size,
                                      ref_genome = ref_genome,
                                      group_by = group_by,
                                      chrom_info = chrom_info)
            new_gr
          })

#' @keywords internal
#' @noRd
.setGenDistToRowRanges <- function(co_count,mapping_fun){
  cm_dist <- .rb_to_dist(rowRanges(co_count)$raw_rate, mapping_fun=mapping_fun)

  if(mapping_fun == "k")
    rowRanges(co_count)$kosambi <- cm_dist
  else
    rowRanges(co_count)$haldane <- cm_dist

  co_count
}

#'@rdname calGeneticDist
#'@importFrom SummarizedExperiment rowRanges rowRanges<-
#'
setMethod("calGeneticDist",
          signature = c(co_count = 'RangedSummarizedExperiment',
                        bin_size='missing',
                        group_by='missing'),

          function(co_count,
                   bin_size=NULL,
                   mapping_fun="k",
                   ref_genome="mm10",
                   group_by = NULL,
                   chrom_info){
            .check_mapping_fun(mapping_fun)
            co_rate <- rowMeans(as.matrix(assay(co_count)))
            rowRanges(co_count)$raw_rate <- co_rate
            # cm_dist <- .rb_to_dist(co_rate, mapping_fun=mapping_fun)
            # if(mapping_fun == "k")
            #   rowRanges(co_count)$kosambi <- cm_dist
            # else
            #   rowRanges(co_count)$haldane <-  cm_dist
            #
            # co_count
            co_count <- .setGenDistToRowRanges(co_count,
                                               mapping_fun = mapping_fun)
            co_count
          })

#'@rdname calGeneticDist
setMethod("calGeneticDist",
          signature = c(co_count = 'RangedSummarizedExperiment',
                        bin_size='missing',
                        group_by='character'),

          function(co_count,
                   bin_size=NULL,
                   mapping_fun="k",
                   ref_genome="mm10",
                   group_by,
                   chrom_info){
            .check_mapping_fun(mapping_fun)
            stopifnot(length(group_by)==1)
            ##group_by should refer to a column in colData(co_count)
            stopifnot(group_by %in% colnames(colData(co_count)))
            fcts <- colData(co_count)[,group_by]

            co_rate <- vapply(unique(as.character(fcts)),
                              function(fct){
                              rowMeans(as.matrix(assay(co_count)[,fcts==fct]))},
                              FUN.VALUE = numeric(nrow(assay(co_count))))

            rowRanges(co_count)$raw_rate <- co_rate
            # cm_dist <- .rb_to_dist(co_rate, mapping_fun=mapping_fun)
            #
            # if(mapping_fun == "k")
            #   rowRanges(co_count)$kosambi <- cm_dist
            # else
            #   rowRanges(co_count)$haldane <- cm_dist

            co_count <- .setGenDistToRowRanges(co_count,
                                               mapping_fun = mapping_fun)
            co_count
          })


#'@rdname calGeneticDist
setMethod("calGeneticDist",
          signature = c(co_count = 'RangedSummarizedExperiment',
                        bin_size='numeric',
                        group_by='character'),

          function(co_count,
                   bin_size,
                   mapping_fun="k",
                   ref_genome="mm10",
                   group_by,
                   chrom_info){

            .check_mapping_fun(mapping_fun)
            stopifnot(length(group_by)==1)
            ##group_by should refer to a column in colData(co_count)
            stopifnot(group_by %in% colnames(colData(co_count)))
            fcts <- colData(co_count)[,group_by]

            #####
            # divide the cos into bins
            ####
            co_rate <- vapply(unique(as.character(fcts)),
                              function(fct){
                              rowMeans(as.matrix(assay(co_count)[,fcts==fct]))},
                              FUN.VALUE = numeric(nrow(assay(co_count))))

            rowRanges(co_count)$raw_rate <- co_rate
            ## divide into bins
            new_gr <- granges(co_count)
            mcols(new_gr) <- .rb_to_dist(co_rate, mapping_fun = mapping_fun)
            dist_gr <- cal_bin_dist(new_gr,bin_size = bin_size,
                                    ref_genome = ref_genome,
                                    chrom_info = chrom_info)
            dist_gr

          })



#'@rdname calGeneticDist
setMethod("calGeneticDist",
          signature = c(co_count = 'RangedSummarizedExperiment',
                        bin_size='numeric',
                        group_by='missing'),

          function(co_count,
                   bin_size,
                   mapping_fun="k",
                   ref_genome="mm10",
                   group_by =NULL,
                   chrom_info){

            .check_mapping_fun(mapping_fun)

            #####
            # divide the cos into bins
            ####
            co_rate <- rowMeans(as.matrix(assay(co_count)))

            rowRanges(co_count)$raw_rate <- co_rate
            ## divide into bins
            new_gr <- granges(co_count)
            mcols(new_gr) <- .rb_to_dist(co_rate, mapping_fun = mapping_fun)
            dist_gr <- cal_bin_dist(new_gr,bin_size = bin_size,
                                    ref_genome = ref_genome,
                                    chrom_info = chrom_info)
            dist_gr

        })

