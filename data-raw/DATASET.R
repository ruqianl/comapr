## code to prepare `DATASET` dataset goes here

#usethis::use_data("DATASET")
## making example data set for examples and test cases

## 
library(SummarizedExperiment)
library(BiocParallele)

mutate_inform1 <- read_excel(path = "/mnt/mcfiles/rlyu/Projects/Snakemake_projects/yeln_2019_spermtyping/analysis/All_data_May_to_August_2019.xlsx",
                             sheet = "mutant_25-6-19_informative")
mutate_inform1 <- data.frame(mutate_inform1)
snp_geno <- mutate_inform1[c(1:10,30:41,60:70),]
parents_geno <- data.frame(ref = snp_geno$C57BL.6J,
                           alt = snp_geno$FVB.NJ..i.)
usethis::use_data(parents_geno)

snp_geno_gr <- GenomicRanges::GRanges(
  seqnames = snp_geno$CHR,
  ranges = IRanges::IRanges(start = snp_geno$POS,
                            width = 1))

GenomicRanges::mcols(snp_geno_gr) <- snp_geno[,grep("X",colnames(snp_geno))]
                               
usethis::use_data(snp_geno_gr)                                     


###### construct RangedSummarizedExperiment

co_counts_list <- bplapply(paste0("chr",1:2), function(chr){
  cocount1 <- readRDS(file= paste0("/mnt/mcscratch/rlyu/Projects/Snakemake_projects/yeln_2019_spermtyping/output/scCNV/HMM/WC_522/",
                                   chr,"-countCO.rds"))
  cocount2 <- readRDS(file= paste0("/mnt/mcscratch/rlyu/Projects/Snakemake_projects/yeln_2019_spermtyping/output/scCNV/HMM/WC_526/",
                                   chr,"-countCO.rds"))
  stopifnot(identical(cocount1[,0],cocount2[,0]))
  colnames(mcols(cocount1)) <- paste0("WC-522_",colnames(mcols(cocount1))) 
  colnames(mcols(cocount2)) <- paste0("WC-526_",colnames(mcols(cocount2)))
  mcols(cocount1) <- cbind(mcols(cocount1),mcols(cocount2))
  cocount1 <- cocount1[rowSums(as.matrix(mcols(cocount1)))>0,]
  cocount1
})

wg_co_gr <- suppressWarnings(do.call(c,co_counts_list))


wg_co_gr_sprs <- Matrix::Matrix(as.matrix(mcols(wg_co_gr)),sparse=T)

wg_co_ranges <- wg_co_gr

mcols(wg_co_ranges) <- NULL

# nrows <- 200; 
# ncols <- 6
# counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
# rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
#                      IRanges(floor(runif(200, 1e5, 1e6)), width=100),
#                      strand=sample(c("+", "-"), 200, TRUE),
#                      feature_id=sprintf("ID%03d", 1:200))
# 
# colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
#                      row.names=LETTERS[1:6])
# rse <- SummarizedExperiment(assays=SimpleList(counts=counts),
#                             rowRanges=rowRanges, colData=colData)
#rowRanges(rse)


col_data <- DataFrame(SampleName=substr(colnames(wg_co_gr_sprs),start = 1,
                                                  stop=6),
                                totalReads=rbinom(ncol(wg_co_gr_sprs),size=1e7,
                                                  prob = .8),
                      totalSNPcovered=rbinom(ncol(wg_co_gr_sprs),size=1e5,
                                             prob = .8))

sampleRse <- SummarizedExperiment(assays=SimpleList(counts=wg_co_gr_sprs),
                                  rowRanges=wg_co_ranges,colData = col_data)

sampleRse
rowData(sampleRse)


### Generate vi.mtx for demo 
library(Matrix)
barcode_s1 <- paste0("BC",1:5)
barcode_s2 <- paste0("BC",letters[1:5])

write.table(barcode_s1,file = paste0("inst/extdata/s1_barcodes.txt"),
            quote = F,row.names = F,col.names = F)
write.table(barcode_s2,file = paste0("inst/extdata/s2_barcodes.txt"),
            quote = F,row.names = F,col.names = F)

snp_anno <- data.frame(POS = c(3000,
                               3200,
                               4000,
                               4500,
                               5500,
                               6000),
                       REF="G",
                       ALT="A")

for(chr in paste0("chr",seq(1:5))){
  write.table(snp_anno,file = paste0("inst/extdata/s1_",chr,"_snpAnnot.txt"),
              quote = F,row.names = F,col.names = T)
  
}

snp_anno <- data.frame(POS = c(3000,
                               3100,
                               4000,
                               4500,
                               5500,
                               6000),
                       REF="G",
                       ALT="A")

for(chr in paste0("chr",seq(1:5))){
  write.table(snp_anno,file = paste0("inst/extdata/s2","_",chr,"_snpAnnot.txt"),
              quote = F,row.names = F,col.names = T)
  
}



vi_state <- Matrix(data=0,nrow = 6,ncol=5)
vi_state[c(1,2),1] <- 1
vi_state[c(5,6),1] <- 2

vi_state[c(1,3),2] <- 2
vi_state[c(5,6),2] <- 1
vi_state[c(1,2,3,5,6),3] <- 1

vi_state[c(1,2,6),4] <-2


vi_state[c(1,2,3,4),5] <-1
vi_state[6,5] <-2
vi_state
for(chr in paste0("chr",seq(1:5))){
  writeMM(vi_state,file=paste0("inst/extdata/s1_",chr,"_vi.mtx"))
}

vi_state <- Matrix(data=0,nrow = 6,ncol=5)
vi_state[c(1,2,6),1] <- 1

vi_state[c(1,3),2] <- 2
vi_state[c(5,6),2] <- 1
vi_state[c(1,2,3,5,6),3] <- 1

vi_state[c(1,2,6),4] <-2


vi_state[c(1,2,3,4),5] <-1
vi_state[6,5] <-2
vi_state
for(chr in paste0("chr",seq(1:5))){
  writeMM(vi_state,file=paste0("inst/extdata/s2_",chr,"_vi.mtx"))
}

segInfo <- data.frame(ithSperm = c("ithSperm0",
                                   "ithSperm0",
                                   "ithSperm1",
                                   "ithSperm1",
                                   "ithSperm2",
                                   "ithSperm3",
                                   "ithSperm4",
                                   "ithSperm4"),
                      Seg_start = c(3000,
                                    5500,
                                    3000,
                                    5500,
                                    3000,
                                    3000,
                                    3000,
                                    6000),
                      Seg_end = c(3200,6000,4000,
                                  6000,6000,6000,
                                  4500,6000),
                      logllRatio = c(100,150,150,150,
                                     350,200,250,50),
                      cSNP = c(2,2,2,2,5,3,4,1),
                      State = c(1,2,2,1,1,2,1,2))
for(chr in paste0("chr",seq(1:5))){
  write.table(segInfo,col.names = F,row.names = F,quote = F,
              file = paste0("inst/extdata/s1_",chr,"_viSegInfo.txt"))
}

segInfo <- data.frame(ithSperm = c("ithSperm0",
                                   "ithSperm1",
                                   "ithSperm1",
                                   "ithSperm2",
                                   "ithSperm3",
                                   "ithSperm4",
                                   "ithSperm4"),
                      Seg_start = c(3000,
                                    3000,
                                    5500,
                                    3000,
                                    3000,
                                    3000,
                                    6000),
                      Seg_end = c(6000,4000,
                                  6000,6000,6000,
                                  4500,6000),
                      logllRatio = c(200,150,150,350,200,250,50),
                      cSNP = c(3,2,2,5,3,4,1),
                      State = c(1,2,1,1,2,1,2))
for(chr in paste0("chr",seq(1:5))){
  write.table(segInfo,col.names = F,row.names = F,quote = F,
              file = paste0("inst/extdata/s2_",chr,"_viSegInfo.txt"))
}



