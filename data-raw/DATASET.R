## code to prepare `DATASET` dataset goes here

#usethis::use_data("DATASET")
## making example data set for examples and test cases

##
library(SummarizedExperiment)
library(BiocParallele)
geno_file <- "/mnt/mcfiles/rlyu/Projects/Snakemake_projects/yeln_2019_spermtyping/analysis/All_data_May_to_August_2019.xlsx"
mutate_inform1 <- read_excel(path = geno_file,
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




### Generate vi.mtx for demo
library(Matrix)
barcode_s1 <- paste0("BC",1:5)
barcode_s2 <- paste0("BC",letters[1:5])

write.table(barcode_s1,file = paste0("inst/extdata/s1_barcodes.txt"),
            quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(barcode_s2,file = paste0("inst/extdata/s2_barcodes.txt"),
            quote = FALSE,row.names = FALSE,col.names = FALSE)

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
              quote = FALSE,row.names = FALSE,col.names = TRUE)

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
              quote = FALSE,row.names = FALSE,col.names = TRUE)

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
  write.table(segInfo,col.names = FALSE,row.names = FALSE,quote = FALSE,
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
  write.table(segInfo,col.names = FALSE,row.names = FALSE,quote = FALSE,
              file = paste0("inst/extdata/s2_",chr,"_viSegInfo.txt"))
}

### Generate cocount for examples

demo_path <-paste0(system.file("extdata",package = "comapr"),"/")
s1_rse_state <- readHapState("s1",chroms=c("chr1"),
                             path=demo_path,barcodeFile=NULL,minSNP = 0,
                             minlogllRatio = 50,
                             bpDist = 100,maxRawCO=10,
                             minCellSNP = 1)

s2_rse_state <- readHapState("s2",chroms=c("chr1"),
                             path=demo_path,
                             barcodeFile=paste0(demo_path,"s2_barcodes.txt"),
                             minSNP = 0,
                             minlogllRatio = 50,
                             bpDist = 100,maxRawCO=10,
                             minCellSNP = 1)
colData(s1_rse_state)$sampleGroup <- "s1"

colData(s2_rse_state)$sampleGroup <- "s2"
twoSamples <- combineHapState(s1_rse_state,s2_rse_state)
#colData(s1_rse_state)
usethis::use_data(twoSamples)

coCount <- countCOs(twoSamples)

usethis::use_data(coCount)


