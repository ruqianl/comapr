## code to prepare `DATASET` dataset goes here

#usethis::use_data("DATASET")
## making example data set for examples and test cases

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
