#' plotGTFreq
#'
#' Function to plot the genotypes for all samples faceted by genotype

#' @importFrom plotly plot_ly subplot
#' @importFrom reshape2 melt
#' @inheritParams countGT
#' @author Ruqian Lyu
#' @export
#' @return A ggplot object
#' @examples
#' or_geno <- snp_geno[,grep("X",colnames(snp_geno))]
#' rownames(or_geno) <- paste0(snp_geno$CHR,"_",snp_geno$POS)
#' or_geno[1,] <- rep("Fail",dim(or_geno)[2])
#' cr_geno <- correctGT(or_geno,ref = snp_geno$C57BL.6J,
#'                     alt = snp_geno$FVB.NJ..i.)
#' ft_gt <- filterGT(cr_geno)
#' plotGTFreq(ft_gt)

plotGTFreq <- function(geno){

  # geno_table <- sapply(colnames(geno), function(sample){
  #   list(Het = ifelse(is.na(table(geno[,sample],useNA = "no")["Het"]),
  #                     0,
  #                     table(geno[,sample],useNA = "no")["Het"]),
  #        Homo_alt = ifelse(is.na(table(geno[,sample],useNA = "no")["Homo_alt"])
  #                          ,0,
  #                          table(geno[,sample],useNA = "no")["Homo_alt"]))
  # })
  #
  # geno_table <- data.frame(samples = colnames(geno_table),
  #                          geno1=  as.numeric(unlist(geno_table["Het",]))/nrow(geno),
  #                          geno2 = as.numeric(unlist(geno_table["Homo_alt",])/nrow(geno)),
  #                          stringsAsFactors = FALSE)
  # pltdf <- melt(geno_table)
  pltdf <- data.frame(geno) %>% tidyr::pivot_longer(colnames(geno),
                                                    names_to = "sample",
                                                    values_to="geno")

  # if(interactive){
  #   ply1 <- plot_ly(pltdf, x=~.data$samples,y=~.data$value,type = "scatter",
  #                   color = ~.data$variable,mode = "markers",colors = color_set)
  #   return(ply1)
  # } else {

  stplt1 <- ggplot(data = pltdf)+
    geom_bar(mapping = aes(x = sample, fill = geno),position = "fill")+
    theme_classic(base_size = 11)+
    ylab("Genotype Frequecies for each sample")+labs(color ="Genotype")+
    theme(axis.text.x = element_text(angle =-90))

  return(stplt1)

}
