BiocManager::install("GenomicAlignments")
install.packages("UniprotR")
install.packages("protti")
library(UniprotR)
library(protti)
file_content <- readLines("5assenfb.txt")
tab_separated_content <- gsub(" +", "\t", file_content)
writeLines(tab_separated_content, "5assenfb_tab_separated.txt")
data <- read.table("5assenfb.txt", header = TRUE, sep = "\t")  # Assuming tab-separated file
accession_numbers <- c("Q1HK80", "Q9ZZ57", "Q08GU4", "A0A8U0SPQ4", "F4YPU8")
library(biomaRt)
get_GO_terms <- function(accession_number) {
  mart <- useMart("unimart")
  dataset <- useDataset("unimart", "uniprot")
library(biomaRt)
 listMarts()
  
  # Define function to get GO terms
  get_GO_terms <- function(accession_number) {
    mart <- useMart("ensembl")
    listDatasets(mart)
    dataset <- useDataset("unimart", mart=mart)
    
    go_info <- GetProteinGOInfo(accession_number, mart = mart, dataset = dataset)
    go_terms <- go_info$go_id
    
    return(go_terms)
  }
  
  accession_numbers <- c("Q1HK80", "Q9ZZ57", "Q08GU4", "A0A8U0SPQ4", "F4YPU8")
  
  all_go_terms <- list()
  
  for (accession in accession_numbers) {
    go_terms <- get_GO_terms(accession)
    all_go_terms[[accession]] <- go_terms
  }
  
  print(all_go_terms)
  
  install.packages("ggplot2")
  library(ggplot2)
  
  PlotGoInfo(all_go_terms)
  
  library(ggplot2)
  plot <- ggplot(data = your_data) + geom_bar(aes(x = your_x_variable, y = your_y_variable), stat = "identity")
  plot
  ggsave("plot.png", plot)
  getwd()  
  
  