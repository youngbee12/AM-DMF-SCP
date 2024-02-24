#############################DIA-NN计算肽段数的方法########################


              #######安装DIA-NN R包###########
install.packages("devtools")
library(devtools)
install_github("https://github.com/vdemichev/diann-rpackage")
library(diann)

library(dplyr)
setwd("/Volumes/Yzc_data/2023中期报告/中期/23_DIA/")#设置工作目录

df <- diann_load("23_DIA.tsv")  #import data from diann_report.tsv




name_sample<-c(paste0("ACN_Trypsin_rep",1:3),paste0("ACN_Trypsin.LysC_rep",1:3),paste0("TFE_Trypsin_rep",1:3))

#precursors
data_precursors <- data.frame(diann_matrix(df, pg.q = 0.01))
colnames(data_precursors)<-name_sample
data_precursors_correct<-median_correct(data_precursors)

#Peptides
data_pep<-data.frame(diann_matrix(df, id.header="Stripped.Sequence", pg.q = 0.01))
colnames(data_pep)<-name_sample
data_pep_correct<-median_correct(data_pep)


#proteins
data_genes <- data.frame(diann_matrix(df, id.header="Genes", quantity.header="Genes.MaxLFQ.Unique", proteotypic.only = T, pg.q = 0.01))
colnames(data_genes)<-name_sample
data_genes_correct<-median_correct(data_genes)
data_genes_correct$Genename<-rownames(data_genes_correct)
#强度
newprotein_after<-gather(data_genes_correct,key="group",value="intensity",-Genename)#宽数据转长数据
rownames(newprotein_after)<-NULL
rownames(newprotein_after)<-newprotein_after$Genename
newprotein_after$group1<-split_string(newprotein_after$group)
split_at_second_underscore <- function(text) {
  #text<-"acn_tfe_s"
  underscores <- gregexpr("_", text)[[1]]
  
  if (length(underscores) >= 2) {
    split_point <- underscores[2]
    return(c(substr(text, 1, split_point-1), substr(text, split_point+1, nchar(text)))[1])
  } else {
    return(text[1])
  }
}

newprotein_after$group1 <- sapply(newprotein_after$group, split_at_second_underscore)




#
fasta_file_path<-"E:\\UP000005640_9606.fasta"

get_hydrovalue_kyte_doolittle_genelist<-function(fasta_file_path,geneList){
get_seq_from_fasta_genelist<-function(fasta_file_path=fasta_file_path,geneList=c("TP53")){
  library(Biostrings)
  library(dplyr)
  sequences <- readBStringSet(fasta_file_path)
  sq<-data.frame(sequences)
  sq$name<-rownames(sq)
  se1<-sq
  extract_gene_name <- function(str) {
    return(sub(".*GN=([^ ]+).*", "\\1", str))
  }
  
  gene_names <- sapply(sq$name, extract_gene_name)
  se1$gene_names<-gene_names
  se2<-se1[se1$gene_names %in% geneList,]
  
  
  return(se2)
}
a<-get_seq_from_fasta_genelist(fasta_file_path=fasta_file_path,geneList =geneList)
hydrophobicity_kyte_doolittle <- list(
  'A' = 1.8,  'C' = 2.5,  'D' = -3.5, 'E' = -3.5, 
  'F' = 2.8,  'G' = -0.4, 'H' = -3.2, 'I' = 4.5,
  'K' = -3.9, 'L' = 3.8,  'M' = 1.9,  'N' = -3.5,
  'P' = -1.6, 'Q' = -3.5, 'R' = -4.5, 'S' = -0.8,
  'T' = -0.7, 'V' = 4.2,  'W' = -0.9, 'Y' = -1.3
)
calculate_protein_hydrophobicity <- function(protein_sequence) {
  total_score <- 0
  for (aa in strsplit(protein_sequence, '')[[1]]) {
    if (aa %in% names(hydrophobicity_kyte_doolittle)) {
      total_score <- total_score + hydrophobicity_kyte_doolittle[[aa]]
    }
  }
  return(total_score / nchar(protein_sequence))
}
results <- sapply(a$sequences, calculate_protein_hydrophobicity)
value<-data.frame(results)

return(value)
}

get_hydrovalue_ECS_genelist<-function(fasta_file_path,geneList){
  get_seq_from_fasta_genelist<-function(fasta_file_path=fasta_file_path,geneList=c("TP53")){
    library(Biostrings)
    library(dplyr)
    sequences <- readBStringSet(fasta_file_path)
    sq<-data.frame(sequences)
    sq$name<-rownames(sq)
    se1<-sq
    extract_gene_name <- function(str) {
      return(sub(".*GN=([^ ]+).*", "\\1", str))
    }
    
    gene_names <- sapply(sq$name, extract_gene_name)
    se1$gene_names<-gene_names
    se2<-se1[se1$gene_names %in% geneList,]
    
    
    return(se2)
  }
  a<-get_seq_from_fasta_genelist(fasta_file_path=fasta_file_path,geneList =geneList)
  hydrophobicity_Eisenberg_consensus_scale <- list(
    'A' = 0.62,  'C' = 0.29,  'D' = -0.9, 'E' = -0.74, 
    'F' = 1.2,  'G' = 0.48, 'H' = -0.4, 'I' = 1.4,
    'K' = -1.5, 'L' = 1.1,  'M' = 0.64,  'N' = -0.78,
    'P' = -0.12, 'Q' = -0.85, 'R' = -2.5, 'S' = -0.18,
    'T' = -0.05, 'V' = 1.1,  'W' = 0.81, 'Y' = 0.26
  )
  calculate_protein_hydrophobicity <- function(protein_sequence) {
    total_score <- 0
    for (aa in strsplit(protein_sequence, '')[[1]]) {
      if (aa %in% names(hydrophobicity_Eisenberg_consensus_scale)) {
        total_score <- total_score + hydrophobicity_Eisenberg_consensus_scale[[aa]]
      }
    }
    return(total_score / nchar(protein_sequence))
  }
  results <- sapply(a$sequences, calculate_protein_hydrophobicity)
  value<-data.frame(results)
  
  return(value)
}


v1<-get_hydrovalue_kyte_doolittle_genelist(fasta_file_path=fasta_file_path,geneList=Depro_1_3$DEG_up)
rownames(v1)<-NULL
v11<-unique(v1$results)

v2<-get_hydrovalue_kyte_doolittle_genelist(fasta_file_path=fasta_file_path,geneList=Depro_1_3$DEG_down)
rownames(v2)<-NULL
plot(v1$results)
plot(v2$results)


df1 <- data.frame(Group="Group1", Value=v1$results)
df2 <- data.frame(Group="Group2", Value=v2$results)
df <- rbind(df1, df2)
# 使用ggplot2绘制箱线图
library(ggplot2)
ggplot(df, aes(x=Group, y=Value)) +
  geom_boxplot() +
  ggtitle("Boxplots of Two Groups") +
  ylab("Value")
v1<-get_hydrovalue_ECS_genelist(fasta_file_path=fasta_file_path,geneList=Depro_1_3$DEG_up)
rownames(v1)<-NULL

v2<-get_hydrovalue_ECS_genelist(fasta_file_path=fasta_file_path,geneList=Depro_1_3$DEG_down)
rownames(v2)<-NULL

df1 <- data.frame(Group="Group1", Value=v1$results)
df2 <- data.frame(Group="Group2", Value=v2$results)
df <- rbind(df1, df2)
# 使用ggplot2绘制箱线图
library(ggplot2)
ggplot(df, aes(x=Group, y=Value)) +
  geom_boxplot() +
  ggtitle("Boxplots of Two Groups") +
  ylab("Value")
