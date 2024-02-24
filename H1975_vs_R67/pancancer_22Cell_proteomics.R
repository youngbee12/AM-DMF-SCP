library(readxl)
library(dplyr)
cell_line_define<-read_xlsx("analysis/external_data/Pan_cancer_pro_2022_Cell/Cell_Lines_Details.xlsx",sheet = "Cell line details")
pan_pro_express<-read_xlsx("analysis/external_data/Pan_cancer_pro_2022_Cell/mmc3.xlsx",sheet = "Full protein matrix")
drug_response_data<-read_xlsx("analysis/external_data/Pan_cancer_pro_2022_Cell/GDSC2_fitted_dose_response_24Jul22.xlsx",sheet = "Sheet 1")
muta_luad_info<-read.csv("analysis/external_data/Pan_cancer_pro_2022_Cell/LUAD_Genetic_features_variant_Sat Nov 25 06_39_19 2023.csv")
#data_tranform
LUAD_cellline<-cell_line_define %>% filter(TYPE=="LUAD")


pan_pro_express$Project_Identifier<-unlist(lapply(pan_pro_express$Project_Identifier,function(x)unlist(strsplit(x, ";"))[2]))
#colnames(pan_pro_express)[-1]<-unlist(lapply(colnames(pan_pro_express)[-1],function(x)unlist(strsplit(x, ";|_"))[2]))
#ttt<-unlist(lapply(colnames(pan_pro_express)[-1],function(x)unlist(strsplit(x, ";|_"))[2]))
LUAD_pro_express<-pan_pro_express %>% dplyr::filter( Project_Identifier %in% LUAD_cellline$`Sample Name`)

LUAD_drug_response_data<-drug_response_data %>% filter(TCGA_DESC=="LUAD")

#VIM express
VIM_LUAN_pro_Express<-LUAD_pro_express %>% dplyr::select(Project_Identifier,contains("VIM")) 
colnames(VIM_LUAN_pro_Express) <-c("CELL_LINE_NAME","VIM")

#EGFR response PUTATIVE_TARGET
relate_drug<-lapply(LUAD_drug_response_data$PUTATIVE_TARGET,function(x)unlist(strsplit(x, ", ")))
unlist(lapply(lapply(LUAD_drug_response_data$PUTATIVE_TARGET,function(x)unlist(strsplit(x, ", "))),function(x)any(x %in% "EGFR" )))


EGFR_LUAD_drug_response_data <- LUAD_drug_response_data %>% 
  filter( unlist(lapply(lapply(LUAD_drug_response_data$PUTATIVE_TARGET,function(x)unlist(strsplit(x, ", "))),function(x)any(x %in% "EGFR" )))
) %>% 
  dplyr::select(CELL_LINE_NAME,DRUG_NAME,Z_SCORE,LN_IC50)
split_EGFR_LUAD_drug_response_data <- split(EGFR_LUAD_drug_response_data, EGFR_LUAD_drug_response_data$DRUG_NAME)
test<-left_join(VIM_LUAN_pro_Express,EGFR_LUAD_drug_response_data) 
TT<-na.omit(test)
cor.test(TT$VIM,TT$Z_SCORE)
library(ggplot2)
summary_data <- TT %>%
  group_by(DRUG_NAME) %>%
  summarize(r = cor.test(Z_SCORE , VIM)$estimate,
            pvalue = cor.test(Z_SCORE , VIM)$p.value)
summary_data1 <- TT %>%
  group_by(DRUG_NAME) %>%
  summarize(r = cor.test(Z_SCORE , VIM,alternative="greater")$estimate,
            pvalue = cor.test(Z_SCORE , VIM,alternative="greater")$p.value)
pdf("analysis/H1975_vs_R67/20231201_pancancer_LUAD_EGFR_VIM.pdf",width = 12,height = )
TT %>%
  ggplot(aes(x = VIM, y = Z_SCORE)) +
  #geom_line(aes(x=sample,y = expression,group = gene), color = "grey", size = 0.5) +
  geom_point() +
  scale_color_manual(values=c("#4682B4" ,"#AF46B4","#B47846","black"))+
  scale_fill_manual(values=c("#4682B4" ,"#AF46B4","#B47846","black"))+
  ylab("Resistant_Score")+
  xlab("VIM_Expression_Value")+
  facet_wrap(~ DRUG_NAME, scales = "free_x", nrow = 2) +
  geom_smooth(method = "lm",  col = "red", linetype = "dashed",linewidth=0.5, fill="#69b3a2", se=TRUE) +
  stat_smooth(method = lm, se = FALSE, formula = y ~ x) +
  geom_text(data = summary_data1, aes(x = 3, y = -4, label = paste0("r = ", round(r, 2), ", p = ", formatC(pvalue,format = "e", digits = 3))),
            hjust = 0, vjust = 1, color = "black", size = 4) +  
  theme_bw()+
  theme(
    legend.position = "none",
    #axis.line = element_line(size=1.5),
    text = element_text(size = 16, face="plain"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0),# Values for face are one of "plain", "italic", "bold" and "bold.italic"
    axis.text.y = element_text(color="black", size=12, face="plain"),
    axis.title.x = element_text(color="black", size=16, face="bold"),
    axis.title.y = element_text(color="black", size=16, face="bold"),
    legend.text=element_text(color="black", size=16, face="bold"),
    legend.title=element_text(color="black", size=16, face="bold"),
    #title = element_text(color="black", size=16, face="bold"),
    #axis.text.x = element_blank(),
    #axis.title.x = element_blank()
    #axis.ticks.x = element_blank(),
    #strip.text = element_blank(),
    #strip.background = element_blank()
  )

dev.off()

tt1<-TT %>% filter(DRUG_NAME=="Osimertinib")
tt2<-tt1[tt1$CELL_LINE_NAME %in%muta_luad_info_EGFR$Cell.Line.Name ,]
tt2 %>%
  ggplot(aes(x = VIM, y = Z_SCORE)) +
  #geom_line(aes(x=sample,y = expression,group = gene), color = "grey", size = 0.5) +
  geom_point() +
  scale_color_manual(values=c("#4682B4" ,"#AF46B4","#B47846","black"))+
  scale_fill_manual(values=c("#4682B4" ,"#AF46B4","#B47846","black"))+
  ylab("Resistant_Score")+
  xlab("VIM_Expression_Value")+
  #facet_wrap(~ DRUG_NAME, scales = "free_x", nrow = 2) +
  geom_smooth(method = "lm",  col = "red", linetype = "dashed",linewidth=0.5, fill="#69b3a2", se=TRUE) +
  stat_smooth(method = lm, se = FALSE, formula = y ~ x) +
  geom_text(data = summary_data, aes(x = 3, y = -4, label = paste0("r = ", round(r, 2), ", p = ", formatC(pvalue,format = "e", digits = 3))),
            hjust = 0, vjust = 1, color = "black", size = 4) +  
  theme_bw()+
  theme(
    legend.position = "none",
    #axis.line = element_line(size=1.5),
    text = element_text(size = 16, face="plain"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0),# Values for face are one of "plain", "italic", "bold" and "bold.italic"
    axis.text.y = element_text(color="black", size=12, face="plain"),
    axis.title.x = element_text(color="black", size=16, face="bold"),
    axis.title.y = element_text(color="black", size=16, face="bold"),
    legend.text=element_text(color="black", size=16, face="bold"),
    legend.title=element_text(color="black", size=16, face="bold"),
    #title = element_text(color="black", size=16, face="bold"),
    #axis.text.x = element_blank(),
    #axis.title.x = element_blank()
    #axis.ticks.x = element_blank(),
    #strip.text = element_blank(),
    #strip.background = element_blank()
  )

LUAD_pro_express<-data.frame(LUAD_pro_express)
rownames(LUAD_pro_express)<-LUAD_pro_express$Project_Identifier
LUAD_pro_express<-LUAD_pro_express[,-1]
LUAD_pro_express1 <- data.table(t(LUAD_pro_express))
colnames(LUAD_pro_express1) <-rownames(LUAD_pro_express)
LUAD_pro_express1$gene<-colnames(LUAD_pro_express)

LUAD_pro_express1$gene<-unlist(lapply(LUAD_pro_express1$gene,function(x)unlist(strsplit(x, "\\.|_"))[2]))

LUAD_pro_express2 <- LUAD_pro_express1 %>%
  group_by(gene) %>%
  summarise_if(is.numeric, funs(mean), na.rm = TRUE)%>%
  data.table()
LUAD_pro_express3<-data.frame(LUAD_pro_express2)
rownames(LUAD_pro_express3)<-LUAD_pro_express3$gene
LUAD_pro_express3[LUAD_pro_express3=="NaN"]<-min(LUAD_pro_express3[,-1],na.rm=T)/10
######HALLMARK PATHWAY SCORE######
library(GSVA)
#计算markers的score
library(msigdbr)
hall <- msigdbr(species = "Homo sapiens", category = "H")
kegg <-msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:KEGG")
reactome <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:REACTOME")
gobp<-msigdbr(species = "Homo sapiens", category = "C5",subcategory = "BP")
EMT_marker=hall %>% filter(gs_name=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") %>%.[,"gene_symbol"] %>% unlist(.)
library(GSVA)
gsva_out<-gsva(as.matrix(LUAD_pro_express3[,-1]), list(EMT_scores=EMT_marker),method="ssgsea", min.sz = 1,max.sz=Inf,verbose=T,parallel.sz=1)

gsva_out<-data.frame(t(gsva_out))
gsva_out$CELL_LINE_NAME<-colnames(LUAD_pro_express1[,1:62])

test2<-left_join(gsva_out,EGFR_LUAD_drug_response_data) 
TT2<-na.omit(test2)
cor.test(TT2$EMT_scores,TT2$Z_SCORE)
library(ggplot2)
summary_data2 <- TT2 %>%
  group_by(DRUG_NAME) %>%
  summarize(r = cor.test(Z_SCORE , EMT_scores)$estimate,
            pvalue = cor.test(Z_SCORE , EMT_scores)$p.value)
pdf("analysis/H1975_vs_R67/20231030_pancancer_LUAD_EGFR_EMT.pdf",width = 12,height = )
TT2 %>%
  ggplot(aes(x = EMT_scores, y = Z_SCORE)) +
  #geom_line(aes(x=sample,y = expression,group = gene), color = "grey", size = 0.5) +
  geom_point() +
  scale_color_manual(values=c("#4682B4" ,"#AF46B4","#B47846","black"))+
  scale_fill_manual(values=c("#4682B4" ,"#AF46B4","#B47846","black"))+
  ylab("Resistant_Score")+
  xlab("EMT_Score")+
  facet_wrap(~ DRUG_NAME, scales = "free_x", nrow = 2) +
  geom_smooth(method = "lm",  col = "red", linetype = "dashed",linewidth=0.5, fill="#69b3a2", se=TRUE) +
  stat_smooth(method = lm, se = FALSE, formula = y ~ x) +
  geom_text(data = summary_data2, aes(x = -0.25, y = -4, label = paste0("r = ", round(r, 2), ", p = ", formatC(pvalue,format = "e", digits = 3))),
            hjust = 0, vjust = 1, color = "black", size = 4) +  
  theme_bw()+
  theme(
    legend.position = "none",
    #axis.line = element_line(size=1.5),
    text = element_text(size = 16, face="plain"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),# Values for face are one of "plain", "italic", "bold" and "bold.italic"
    axis.text.y = element_text(color="black", size=12, face="plain"),
    axis.title.x = element_text(color="black", size=16, face="bold"),
    axis.title.y = element_text(color="black", size=16, face="bold"),
    legend.text=element_text(color="black", size=16, face="bold"),
    legend.title=element_text(color="black", size=16, face="bold"),
    #title = element_text(color="black", size=16, face="bold"),
    #axis.text.x = element_blank(),
    #axis.title.x = element_blank()
    #axis.ticks.x = element_blank(),
    #strip.text = element_blank(),
    #strip.background = element_blank()
  )

dev.off()
