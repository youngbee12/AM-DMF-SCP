#############################DIA-NN计算肽段数的方法########################


              #######安装DIA-NN R包###########
#install.packages("devtools")
#library(devtools)
#install_github("https://github.com/vdemichev/diann-rpackage")
library(diann)


df <- diann_load("data/ALL/all_libbase/all_libbase_report-first-pass.tsv")  #import data from diann_report.tsv

polu_df<-df %>% filter(grepl("CON_",Protein.Group))

#Peptides
data_pep<-data.frame(diann_matrix(df %>% filter(!grepl("CON_",Protein.Group)), id.header="Stripped.Sequence", pg.q = 0.01))
data_pep_polu<-data.frame(diann_matrix(df %>% filter(grepl("CON_",Protein.Group)), id.header="Stripped.Sequence", pg.q = 0.01))

#split_string <- unlist(lapply(colnames(data_pep),function(x)unlist(strsplit(x, "yzc_"))[2]))
#split_string1 <- unlist(lapply(split_string,function(x)unlist(strsplit(x, "\\."))[1]))
#colnames(data1) <- split_string1


pep_group_info <- group_info %>% arrange(sample)
pep_group_info$rank <- c(rep(1,18),rep(2,18),rep(3,18),rep(4,18),rep(5,17),rep(7,18),rep(6,6),rep(8,6))
pep_group_info<-pep_group_info %>% arrange(rank)
colnames(data_pep)<-pep_group_info$sample
colnames(data_pep_polu)<-pep_group_info$sample

data_pep_correct<-median_correct(data_pep)
data_pep_polu_correct<-median_correct(data_pep_polu)


#proteins
data_genes_polu <- data.frame(diann_matrix(polu_df, id.header="Protein.Group", quantity.header="Genes.MaxLFQ.Unique", proteotypic.only = T, pg.q = 0.01))
colnames(data_genes_polu)<-pep_group_info$sample
data_genes_correct<-median_correct(data_genes)

#proteins
data_genes <- data.frame(diann_matrix(df, id.header="Protein.Group", quantity.header="Genes.MaxLFQ.Unique", proteotypic.only = T, pg.q = 0.01))
colnames(data_genes)<-pep_group_info$sample
data_genes_correct<-median_correct(data_genes)






#强度
newprotein_after<-gather(data_genes_correct,key="group",value="intensity",-Genename)#宽数据转长数据
rownames(protein_result_after)<-NULL
rownames(protein_result_after)<-protein_result_after$Genename

ggplot(newprotein_after,aes(color=group,x=group,y=log2(as.numeric(intensity))))+
  geom_boxplot()+
  labs(title="校正后")+
  ylab("Log2(intensity)")+
  theme_bw() +
  
  theme(  legend.position = "none",
          axis.line = element_line(size=1.5),
          axis.text.x = element_text(color="black", size=12, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
          axis.text.y = element_text(color="black", size=12, face="plain"),
          axis.title.x = element_text(color="black", size=16, face="bold"),
          axis.title.y = element_text(color="black", size=16, face="bold"),
          legend.text=element_text(color="black", size=16, face="bold"),
          legend.title=element_text(color="black", size=16, face="bold"),
          title = element_text(color="black", size=16, face="bold"))

cv_ACN_Trypsin<-apply(data_precursors_correct[,1:3],1,function(x)cv.coef(as.numeric(x)))
cv_ACN_Trypsin.LysC<-apply(data_precursors_correct[,4:6],1,function(x)cv.coef(as.numeric(x)))
cv_TFE_Trypsin<-apply(data_precursors_correct[,7:9],1,function(x)cv.coef(as.numeric(x)))

cv_ACN_pro<-apply(data_genes_correct[,1:5],1,function(x)cv.coef(as.numeric(x)))
cv_ACN_Trypsin.LysC_pro<-apply(data_genes_correct[,4:6],1,function(x)cv.coef(as.numeric(x)))
cv_TFE_Trypsin_pro<-apply(data_genes_correct[,7:9],1,function(x)cv.coef(as.numeric(x)))

median(cv_TFE_Trypsin_pro,na.rm = T)
median(cv_ACN_Trypsin.LysC_pro,na.rm = T)

data_cv<-data.frame(cv_ACN_Trypsin,cv_ACN_Trypsin.LysC,cv_TFE_Trypsin)
data_cv_pro<-data.frame(cv_ACN_pro,cv_ACN_Trypsin.LysC_pro,cv_TFE_Trypsin_pro)

data_cv$id<-rownames(data_cv)
data_cv_pro$id<-rownames(data_cv_pro)

data_cv_pro_median <- data_cv_pro %>%
  gather(key = "variable", value = "value", -id) %>%
  group_by(id, variable) %>%
  summarise(median_value = median(value, na.rm = TRUE))

data_cv_box<-gather(data_cv,key="group",value = "CV",-id)

data_cv_box_medians <- data_cv_box %>%
  group_by(group) %>%
  summarise(median_value = median(CV, na.rm = TRUE))

RColorBrewer::brewer.pal(5,"RdGy")
pdf(file="20230726_cv_pep.pdf")
ggplot(data_cv_box,aes(x=group,fill=group,y=as.numeric(100*CV)))+
  geom_violin()+
  geom_boxplot(width=0.2)+
  geom_text(data = data_cv_box_medians,
            aes(x=group,y = median_value*100, 
                label = round(median_value, 4)*100),
            vjust = 0, color = "white") +
   labs(title="violin_plot_CV_Pep")+
  ylab("CV(%)")+
  xlab("")+
  ylim(0,50)+
  scale_fill_manual(values=c("#24C9D4","#324B4D","#94B0B3","#B3AFF7","#7C7BBE"))+
  theme_bw() +
  
  theme(  legend.position = "none",
    axis.line = element_line(size=1.5),
    axis.text.x = element_text(color="black", size=12, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
    axis.text.y = element_text(color="black", size=12, face="plain"),
          axis.title.x = element_text(color="black", size=16, face="bold"),
          axis.title.y = element_text(color="black", size=16, face="bold"),
          legend.text=element_text(color="black", size=16, face="bold"),
          legend.title=element_text(color="black", size=16, face="bold"),
          title = element_text(color="black", size=16, face="bold"))
dev.off()
data_cv_box_pro<-gather(data_cv_pro,key="group",value = "CV",-id)

data_cv_box_pro_medians <- data_cv_box_pro %>%
  group_by(group) %>%
  summarise(median_value = median(CV, na.rm = TRUE))
pdf(file="20230726_cv_pro.pdf")
ggplot(data_cv_box_pro,aes(x=group,fill=group,y=as.numeric(100*CV)))+
  geom_violin()+
  geom_boxplot(width=0.2)+
  geom_text(data = data_cv_box_pro_medians,
            aes(x=group,y = median_value*100, 
                label = round(median_value, 4)*100),
            vjust = 0, color = "white") +
  labs(title="violin_plot_cv_pro")+
  ylab("CV(%)")+
  xlab("")+
  ylim(0,50)+
  scale_fill_manual(values=c("#24C9D4","#324B4D","#94B0B3","#B3AFF7","#7C7BBE"))+
  theme_bw() +
  
  theme(  legend.position = "none",
          axis.line = element_line(size=1.5),
          axis.text.x = element_text(color="black", size=12, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
          axis.text.y = element_text(color="black", size=12, face="plain"),
          axis.title.x = element_text(color="black", size=16, face="bold"),
          axis.title.y = element_text(color="black", size=16, face="bold"),
          legend.text=element_text(color="black", size=16, face="bold"),
          legend.title=element_text(color="black", size=16, face="bold"),
          title = element_text(color="black", size=16, face="bold"))
dev.off()
#蛋白鉴定量
name_sample<-c(paste0("ACN_Trypsin_rep",1:3),paste0("ACN_Trypsin.LysC_rep",1:3),paste0("TFE_Trypsin_rep",1:3))

nacol<-data.frame(nrow(data_genes_correct) - colSums ( is.na ( data_genes_correct ) ) )
colnames(nacol)<-"p_num"
p_num <- data.frame(matrix(nacol$p_num, ncol = 3, byrow = TRUE))
colnames(p_num)<-paste0("rep",1:3)
rownames(p_num)<-c("ACN_Trypsin","ACN_Trypsin.LysC","TFE_Trypsin")
mean_num<-apply(p_num,1,mean)
sd_num<-apply(p_num,1,sd)
p_num$Group<-rownames(p_num)
p_num$mean<-mean_num
p_num$sd<-sd_num

p_num_box<-gather(p_num[,c(1:4)],key="group",value = "p_num",-Group)
pdf(file="20230726_num_pro.pdf")
ggplot(p_num, aes(x=Group, y=mean,color=Group,fill=Group)) +
  geom_bar(position=position_dodge(), stat="identity",colour='white',width=.6) +
  # scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.4)+
  scale_color_manual(values=c("#24C9D4","#324B4D","#94B0B3","#B3AFF7","#7C7BBE"))+
  scale_fill_manual(values=c("#24C9D4","#324B4D","#94B0B3","#B3AFF7","#7C7BBE"))+
  geom_jitter(data = p_num_box,aes(x=Group,y=p_num),width=.1)+
  scale_y_continuous(limits=c(0,8500),breaks=seq(0,8000,1000))+
  geom_text(aes(label = ceiling(mean)), vjust = 5, colour = "white")+
  ylab("Proteins")+
  xlab("")+
  theme_bw() +
  
  theme(  legend.position = "none",
          axis.line = element_line(size=1.5),
          axis.text.x = element_text(color="black", size=12, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
          axis.text.y = element_text(color="black", size=12, face="plain"),
          axis.title.x = element_text(color="black", size=16, face="bold"),
          axis.title.y = element_text(color="black", size=16, face="bold"),
          legend.text=element_text(color="black", size=16, face="bold"),
          legend.title=element_text(color="black", size=16, face="bold"),
          title = element_text(color="black", size=16, face="bold"))

dev.off()
#肽段鉴定量

nacol_pep<-data.frame(nrow(data_pep_correct) - colSums ( is.na ( data_pep_correct ) ) )
colnames(nacol_pep)<-"pep_num"
pep_num <- data.frame(matrix(nacol_pep$pep_num, ncol = 3, byrow = TRUE))
colnames(pep_num)<-paste0("rep",1:3)
rownames(pep_num)<-c("ACN_Trypsin","ACN_Trypsin.LysC","TFE_Trypsin")
pepmean_num<-apply(pep_num,1,mean)
pepsd_num<-apply(pep_num,1,sd)
pep_num$Group<-rownames(pep_num)
pep_num$mean<-pepmean_num
pep_num$sd<-pepsd_num
pep_num$Group

pep_num_bar<-gather(pep_num[,c(1:4)],key="group",value = "pep_num",-Group)
pdf(file="20230726_num_pep.pdf")
ggplot(pep_num, aes(x=Group, y=mean,color=Group,fill=Group)) +
  geom_bar(position=position_dodge(), stat="identity",colour='white',width=.6) +
  # scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.4)+
  scale_color_manual(values=c("#24C9D4","#324B4D","#94B0B3","#B3AFF7","#7C7BBE"))+
  scale_fill_manual(values=c("#24C9D4","#324B4D","#94B0B3","#B3AFF7","#7C7BBE"))+
  geom_jitter(data = pep_num_bar,aes(x=Group,y=pep_num),width=.1)+
  scale_x_discrete(labels = c("value1" = "Metric 1", "value2" = "Metric 2", "value3" = "Metric 3")) +
  # scale_y_continuous(limits=c(0,8500),breaks=seq(0,8000,500))+
  geom_text(aes(label = ceiling(mean)), vjust = 5, colour = "white")+
  ylab("Peptides")+
  xlab("")+
  theme_bw() +
  
  theme(  legend.position = "none",
          axis.line = element_line(size=1.5),
          axis.text.x = element_text(color="black", size=12, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
          axis.text.y = element_text(color="black", size=12, face="plain"),
          axis.title.x = element_text(color="black", size=16, face="bold"),
          axis.title.y = element_text(color="black", size=16, face="bold"),
          legend.text=element_text(color="black", size=16, face="bold"),
          legend.title=element_text(color="black", size=16, face="bold"),
          title = element_text(color="black", size=16, face="bold"))
dev.off()

#漏切率

miss_cle<-c(0.22,0.22,0.22,0.20,0.20,0.20,0.30,0.30,0.30)
miss_cle_bar <- data.frame(matrix(miss_cle, ncol = 3, byrow = TRUE))
colnames(miss_cle_bar)<-paste0("rep",1:3)
rownames(miss_cle_bar)<-c("ACN_Trypsin","ACN_Trypsin.LysC","TFE_Trypsin")

mismean_num<-apply(miss_cle_bar,1,mean)
missd_num<-apply(miss_cle_bar,1,sd)
miss_cle_bar$Group<-rownames(miss_cle_bar)
miss_cle_bar$mean<-mismean_num
miss_cle_bar$sd<-missd_num

miss_cle_bar_pl<-gather(miss_cle_bar[,c(1:4)],key="group",value = "miss",-Group)
pdf(file="20230726_miss_pep.pdf")
ggplot(miss_cle_bar, aes(x=Group, y=mean*100,color=Group,fill=Group)) +
  geom_bar(position=position_dodge(), stat="identity",colour='white',width=.6) +
  # scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  geom_errorbar(aes(ymin=(mean-sd)*100, ymax=(mean+sd)*100, width=.4))+
  scale_color_manual(values=c("#24C9D4","#324B4D","#94B0B3","#B3AFF7","#7C7BBE"))+
  scale_fill_manual(values=c("#24C9D4","#324B4D","#94B0B3","#B3AFF7","#7C7BBE"))+
  geom_jitter(data = miss_cle_bar_pl,aes(x=Group,y=miss*100),width=.1)+
  # scale_y_continuous(limits=c(0,8500),breaks=seq(0,8000,500))+
  geom_text(aes(label = mean*100), vjust = 5, colour = "white",)+
  ylab("Average missed tryptic cleavages(%)")+
  xlab("")+
  theme_bw() +
  
  theme(  legend.position = "none",
          axis.line = element_line(size=1.5),
          axis.text.x = element_text(color="black", size=12, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
          axis.text.y = element_text(color="black", size=12, face="plain"),
          axis.title.x = element_text(color="black", size=16, face="bold"),
          axis.title.y = element_text(color="black", size=16, face="bold"),
          legend.text=element_text(color="black", size=16, face="bold"),
          legend.title=element_text(color="black", size=16, face="bold"),
          title = element_text(color="black", size=16, face="bold"))
dev.off()
#protein upset plot

data_genes$Genes<-rownames(data_genes)
pro_overlap_bar<-gather(data_genes,key="group",value = "quantity",-Genes)
ggplot(pep_num_bar,aes(x=Genres)) +
  geom_bar() +
  scale_x_upset(n_intersections = 20)
x<-c(1,2,NA,NA,3)
rownames(miss_cle_bar)<-c("ACN_Trypsin","ACN_Trypsin.LysC","TFE_Trypsin")
ACN_Trypsin_p_na<-apply(data_genes[,1:3],1,function(x)ifelse(sum(is.na(x))>2,NA,"ACN_Trypsin"))
ACN_Trypsin.LysC_p_na<-apply(data_genes[,4:6],1,function(x)ifelse(sum(is.na(x))>2,NA,"ACN_Trypsin.LysC"))
TFE_Trypsin_p_na<-apply(data_genes[,7:9],1,function(x)ifelse(sum(is.na(x))>2,NA,"TFE_Trypsin"))


data_na_gene<-data.frame(ACN_Trypsin_p_na,ACN_Trypsin.LysC_p_na,TFE_Trypsin_p_na)
colnames(data_na_gene)<-c("ACN_Trypsin","ACN_Trypsin.LysC","TFE_Trypsin")
#删除全是NA的列

pro_upset_list<-apply(data_na_gene,1,function(x)list(na.omit(x))) #准备成list格式
pro_upset<-as.tibble(pro_upset_list)#转化成tibble格式
pro_upset_final<-data.frame(t(pro_upset))#输入格式
colnames(pro_upset_final)<-"upsetR"
pdf(file="20230726_upset_pro.pdf")
ggplot(pro_upset_final,aes(x=upsetR)) +
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
  geom_bar() +
  ylim(0,8800)+
  scale_x_upset(n_intersections = 50)+
  ylab("")+
  xlab("")+
  
  theme_bw() +
  
  theme(  legend.position = "none",
          axis.line = element_line(size=1.5),
          axis.text.x = element_text(color="black", size=12, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
          axis.text.y = element_text(color="black", size=12, face="plain"),
          axis.title.x = element_text(color="black", size=16, face="bold"),
          axis.title.y = element_text(color="black", size=16, face="bold"),
          legend.text=element_text(color="black", size=16, face="bold"),
          legend.title=element_text(color="black", size=16, face="bold"),
          title = element_text(color="black", size=16, face="bold"))
  #theme(text = element_text(family = "Arial"))


dev.off()
#Venn to extract Gene
ACV_venn<-names(na.omit(ACN_p_na))
DDM_venn<-names(na.omit(DDM_p_na))
Rapigest_venn<-names(na.omit(Rapigest_p_na))
SDC_venn<-names(na.omit(SDC_p_na))
TFE_venn<-names(na.omit(TFE_p_na))
list_veen<-list(ACV_venn,DDM_venn,Rapigest_venn,SDC_venn,TFE_venn)
library(Vennerable)
temp=Venn(list_veen)  #provide all your groups as list


ACN_unique_pro<-temp@IntersectionSets[["10000"]]
DDM_unique_pro<-temp@IntersectionSets[["01000"]]
Rapigest_unique_pro<-temp@IntersectionSets[["00100"]]
SDC_unique_pro<-temp@IntersectionSets[["00010"]]
TFE_unique_pro<-temp@IntersectionSets[["00001"]]
all_pro<-temp@IntersectionSets[["11111"]]

DDM_unique_pro_clus<-clustern_kegg_go_list(DDM_unique_pro)


#peptides upset plot

data_pep_correct$peptides<-rownames(data_pep_correct)
pep_overlap_bar<-gather(data_pep_correct,key="group",value = "quantity",-peptides)


ACN_Trypsin_pep_na<-apply(data_pep_correct[,1:3],1,function(x)ifelse(sum(is.na(x))>2,NA,"ACN_Trypsin"))
ACN_Trypsin.LysC_pep_na<-apply(data_pep_correct[,4:6],1,function(x)ifelse(sum(is.na(x))>2,NA,"ACN_Trypsin.LysC"))
TFE_Trypsin_pep_na<-apply(data_pep_correct[,7:9],1,function(x)ifelse(sum(is.na(x))>2,NA,"TFE_Trypsin"))


data_na_pep<-data.frame(ACN_Trypsin_pep_na,ACN_Trypsin.LysC_pep_na,TFE_Trypsin_pep_na)
colnames(data_na_pep)<-c("ACN_Trypsin","ACN_Trypsin.LysC","TFE_Trypsin")

#删除全是NA的列

pep_upset_list<-apply(data_na_pep,1,function(x)list(na.omit(x))) #准备成list格式
pep_upset<-as.tibble(pep_upset_list)#转化成tibble格式
pep_upset_final<-data.frame(t(pep_upset))#输入格式
colnames(pep_upset_final)<-"upsetR"
pdf(file="20230726_upset_pep.pdf")
ggplot(pep_upset_final,aes(x=upsetR)) +
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
  geom_bar() +
   ylim(0,80000)+
  scale_x_upset(n_intersections = 50)+
  ylab("")+
  xlab("")+
  
  theme_bw() +
  
  theme(  legend.position = "none",
          axis.line = element_line(size=1.5),
          axis.text.x = element_text(color="black", size=12, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
          axis.text.y = element_text(color="black", size=12, face="plain"),
          axis.title.x = element_text(color="black", size=16, face="bold"),
          axis.title.y = element_text(color="black", size=16, face="bold"),
          legend.text=element_text(color="black", size=16, face="bold"),
          legend.title=element_text(color="black", size=16, face="bold"),
          title = element_text(color="black", size=16, face="bold"))
  #theme(text = element_text(family = "Arial"))

dev.off()

#热图colnames(data_na_pep)<-c("ACN_Trypsin","ACN_Trypsin.LysC","TFE_Trypsin")

heatmap_plot(data=log2(data_genes_correct))
group_info<-data.frame(sample=colnames(data_genes_correct),
                       group=c(rep("ACN_Trypsin",3),rep("ACN_Trypsin.LysC",3),rep("TFE_Trypsin",3)))
data_heat<-filter_na_bygroup(data = log2(data_genes_correct),group_info,x=0.2 )
pdf(file="20230726_heatmap_pro.pdf")
a<-heatmap_plot(data_heat)
a
dev.off()
heatmap_plot<-function(data,annotation_col=data.frame(Groups=factor(c(rep("ACN_Trypsin",3),rep("ACN_Trypsin.LysC",3),rep("TFE_Trypsin",3))))){
  library(pheatmap)
  library(dplyr)
  annotation_col = annotation_col

  ann_colors = list(Groups = c(ACN_Trypsin = '#24C9D4',ACN_Trypsin.LysC = '#324B4D',TFE_Trypsin="#94B0B3"))
  rownames(annotation_col)<-colnames(data)
  bks<-seq(-1,1,length.out=10000)
  
                        
  heatmap<-pheatmap(data,
                    breaks=bks,
                    na_col = "grey",
                    scale = "row",
                    show_rownames=F,
                    show_colnames=TRUE,
                    #clustering_distance_cols = "correlation",
                    #clustering_method = "median",
                    cluster_rows = T,
                    cluster_cols = F,
                    cellwidth = 16, 
                    cellheight = 0.02,
                    display_numbers = FALSE,
                    annotation_col = annotation_col,
                    annotation_colors = ann_colors,
                    fontsize =12,
                    color = colorRampPalette(c("royalblue4","white","firebrick3"))(10000),
                    
                    # gaps_col = 3,
                    # gaps_row = 3,
  )
  return(heatmap)
}

data_heat<-filter_na_bygroup(data = log2(data_genes_correct),group_info,x=0.5 )
Depro_1_2<-wil_cox_nopair_test(data_heat,n1=1,n2=3,n3=4,n4=6,p_value=0.05,Fc_value=2,method="BH")
Depro_1_3<-wil_cox_nopair_test(data_heat,n1=1,n2=3,n3=7,n4=9,p_value=0.05,Fc_value=2,method="BH")
Depro_2_3<-wil_cox_nopair_test(data_heat,n1=4,n2=6,n3=7,n4=9,p_value=0.05,Fc_value=2,method="BH")


write.csv(Depro_2_3$data,file="20230726_depro_2_3.csv")
write.csv(Depro_1_3$data,file="20230726_depro_1_3.csv")
write.csv(Depro_1_2$data,file="20230726_depro_1_2.csv")
pdf(file="20230726_volpro_2_3.pdf")
volpro_2_3<-vocanol_plot(Depro_2_3$data,p.value=0.05,Fc_value=2)
volpro_2_3
dev.off()
pdf(file="20230726_volpro_1_3.pdf")
volpro_1_3<-vocanol_plot(Depro_1_3$data,p.value=0.05,Fc_value=2)
volpro_1_3
dev.off()

pdf(file="20230726_volpro_1_2.pdf")
volpro_1_2<-vocanol_plot(Depro_1_2$data,p.value=0.05,Fc_value=2)+
  scale_color_manual(values=c("royalblue4","grey70"))
volpro_1_2
dev.off()


overlap_two_vector(Depro_2_3$DEG_down,Depro_1_3$DEG_down)
  wil_cox_nopair_test<-function(data_log2,n1=1,n2=120,n3=121,n4=240,p_value=0.01,Fc_value=2,method="BH"){
    data<-data_log2
    data_BH <-  p.adjust(apply(data, 1, function(x)
      t.test(x[seq(n1,n2,1)], x[seq(n3,n4,1)],
             paired = F, exact = FALSE)$p.value), method = method)
    data$BH_adjust<-data_BH
    data_FC <- apply(data, 1, function(x)(
      mean(x[seq(n3,n4,1)],na.rm = T)-mean(x[seq(n1,n2,1)],na.rm = T))
    )
    data$Log2FC<-data_FC
    DEG_up <- intersect(names(data_BH[data_BH < p_value]), names(data_FC[data_FC > log2(Fc_value)]))
    DEG_down<- intersect(names(data_BH[data_BH < p_value]), names(data_FC[data_FC < -log2(Fc_value)]))
    listx<-list()
    listx$DEG_up<-DEG_up
    listx$DEG_down<-DEG_down
    listx$data<-data
    cat("上调基因有",length(DEG_up),"个，下调基因有",length(DEG_down),"个。")
    return(listx)
    
  }


  
#回收

  
  rec<-c(41.25,45.99,39.93,31.83,33.06,34.26,28.89,24.54,21.45)
  rec_bar <- data.frame(matrix(rec, ncol = 3, byrow = TRUE))
  colnames(rec_bar)<-paste0("rep",1:3)
  rownames(rec_bar)<-c("ACN_Trypsin","ACN_Trypsin.LysC","TFE_Trypsin")
  
  recmean_num<-apply(rec_bar,1,mean)
  rec_num<-apply(rec_bar,1,sd)
  rec_bar$Group<-rownames(rec_bar)
  rec_bar$mean<-recmean_num
  rec_bar$sd<-rec_num
  
  rec_bar_pl<-gather(rec_bar[,c(1:4)],key="group",value = "miss",-Group)
  pdf(file="20230726_rec_bar.pdf")
  ggplot(rec_bar, aes(x=Group, y=mean,color=Group,fill=Group)) +
    geom_bar(position=position_dodge(), stat="identity",colour='white',width=.6) +
    # scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
    geom_errorbar(aes(ymin=(mean-sd), ymax=(mean+sd), width=.4))+
    scale_color_manual(values=c("#24C9D4","#324B4D","#94B0B3","#B3AFF7","#7C7BBE"))+
    scale_fill_manual(values=c("#24C9D4","#324B4D","#94B0B3","#B3AFF7","#7C7BBE"))+
    geom_jitter(data = miss_cle_bar_pl,aes(x=Group,y=miss*100),width=.1)+
    # scale_y_continuous(limits=c(0,8500),breaks=seq(0,8000,500))+
    #geom_text(aes(label = mean), vjust = 5, colour = "white",)+
    ylab("Recovery amounts (μg)")+
    xlab("")+
    theme_bw() +
    
    theme(  legend.position = "none",
            axis.line = element_line(size=1.5),
            axis.text.x = element_text(color="black", size=12, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
            axis.text.y = element_text(color="black", size=12, face="plain"),
            axis.title.x = element_text(color="black", size=16, face="bold"),
            axis.title.y = element_text(color="black", size=16, face="bold"),
            legend.text=element_text(color="black", size=16, face="bold"),
            legend.title=element_text(color="black", size=16, face="bold"),
            title = element_text(color="black", size=16, face="bold"))
  dev.off()
  















#Venn to extract Gene
ACN_venn_pep<-names(na.omit(ACN_pep_na))
DDM_venn_pep<-names(na.omit(DDM_pep_na))
Rapigest_venn_pep<-names(na.omit(Rapigest_pep_na))
SDC_venn_pep<-names(na.omit(SDC_pep_na))
TFE_venn_pep<-names(na.omit(TFE_pep_na))
list_veen_pep<-list(ACN_venn_pep,DDM_venn_pep,Rapigest_venn_pep,SDC_venn_pep,TFE_venn_pep)
library(Vennerable)
temp_pep=Venn(list_veen_pep)  #pepvide all your groups as list
ACN_unique_pep<-temp_pep@IntersectionSets[["10000"]]
DDM_unique_pep<-temp_pep@IntersectionSets[["01000"]]
Rapigest_unique_pep<-temp_pep@IntersectionSets[["00100"]]
SDC_unique_pep<-temp_pep@IntersectionSets[["00010"]]
TFE_unique_pep<-temp_pep@IntersectionSets[["00001"]]
all_pep<-temp_pep@IntersectionSets[["11111"]]
write.table(ACN_unique_pep,file="20230421_ACN_unique_pep.txt",sep = "")
write.table(DDM_unique_pep,file="20230421_DDM_unique_pep.txt",sep = "")
write.table(Rapigest_unique_pep,file="20230421_Rapigest_unique_pep.txt",sep = "")
write.table(SDC_unique_pep,file="20230421_SDC_unique_pep.txt",sep = "")
write.table(TFE_unique_pep,file="20230421_TFE_unique_pep.txt",sep = "")
getwd()
save.image(file = "/Volumes/YZC_X6/experiment/工作汇报/2023/20230425/20230419.Rdata")







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

df1 <- data.frame(Group="DEG_up", Value=v1$results)
df2 <- data.frame(Group="DEG_down", Value=v2$results)
df <- rbind(df1, df2)
compare_means( Value~Group,df, method = "t.test")
# 使用ggplot2绘制箱线图
library(ggplot2)
pdf(file="20230726_hydro_DE1_3.pdf")
ggplot(df, aes(x=Group, y=Value,color=Group)) +
  geom_boxplot() +
  #labs(title=paste0('boxplot_',gene))+
  geom_jitter( size=1, alpha=0.9)+
  scale_color_manual(values=c("#324B4D","#7C7BBE"))+
  ylim(-0.8,0.8)+
  ylab("ECS score")+
  xlab("")+
  ggtitle("Depro_1_3")+
  stat_compare_means(method = "t.test")+      # Add global p-value
  #stat_compare_means(label = "p.signif", method = "t.test")+                  # Pairwise comparison against all
  theme_bw() +
  
  theme(  legend.position = "none",
          axis.line = element_line(size=1.5),
          axis.text.x = element_text(color="black", size=12, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
          axis.text.y = element_text(color="black", size=12, face="plain"),
          axis.title.x = element_text(color="black", size=16, face="bold"),
          axis.title.y = element_text(color="black", size=16, face="bold"),
          legend.text=element_text(color="black", size=16, face="bold"),
          legend.title=element_text(color="black", size=16, face="bold"),
          title = element_text(color="black", size=16, face="bold"))
dev.off()
v1<-get_hydrovalue_ECS_genelist(fasta_file_path=fasta_file_path,geneList=Depro_2_3$DEG_up)
rownames(v1)<-NULL

v2<-get_hydrovalue_ECS_genelist(fasta_file_path=fasta_file_path,geneList=Depro_2_3$DEG_down)
rownames(v2)<-NULL

df1 <- data.frame(Group="DEG_up", Value=v1$results)
df2 <- data.frame(Group="DEG_down", Value=v2$results)
df <- rbind(df1, df2)


pdf(file="20230726_hydro_DE2_3.pdf")
ggplot(df, aes(x=Group, y=Value,color=Group)) +
  geom_boxplot() +
  #labs(title=paste0('boxplot_',gene))+
  geom_jitter( size=1, alpha=0.9)+
  scale_color_manual(values=c("#324B4D","#7C7BBE"))+
  ylim(-0.8,0.8)+
  ylab("ECS score")+
  xlab("")+
  ggtitle("Depro_2_3")+
  stat_compare_means(method = "t.test")+      # Add global p-value
  #stat_compare_means(label = "p.signif", method = "t.test")+                  # Pairwise comparison against all
  theme_bw() +
  
  theme(  legend.position = "none",
          axis.line = element_line(size=1.5),
          axis.text.x = element_text(color="black", size=12, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
          axis.text.y = element_text(color="black", size=12, face="plain"),
          axis.title.x = element_text(color="black", size=16, face="bold"),
          axis.title.y = element_text(color="black", size=16, face="bold"),
          legend.text=element_text(color="black", size=16, face="bold"),
          legend.title=element_text(color="black", size=16, face="bold"),
          title = element_text(color="black", size=16, face="bold"))
dev.off()



