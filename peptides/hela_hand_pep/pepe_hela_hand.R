########################
######导入数据##########
########################
library(dplyr)
library(PerformanceAnalytics)
library(corrplot)
library(RColorBrewer)
library(colortools)
library(ComplexHeatmap)
library(circlize)


load(file="./data/20231011_scdata_pre.Rdata")
sc_data_raw<-sc_data
load("./data/20231011_data_clean_sc_filterdatabygroup_80.Rdata")
data_pep




pep_hela_hand_sample<-pep_group_info %>% arrange(Group)%>%filter(Group=="outchip"|Group=="Hela"&batch=="batch1") 
hela_hand_sample[hela_hand_sample=="Hela"]<-"inchip"

data_pep_polu_select<-data_pep_polu[,hela_hand_sample$sample]
data_pep_select<-data_pep[,hela_hand_sample$sample]
com<-overlap_two_vector(listn[[1]],listn[[2]])$common
out_spe<-overlap_two_vector(listn[[1]],listn[[2]])$uni_y
in_spe<-overlap_two_vector(listn[[1]],listn[[2]])$uni_x

data_pep_select_inchip <- data_pep_select[com,1:6]
data_pep_select_outchip <- data_pep_select[com,7:12]
data_pep_select_inchip$mean_int<-apply(data_pep_select_inchip,1,function(x)mean(x,na.rm=T))
data_pep_select_outchip$mean_int<-apply(data_pep_select_outchip,1,function(x)mean(x,na.rm=T))

outchip_up_misc<-unlist(lapply(listn[[2]],function(x)pep_miscleva_num_calc(input = x)))
inchip_up_misc<-unlist(lapply(listn[[1]],function(x)pep_miscleva_num_calc(input = x)))
com_misc<-unlist(lapply(com,function(x)pep_miscleva_num_calc(input = x)))
out_spe_mis<-unlist(lapply(out_spe,function(x)pep_miscleva_num_calc(input = x)))
in_spe_mis<-unlist(lapply(in_spe,function(x)pep_miscleva_num_calc(input = x)))

cor(outchip_up_misc,inchip_up_misc)

length(inchip_up_misc[inchip_up_misc==0])/length(inchip_up_misc)
length(outchip_up_misc[outchip_up_misc==0])/length(outchip_up_misc)
length(com_misc[com_misc==0])/length(com_misc)
length(out_spe_mis[out_spe_mis==0])/length(out_spe_mis)
length(in_spe_mis[in_spe_mis==0])/length(in_spe_mis)

t.test(log2(data_pep_select_inchip$mean_int),log2(data_pep_select_outchip$mean_int),paired = F)


listn<-list()

nacol_pep<-data.frame(nrow(data_pep_select) - colSums ( is.na ( data_pep_select ) ) )
colnames(nacol_pep)<-"pep_num"
hela_hand_sample$pep_num <- nacol_pep$pep_num
wide_df_pep <- hela_hand_sample %>% 
  pivot_wider(names_from = Group, values_from = pep_num) %>% dplyr::select(-Rep,-sample,-batch,-batch_sample,-p_num) %>% t(.) %>% data.frame(.)
colnames(wide_df_pep)<-rownames(hela_hand_sample)
mean_num<-apply(wide_df_pep,1,function(x)mean(x,na.rm=T))
sd_num<-apply(wide_df_pep,1,function(x)sd(x,na.rm=T))
wide_df_pep$Group<-rownames(wide_df_pep)
wide_df_pep$mean<-mean_num
wide_df_pep$sd<-sd_num
pdf(file="analysis/peptides/hela_hand_pep/20231020_identi_pep.pdf",height = 5,width = 5 )
ggplot(wide_df_pep, aes(x=Group, y=mean,color=Group,fill=Group)) +
  geom_bar(position=position_dodge(), stat="identity",colour='white',width=.6,alpha=0.8) +
  # scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.4)+
  scale_color_manual(values=c("#4682B4","#AF46B4"))+
  scale_fill_manual(values=c("#4682B4","#AF46B4"))+
  geom_jitter(data = hela_hand_sample,aes(x=Group,y=pep_num),width=.1)+
 scale_y_continuous(limits=c(0,12000),breaks=seq(0,12000,3000))+
  stat_compare_means(data = hela_hand_sample,aes(x=Group, y=pep_num),method = "t.test",na.rm = T)+      # Add global p-value
  geom_text(aes(label = ceiling(mean)), vjust = -9, colour = "black")+
  ylab("Peptides")+
  xlab("")+
  theme_bw() +
  
  theme(  legend.position = "none",
          axis.line = element_line(size=1.5),
          #axis.text.x = element_text(color="black", size=12, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
          axis.text.x = element_text(color="black", size=16, face="plain") ,
          axis.text.y = element_text(color="black", size=12, face="plain"),
          axis.title.x = element_text(color="black", size=16, face="bold"),
          axis.title.y = element_text(color="black", size=16, face="bold"),
          legend.text=element_text(color="black", size=16, face="bold"),
          legend.title=element_text(color="black", size=16, face="bold"),
          title = element_text(color="black", size=16, face="bold"))

dev.off()







#data=data_total
n<-length(unique(hela_hand_sample$Group))
for (i in 1:n){
  t_sm<-hela_hand_sample[hela_hand_sample$Group==unique(hela_hand_sample$Group)[i],]$sample
  listn[i]<-list(data_pep_polu_select[,colnames(data_pep_polu_select) %in% t_sm] %>%
                   .[ rowSums ( is.na ( . ) )  < ncol(.)*0.5,] %>% rownames(.))
}
overlap_two_vector(listn[[1]],listn[[2]])
com<-overlap_two_vector(listn[[1]],listn[[2]])$common

pair_line_boxplot_data<-data_pep_polu_select[com,]
pair_line_boxplot_data$mean_inchip<-apply(pair_line_boxplot_data[,1:6],1,function(x)mean(x,na.rm=T))
pair_line_boxplot_data$mean_outchip<-apply(pair_line_boxplot_data[,7:12],1,function(x)mean(x,na.rm=T))

pair_line_boxplot_data_1<-log2(pair_line_boxplot_data[,13:14])
pair_line_boxplot_data_1$Genename<-rownames(pair_line_boxplot_data_1)
newprotein_after<-gather(pair_line_boxplot_data_1,key="group",value="intensity",-Genename)
pdf(file="analysis/peptides/hela_hand_pep/20231013_pep_Contamients.pdf",width = 5,height = 5)
ggplot(newprotein_after,aes(color=group,x=group,y=as.numeric(intensity)))+
  scale_color_manual(values=c("#4682B4" ,"#AF46B4"))+
  scale_fill_manual(values=c("#4682B4" ,"#AF46B4"))+
  geom_point(size=1)+
  geom_line(aes(group=Genename),color="grey",alpha=0.3)+
  geom_boxplot(width=0.3,alpha=0.2)+
# labs(title="Contamients_Pep")+
  ylab("Log2(Intensity)")+
  xlab("")+
  stat_compare_means(method = "t.test",paired = F)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "mean_inchip",paired = F)+      
  theme_bw() +
  
  theme(  legend.position = "none",
          #axis.line = element_line(size=1.5),
          axis.text.x = element_text(color="black", size=16, face="plain",angle=0,hjust=0.5) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
          axis.text.y = element_text(color="black", size=16, face="plain"),
          axis.title.x = element_text(color="black", size=16, face="bold"),
          axis.title.y = element_text(color="black", size=16, face="bold"),
          legend.text=element_text(color="black", size=16, face="bold"),
          legend.title=element_text(color="black", size=16, face="bold"),
          title = element_text(color="black", size=16, face="bold"))
dev.off()
t.test(pair_line_boxplot_data_1$mean_inchip,pair_line_boxplot_data_1$mean_outchip)
ggplot(pair_line_boxplot_data_1,)




pdf(file="./analysis/hela_hand/20231011_umap_hela_hand.pdf",height = 5,width = 6 )
overview_umap(data = data_hela_hand,data2 =hela_hand_sample,group = c("Group"),eclipse_level = 0,group_values = wheel("steelblue", num = 2),n_neighbors=10)
dev.off()
res.pca3 = PCA(t(data_hela_hand), scale.unit = TRUE, ncp = 8,  graph = T,axes = 1:3)
pca3<-fviz_pca_ind(res.pca3,
                   label = "none",
                   habillage = factor(hela_hand_sample$Group),
                   #palette =c("#989898", "#E69F00","#56B4E9") ,
                   addEllipses = T,
                   ellipse.level=0.95,
                   pointsize=4,
                   #title="n=9807,PCA T1 vs N "
)+
  theme_bw()+
  theme(
    
    axis.text.x = element_text(color="black", size=16, face="plain",angle=0,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
    axis.text.y = element_text(color="black", size=16, face="plain"),
    axis.title.x = element_text(color="black", size=16, face="bold"),
    axis.title.y = element_text(color="black", size=16, face="bold"),
    legend.text=element_text(color="black", size=16, face="bold"),
    legend.title=element_text(color="black", size=16, face="bold"),
    title = element_text(color="black", size=16, face="bold")
  ) 

data_genes_polu
data_pep_polu_select<-data_genes_polu[,hela_hand_sample$sample]
data_pep_select<-data_pep[,hela_hand_sample$sample]

listn<-list()
#data=data_total
n<-length(unique(hela_hand_sample$Group))
for (i in 1:n){
  t_sm<-hela_hand_sample[hela_hand_sample$Group==unique(hela_hand_sample$Group)[i],]$sample
  listn[i]<-list(data_pep_polu_select[,colnames(data_pep_polu_select) %in% t_sm] %>%
                   .[ rowSums ( is.na ( . ) )  < ncol(.)*0.5,] %>% rownames(.))
}
overlap_two_vector(listn[[1]],listn[[2]])
com<-overlap_two_vector(listn[[1]],listn[[2]])$common

pair_line_boxplot_data<-data_pep_polu_select[com,]
pair_line_boxplot_data$mean_inchip<-apply(pair_line_boxplot_data[,1:6],1,function(x)mean(x,na.rm=T))
pair_line_boxplot_data$mean_outchip<-apply(pair_line_boxplot_data[,7:12],1,function(x)mean(x,na.rm=T))

pair_line_boxplot_data_1<-log10(pair_line_boxplot_data[,13:14])
pair_line_boxplot_data_1$Genename<-rownames(pair_line_boxplot_data_1)
newprotein_after<-gather(pair_line_boxplot_data_1,key="group",value="intensity",-Genename)
  pdf(file="analysis/peptides/hela_hand_pep/20231013_pro_con.pdf",width = 5,height = 5)
ggplot(newprotein_after,aes(color=group,x=group,y=as.numeric(intensity)))+
  scale_color_manual(values=c("#4682B4" ,"#AF46B4"))+
  scale_fill_manual(values=c("#4682B4" ,"#AF46B4"))+
  geom_point(size=1)+
  geom_line(aes(group=Genename),color="grey",alpha=0.3)+
  #geom_jitter(width=0.3,alpha=0.2)+
  labs(title="Con_Protein")+
  ylab("Log2(Intensity)")+
  stat_compare_means(method = "t.test",paired = T)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "mean_inchip",paired = T)+      
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




#polute_pro
library(reshape2)
library(ggplot2)
library(broom)
library(tidyverse)
pair_line_boxplot_data$pep<-rownames(pair_line_boxplot_data[com,])
voli<-melt(pair_line_boxplot_data)
colnames(voli) <- c("gene", "sample", "expression")
voli$group <- c(rep("inchip",42),rep("outchip",42))
voli$gene<-as.character(voli$gene)
voli$sample<-as.character(voli$sample)


# 使用 tidyverse 和 broom 包计算 p 值
library(broom)
p_values <- voli %>%
  group_by(gene)%>%
  do(tidy(t.test(.$expression ~ .$group)))%>%
  mutate(significance = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01 ~ "**",
    p.value < 0.05 ~ "*",
    TRUE ~ ""
  ))


# 将 p 值添加到原始数据中
voli <- inner_join(voli, p_values, by = "gene")

group<-voli$group
c("KRT7",)
# 使用ggplot2绘制小提琴图
pdf(file="analysis/peptides/hela_hand_pep/20231013_pro_con_VIO.pdf",width = 8,height = 6)

voli %>% ggplot(aes(x=gene, y=log2(expression),fill=group)) + 
  geom_violin(scale="width", adjust=0.5) +
  #geom_point(position = position_jitterdodge()) +
  ggtitle("Plu_Pro_vio")+
  ylab("Log2(Intensity)")+
  scale_fill_manual(values=c("#4682B4" ,"#AF46B4"))+
  geom_text(aes(label = significance, y = 8.5), vjust = 0) +
  #coord_flip()+
  theme_bw()+
  theme(  #legend.position = "none",
          axis.line = element_line(size=1.5),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),# Values for face are one of "plain", "italic", "bold" and "bold.italic"
          axis.text.y = element_text(color="black", size=12, face="plain"),
          axis.title.x = element_text(color="black", size=16, face="bold"),
          axis.title.y = element_text(color="black", size=16, face="bold"),
          legend.text=element_text(color="black", size=16, face="bold"),
          legend.title=element_text(color="black", size=16, face="bold"),
          title = element_text(color="black", size=16, face="bold"))
dev.off()


########################
########用到函数########
########################

median_correct<-function(data){
  median_num<-data.frame(t(apply(data,2,function(x)median(x,na.rm = T))))
  datax<-rbind(median_num,data)
  datay<-data.frame(t(datax))
  data_t<-(datay/datay[,1])*10000
  data_final<-data.frame(t(data_t))
  data_final<-data_final[-1,]
  return(data_final)
}
half_quant<-function(data,x=0.5){
  n<-ncol(data)
  data2<-data[ rowSums ( is.na ( data ) )  < n*x , ] 
  return(data2)
}
overview_umap<-function(data,data2,group=c("TNM_stage"),eclipse_level=0.95,group_values=c("#999999", "#E69F00","#56B4E9"),title=c(""),n_neighbors=1){
  library(umap)
  library(ggplot2)
  library(ggrepel)
  blca_data<-data.frame(t(data))
  group_all<-data2
  blca.label<-group_all[,group]
  blca.label2<-group_all[,"sample"]
  blca.umap = umap(blca_data, n_components = 2, random_state = 15,n_neighbors=n_neighbors) 
  layout <- blca.umap[["layout"]] 
  layout <- data.frame(layout) 
  final <- cbind(layout, blca.label,blca.label2)
  
  p<-ggplot(final, aes(x=X1, y=X2,colour=blca.label)) + 
    geom_point()+
    #labs(title=title)+
    #geom_text_repel(label=blca.label2)+
    xlab("UMAP1")+
    ylab("UMAP2")+
    theme_bw() +
    scale_color_manual(values=a,"Group")+
    stat_ellipse(level = 0.95)+
    theme_bw()+   #空白背景
    theme(  #legend.position = "none",
      axis.line = element_line(size=1.5),
      axis.text.x = element_text(color="black", size=12, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
      axis.text.y = element_text(color="black", size=12, face="plain"),
      axis.title.x = element_text(color="black", size=16, face="bold"),
      axis.title.y = element_text(color="black", size=16, face="bold"),
      legend.text=element_text(color="black", size=16, face="bold"),
      legend.title=element_text(color="black", size=16, face="bold"),
      title = element_text(color="black", size=16, face="bold")
    ) 
  p
  return(p)}
scale_neg1_1<-function(data){
  min_num<-data.frame(t(apply(data,2,function(x)min(x,na.rm = T))))
  bottle_num<-data.frame(t(apply(data,2,function(x)(max(x,na.rm = T)-min(x,na.rm = T)))))
  datax<-rbind(min_num,data)
  datay<-data.frame(t(datax))
  data_t<-datay-datay[,1]
  data_tt<-data.frame(t(data_t))
  data_tt<-data_tt[-1,]
  
  datax2<-rbind(bottle_num,data_tt)
  datay2<-data.frame(t(datax2))
  data_t2<-2*datay2/datay2[,1]-1
  data_tt2<-data.frame(t(data_t2))
  data_tt2<-data_tt2[-1,]
  return(data_tt2)
}


data_trans_binary<-function(data){
  medians <- apply(data, 1, median)
  binary_matrix <- ifelse(data > medians, paste0(rownames(data), "_high"), paste0(rownames(data), "_low"))
  return(binary_matrix)
}

overlap_two_vector<-function(x=c("1","2","3","4","5"),y=c("4","6","5","7","8")){
  library(dplyr)
  commons<-intersect(x,y)
  uni_x<-setdiff(x,y)
  uni_y<-setdiff(y,x)
  cat("俩共有的特征有",length(commons),"个")
  
  cat("前面那个独特的特征有",length(uni_x),"个")
  
  cat("后面那个独特的特征有",length(uni_y),"个")
  
  listx<-list(commons=commons,uni_x=uni_x,uni_y=uni_y)
  
}



#1、split_data
##input group_data/data_express
length(aa)
aa<-cell_specific_data(data_express=data_hela_hand,group_data=hela_hand_sample,p_value=0.05,method="BH",Fc_value=1.5)
hand_up_clus<-clustern_kegg_go_list(aa[[1]]$DEG_up)
Hela_up_clus<-clustern_kegg_go_list(aa[[2]]$DEG_up)

HepG2_up_clus<-clustern_kegg_go_list(aa[[3]]$DEG_up)
hand_up<-aa[[1]]$DEG_up
hand_plot_up_ppi<-plot_stringdb(hand_up,layout = "kk",deg_filter=2,fil_nod = 2)



Hela_up<-aa[[2]]$DEG_up
Hela_plot_up_ppi<-plot_stringdb(Hela_up,layout = "kk",deg_filter=5,fil_nod = 20)

HepG2_up<-aa[[3]]$DEG_up

HepG2_plot_up_ppi<-plot_stringdb(HepG2_up,layout = "kk",deg_filter=5,fil_nod = 5)
pdf(file="./analysis/hela_hand/20231011_hand_up_Clus.pdf",width = 6,height = 5)
hand_plot_up_ppi
dev.off()
pdf(file="./analysis/hela_hand/20231011_Hela_up_Clus.pdf",width = 6,height = 10)
Hela_plot_up_ppi
dev.off()

pdf(file="./analysis/hela_hand/20231011_up_Clus.pdf",width = 6,height = 5)
HepG2_plot_up_ppi
dev.off()
write.csv(hand_up,file="./analysis/hela_hand_l02g2/hand_up.csv")
write.csv(Hela_up,file="./analysis/hela_hand_l02g2/Hela_up.csv")
write.csv(L02_up,file="./analysis/hela_hand_l02g2/L02_up.csv")
write.csv(HepG2G2_up,file="./analysis/hela_hand_l02g2/HepG2G2_up.csv")

cell_specific_data<-function(data_express,group_data,p_value=0.01,method="BH",Fc_value=1.5){
  #data_express<-data_hela_hand
  #group_data<-filter_data_hela_hand
  Group_info<-unique(group_data$Group)
  Group_info[1]
  listn <- vector("list",length(Group_info))
  for (i in 1:length(Group_info)){
    listn[i]<-list(assign(paste("data",Group_info[i],sep = "_"),data_express[,unlist(group_data %>% filter(Group==Group_info[i]) %>% .[,"sample"])]))
  }
  names(listn)<-Group_info
  library(data.table)
  listn1 <- vector("list",length(Group_info))
  for (i in 1:length(Group_info)){
    listn1[i]<-list(assign(paste("compare_data",Group_info[i],sep = "_"),cbind(do.call(cbind, listn[-i]),listn[i])))
  }
  
  listn2 <- vector("list",length(Group_info))
  for (i in 1:length(Group_info)){
    n1=1
    n2=ncol(do.call(cbind, listn[-i]))
    n3=ncol(do.call(cbind, listn[-i]))+1
    n4=length(listn1[[i]])
    nn<-paste("DEG_data",Group_info,sep = "_")
    listn2[i]<-list(assign(nn[i],wil_cox_nopair_test(listn1[[i]],n1=n1,n2=n2,n3=n3,n4=n4,p_value = p_value,method = method,Fc_value = Fc_value)))
    
  }
  names(listn2)<-nn
  
  return(listn2)
}
extrat_up_Genes<-function(cell_specific_data){
  n<-length(cell_specific_data)
  DEG<-
    aa[[i]]$DEG_up
  
  
}

clustern_kegg_go_list<-function(genelist,arg1=c("KEGG"),arg2=c("GO_BP"),arg3=c("GO_CC"),arg4=c("GO_MF")){
  library(clusterProfiler)
  suppressMessages(library(org.Hs.eg.db))
  #library(topGO)
  library(Rgraphviz)
  #library(pathview)
  library(enrichplot)
  library(ggplot2)
  library(DOSE)
  library(stringr)
  library(dplyr)
  library(stats)
  library(patchwork)
  
  
  gene<-data.frame(genelist)
  print(dim(gene)[1])
  colnames(gene)<-"gene"
  DEA<- bitr(gene$gene,fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  DEA<-distinct(DEA,ENTREZID,.keep_all = TRUE)
  row.names(DEA)<-DEA$ENTREZID
  
  #kegg
  KEGG<- enrichKEGG(DEA$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,qvalueCutoff =0.05,pAdjustMethod = "BH")
  if (is.null(KEGG)) KEGG_plot<-list()  else (
    KEGG_plot<-barplot(KEGG,showCategory=15,title=arg1)+scale_y_discrete(labels=function(x) str_wrap(x, width=50))+
      theme_bw(base_size = 16) +
      theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+   #空白背景
      theme(
        
        axis.text.x = element_text(color="black", size=16, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
        axis.text.y = element_text(color="black", size=16, face="plain"),
        axis.title.x = element_text(color="black", size=16, face="plain"),
        axis.title.y = element_text(color="black", size=16, face="plain"),
      ) )
  
  
  #gobp
  GOBP<-enrichGO(DEA$ENTREZID,OrgDb="org.Hs.eg.db", keyType = "ENTREZID",ont = "BP", pvalueCutoff = 0.05,qvalueCutoff = 0.05,pAdjustMethod = "BH",readable = T)
  if (is.null(GOBP)) cut_BP<-list() else(
    cut_BP<- simplify(GOBP, cutoff=0.7, by="p.adjust", select_fun=min))
  if (is.null(GOBP)) GO_BP_plot<-list() else(
    GO_BP_plot<-barplot(cut_BP,showCategory=15,title=arg2)+scale_y_discrete(labels=function(x) str_wrap(x, width=50))+
      theme_bw(base_size = 16) +
      theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+   #空白背景
      theme(
        
        axis.text.x = element_text(color="black", size=16, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
        axis.text.y = element_text(color="black", size=16, face="plain"),
        axis.title.x = element_text(color="black", size=16, face="plain"),
        axis.title.y = element_text(color="black", size=16, face="plain"),
      ) 
  )
  #goCC
  GOCC<-enrichGO(DEA$ENTREZID,OrgDb="org.Hs.eg.db", keyType = "ENTREZID",ont = "CC", pvalueCutoff = 0.05,qvalueCutoff = 0.05,pAdjustMethod = "BH",readable = T)
  if (is.null(GOCC)) cut_CC<-list() else(
    cut_CC<- simplify(GOCC, cutoff=0.7, by="p.adjust", select_fun=min))
  if (is.null(GOCC)) GO_CC_plot<-list() else(
    GO_CC_plot<-barplot(cut_CC,showCategory=15,title=arg3)+scale_y_discrete(labels=function(x) str_wrap(x, width=50))+
      theme_bw(base_size = 16) +
      theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+   #空白背景
      theme(
        
        axis.text.x = element_text(color="black", size=16, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
        axis.text.y = element_text(color="black", size=16, face="plain"),
        axis.title.x = element_text(color="black", size=16, face="plain"),
        axis.title.y = element_text(color="black", size=16, face="plain"),
      ) 
  )
  #goMF
  GOMF<-enrichGO(DEA$ENTREZID,OrgDb="org.Hs.eg.db", keyType = "ENTREZID",ont = "MF", pvalueCutoff = 0.05,qvalueCutoff = 0.05,pAdjustMethod = "BH",readable = T)
  if (is.null(GOMF)) cut_MF<-list() else(
    cut_MF<- simplify(GOMF, cutoff=0.7, by="p.adjust", select_fun=min))
  if (is.null(GOMF)) GO_MF_plot<-list() else(
    GO_MF_plot<-barplot(cut_MF,showCategory=15,title=arg4)+scale_y_discrete(labels=function(x) str_wrap(x, width=50))+
      theme_bw(base_size = 16) +
      theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+   #空白背景
      theme(
        
        axis.text.x = element_text(color="black", size=16, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
        axis.text.y = element_text(color="black", size=16, face="plain"),
        axis.title.x = element_text(color="black", size=16, face="plain"),
        axis.title.y = element_text(color="black", size=16, face="plain"),
      ) 
  )
  PLOT<-list(KEGG_plot,GO_BP_plot,GO_CC_plot,GO_MF_plot)
  names(PLOT)<-c("KEGG_plot","GO_BP_plot","GO_CC_plot","GO_MF_plot")
  dataframes<-list(KEGG,GOBP,GOCC,GOMF)
  names(dataframes)<-c("KEGG","GOBP","GOCC","GOMF")
  cut_dataframes<-list(cut_BP,cut_CC,cut_MF)
  names(cut_dataframes)<-c("cut_BP","cut_CC","cut_MF")
  
  ALL<-list(PLOT,dataframes,cut_dataframes)
  names(ALL)<-c("richplot","richterm","cut_richterm")
  print(dim(gene))
  return(ALL)
  
}
plot_stringdb<-function(ge,deg_filter=5,fil_nod=4,layout=c('star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl', 'lgl')){
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(STRINGdb)
  library(igraph)
  library(ggraph)
  # 创建STRINGdb对象
  string_db <- STRINGdb$new( version="11.5", species=9606,
                             score_threshold=400, input_directory="./data/")
  
  # 将Gene Symbol转换为Entrez ID
  gene <- ge %>% clusterProfiler::bitr(fromType = "SYMBOL", 
                                       toType = "ENTREZID", 
                                       OrgDb = "org.Hs.eg.db", 
                                       drop = T)
  data_mapped <- gene %>% string_db$map(my_data_frame_id_col_names = "SYMBOL", 
                                        removeUnmappedRows = TRUE)
  #string_db$plot_network( data_mapped$STRING_id )
  data_links <- data_mapped$STRING_id%>% string_db$get_interactions()
  #data_links %>% distinct(Sepal.Length, Petal.Width, .keep_all = TRUE)
  #使用igraph和ggraph可视化蛋白互作网络图
  # 转换stringID为Symbol，只取前两列和最后一列
  links <- data_links %>% distinct( from,to,combined_score,.keep_all = TRUE) %>%
    dplyr::mutate(from = data_mapped[match(from, data_mapped$STRING_id), "SYMBOL"]) %>% 
    dplyr::mutate(to = data_mapped[match(to, data_mapped$STRING_id), "SYMBOL"]) %>%  
    dplyr::select(from, to , last_col()) %>% 
    dplyr::rename(weight = combined_score)
  # 节点数据
  #nodes <- links %>% { data.frame(gene = c(.$from, .$to)) } %>% dplyr::distinct()
  # 创建网络图
  # 根据links和nodes创建
  #net <- igraph::graph_from_data_frame(d=links,vertices=nodes,directed = F)
  # 添加一些参数信息用于后续绘图
  # V和E是igraph包的函数，分别用于修改网络图的节点（nodes）和连线(links)
  #igraph::V(net)$deg <- igraph::degree(net) # 每个节点连接的节点数
  #igraph::V(net)$size <- igraph::degree(net)/5 #
  #igraph::E(net)$width <- igraph::E(net)$weight/10
  # 去除游离的互作关系
  # 如果links数据框的一个link的from只出现过一次，同时to也只出现一次，则将其去除
  links_2 <- links %>% dplyr::mutate(from_c = dplyr::count(., from)$n[match(from, dplyr::count(., from)$from)]) %>%
    dplyr::mutate(to_c = dplyr::count(., to)$n[match(to, dplyr::count(., to)$to)]) %>%
    dplyr::filter(!(from_c <= fil_nod & to_c <= fil_nod)) %>%
    dplyr::select(1,2,3)
  #fil_nod=5
  # 新的节点数据
  nodes_2 <- links_2 %>% { data.frame(gene = c(.$from, .$to)) } %>% dplyr::distinct()
  # 创建网络图
  net_2 <- igraph::graph_from_data_frame(d=links_2,vertices=nodes_2,directed = F)
  # 添加必要的参数
  igraph::V(net_2)$deg <- igraph::degree(net_2)
  igraph::V(net_2)$size <- igraph::degree(net_2)/5
  igraph::E(net_2)$width <- igraph::E(net_2)$weight/10
  #layout="kk"
  # deg_filter=3
  P<-ggraph(net_2,layout = layout)+
    geom_edge_fan(aes(edge_width=width), color = "lightblue", show.legend = F)+
    geom_node_point(aes(size=size), color="orange", alpha=0.6)+
    geom_node_text(aes(filter=deg>deg_filter, label=name), size = 5, repel = T)+
    scale_edge_width(range = c(0.2,1))+
    scale_size_continuous(range = c(1,10) )+
    guides(size=F)+
    theme_graph()
  P
  return(P)
}

pep_GRVAY_score_calc<-function(input="AAAAKDSADF"){
  library(Biostrings)
  library(dplyr)
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
  results <- calculate_protein_hydrophobicity(input)
  
  return(results)
}
pep_miscleva_num_calc<-function(input="AAAKDSADKR"){
  K_num=length(strsplit(input, '')[[1]][strsplit(input, '')[[1]]=="K"])
  R_num=length(strsplit(input, '')[[1]][strsplit(input, '')[[1]]=="R"])
  
 if( strsplit(input, '')[[1]][length(strsplit(input, '')[[1]])] == "K"|strsplit(input, '')[[1]][length(strsplit(input, '')[[1]])] == "R"){
   Miscleavage=K_num+R_num-1
 } else{
   Miscleavage=K_num+R_num
 }
  return(Miscleavage)
}
pep_miscleva_num_calc()
