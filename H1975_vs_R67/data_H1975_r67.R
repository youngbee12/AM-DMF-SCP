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
load("./data/20231011_data_clean_sc_filterdatabygroup_80.Rdata")

naiyao_sample<-sc_group %>% arrange(Group)%>%filter(Group=="H1975"|Group=="R67") 
data_naiyao <-sc_data[,naiyao_sample$sample]



load(file="./data/20231011_scdata_pre.Rdata")
sc_data_raw<-sc_data





naiyao_sample<-sc_group %>% arrange(Group)%>%filter(Group=="H1975"|Group=="R67") 
data_naiyao <-sc_data[,naiyao_sample$sample]





listn<-list()
#data=data_total
n<-length(unique(naiyao_sample$Group))
for (i in 1:n){
 t_sm<-naiyao_sample[naiyao_sample$Group==unique(naiyao_sample$Group)[i],]$sample
listn[i]<-list(data_naiyao[,colnames(data_naiyao) %in% t_sm] %>%
                 .[ rowSums ( is.na ( . ) )  < ncol(.)*0.5,] %>% rownames(.))
}
#overlap_two_vector(listn[[1]],listn[[2]])
data_50<-filter_na_bygroup(data_naiyao,group_info = dplyr::select(naiyao_sample,c("sample","Group")),x=0.5)
#Knn补缺
library(impute)
half_data_meadian_impute<-impute.knn(as.matrix(log2(data_50)) ,k = 3, rowmax = 0.7, colmax = 0.8, maxp = 1500, rng.seed=362436069)
half_data_meadian_impute<-data.frame(half_data_meadian_impute$data)
colnames(half_data_meadian_impute)<-colnames(data_50)
data_naiyao<-half_data_meadian_impute

pdf(file="./analysis/H1975_vs_R67/20231016_umap_H1975_vs_R67.pdf",height = 5,width = 6 )
overview_umap(data = data_naiyao,data2 =naiyao_sample,group = c("Group"),eclipse_level = 0.95,group_values = wheel("steelblue", num = 4),n_neighbors = 15)
dev.off()
res.pca4 = PCA(t(data_naiyao), scale.unit = TRUE, ncp = 8,  graph = T,axes = 1:3)
pca4<-fviz_pca_ind(res.pca4,
                   label = "none",
                   habillage = factor(naiyao_sample$Group),
                   #palette =c("#989898", "#E69F00","#56B4E9") ,
                   addEllipses = T,
                   ellipse.level=0,
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
wide_df_naiyao<-wide_df %>% filter(Group=="H1975"|Group=="R67")
filter_data_naiyao<-sc_group %>%filter(Group=="H1975"|Group=="R67")
pdf(file="./analysis/H1975_vs_R67/20231009_identi_naiyao.pdf",height = 5,width = 5 )
ggplot(wide_df_naiyao, aes(x=Group, y=mean,color=Group,fill=Group)) +
  geom_bar(position=position_dodge(), stat="identity",colour='white',width=.6,alpha=0.8) +
  # scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.4)+
  scale_color_manual(values=c("#4682B4","#AF46B4"))+
  scale_fill_manual(values=c("#4682B4","#AF46B4"))+
  geom_jitter(data = filter_data_naiyao,aes(x=Group,y=p_num),width=.1)+
  scale_y_continuous(limits=c(0,3100),breaks=seq(0,3000,500))+
  geom_text(aes(label = ceiling(mean)), vjust = -9, colour = "black")+
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

pdf(file = "./analysis/H1975_vs_R67/20231011_corr_H1975_vs_R67.pdf",height = 5,width = 6)


corrplot(cor(log2(sc_data[,naiyao_sample$sample]),method="pearson",use="complete.obs"),
         method="color",tl.col = "black",
         
         col=c(rep("white", 500),rep('grey', 400) ,COL2('RdYlBu', 200)),
         col.lim=c(0.7,1),
         #addCoef.col = "white",
         #type="upper")
)
dev.off()








#PCA载荷
# 运行PCA
pca <- prcomp(t(data_naiyao))

# 查看每个主成分的贡献度
pca$explained_variance_ratio_
# 获取PC1的载荷
loadings <- data.frame(pca$rotation[, 2])





#Umap 载荷



group<-factor(paste0("group",arrange(filter_data_naiyao,Group)$Group))
design<-model.matrix(~0+group)
colnames(design)<-levels(group)
library(limma)
#创建对比矩阵
group_names <-levels(group)
contrasts_names<-combn(group_names,2,paste,collapse="-")
contrast_formulas<-combn(group_names,2,function(x)(paste0(x[1],"-",x[2])))
contrast_matrix<-makeContrasts(contrasts = setNames(contrast_formulas,contrasts_names),levels = design)
fit <-lmFit(data_naiyao,design)
fit2<-contrasts.fit(fit,contrast_matrix)
fit2<-eBayes(fit2)

results<-list()
for(contrast in colnames(contrast_matrix)){
  results[[contrast]] <-topTable(fit2,coef = contrast,number = 200)
  
}
# 创建空列表用于保存所有的基因
all_genes <- list()

# 遍历所有的对比结果
for (contrast in names(results)) {
  # 从对比结果中提取显著差异表达的基因
  significant_genes <- rownames(results[[contrast]][results[[contrast]]$adj.P.Val < 0.05,])
  
  # 将这些基因加入到列表中
  all_genes[[contrast]] <- significant_genes
}

# 创建一个全局基因列表，包含所有对比中出现的基因
heat_gene_list <- unique(unlist(all_genes))
data_heat<-data_naiyao[c(aa[[1]]$DEG_up,aa[[2]]$DEG_up,aa[[3]]$DEG_up),]
pdf(file = "./analysis/H1975_vs_R67/20231012_HEATMAP_hela_a549.pdf",height = 5,width = 6)

anno_df = data.frame(
  
  mod = c(rep("A549",length(aa[[1]]$DEG_up)),rep("Hela",length(aa[[2]]$DEG_up)),rep("HepG2",length(aa[[3]]$DEG_up)))
)
ha = rowAnnotation(df = anno_df,
                   col = list(Group=c("A549"=wheel("steelblue", num = 4)[1],
                                      "Hela"=wheel("steelblue", num = 4)[2],
                                      
                                      "HepG2"=wheel("steelblue", num = 4)[3])))
Heatmap(
  t(scale_neg1_1(t(data_heat))), 
  name = "mat",
  #col = col_fun,
  cluster_columns = F,
  cluster_rows=F,
  show_row_names = F,
  show_column_names = F,
  show_row_dend =F,
  show_column_dend = F,
  #row_km = 4,
  column_split = arrange(naiyao_sample,Group)$Group,
  row_names_gp = gpar(fontsize = 10),
  right_annotation = ha,
  top_annotation = HeatmapAnnotation(df=naiyao_sample %>% arrange(Group)%>%
                                       dplyr::select(Group),
                                     col = list(Group=c("A549"=wheel("steelblue", num = 4)[1],
                                                        "Hela"=wheel("steelblue", num = 4)[2],
                                                        "LO2"=wheel("steelblue", num = 4)[3],
                                                        "HepG2"=wheel("steelblue", num = 4)[4]))
  )
)
dev.off()

write.table(protein_data,file="20231016_sc_Data.txt",sep = "\t")









#####ssGSEA
###hallmark,KEGG,reactome
library(msigdbr)
hall <- msigdbr(species = "Homo sapiens", category = "H")
#kegg <-msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:KEGG")
reactome <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:REACTOME")

merge_term<-rbind(hall,reactome)

merge_term_list <- merge_term %>%
  split(.$gs_name) %>%
  lapply(function(df) {
    df$gene_symbol
  })

names(merge_term_list) <- merge_term$gs_name[!duplicated(merge_term$gs_name)]
########################
######执行ssGSEA########
########################
library(GSVA)
mat<-scale(as.matrix(data_naiyao))
gsva_out_sc<-gsva(mat, merge_term_list,method="gsva", min.sz = 5,max.sz=Inf,verbose=T,parallel.sz=1)
DD<-data.frame(gsva_out_sc)
DD1<-log2(DD)
DEP_sc_ssgsea<-wil_cox_nopair_test(data.frame(gsva_out_sc),n1=1,n2=18,n3=19,n4=36,Fc_value = 0.1,method = "none",p_value = 0.01)
vocanol_DEP_sc_ssGSEA<-vocanol_plot(DEP_sc_ssgsea$data,p.value = 0.01,Fc_value = 0.1)+
  xlab("Protein mean( Log2 fold-change, R67 vs H1975)")+
  xlim(-4,5)



aaa<-msigdb_GSEA(DEP_sc,category = "H")
bbb<-msigdb_GSEA(DEP_total,category = "H")
pdf(file = "./analysis/H1975_vs_R67/20231023_EMT_GSEA_TOTAL.pdf",height = 4,width = 6)
gseaplot2(bbb, geneSetID = 4,title = bbb$Description[4])
dev.off()
pdf(file = "./analysis/H1975_vs_R67/20231023_EMT_GSEA_sc.pdf",height = 4,width = 6)
gseaplot2(aaa, geneSetID = 1,title = aaa$Description[1])
dev.off()
AAA <- setReadable(aaa, 'org.Hs.eg.db', 'ENTREZID')
pdf(file = "./analysis/H1975_vs_R67/20231023_GSEA_sc.pdf",height = 8,width = 5)
cnetplot(AAA, categorySize="pvalue", foldChange=geneList3)

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
overview_umap<-function(data,data2,group=c("TNM_stage"),eclipse_level=0.95,group_values=c("#999999", "#E69F00","#56B4E9"),title=c(""),n_neighbors = 15){
  library(umap)
  library(ggplot2)
  library(ggrepel)
  blca_data<-data.frame(t(data_na))
  group_all<-data2
  blca.label<-group_all[,group]
  blca.label2<-group_all[,"sample"]
  blca.umap = umap(blca_data, n_components = 2, random_state = 15,n_neighbors = n_neighbors) 
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
  p<-list(p,final)
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
aa<-cell_specific_data(data_express=data_naiyao,group_data=naiyao_sample,p_value=0.01,method="BH",Fc_value=2)
A549_up_clus<-clustern_kegg_go_list(aa[[1]]$DEG_up)
Hela_up_clus<-clustern_kegg_go_list(aa[[2]]$DEG_up)

HepG2_up_clus<-clustern_kegg_go_list(aa[[3]]$DEG_up)
A549_up<-aa[[1]]$DEG_up
A549_plot_up_ppi<-plot_stringdb(A549_up,layout = "kk",deg_filter=0,fil_nod = 0)



Hela_up<-aa[[2]]$DEG_up
Hela_plot_up_ppi<-plot_stringdb(Hela_up,layout = "kk",deg_filter=0,fil_nod = 0)

HepG2_up<-aa[[3]]$DEG_up

HepG2_plot_up_ppi<-plot_stringdb(HepG2_up,layout = "kk",deg_filter=1,fil_nod = 1)
pdf(file="./analysis/H1975_vs_R67/20231016_A549_up_Clus.pdf",width = 9,height = 7.5)
A549_plot_up_ppi
dev.off()
pdf(file="./analysis/H1975_vs_R67/20231016_Hela_up_Clus.pdf",width = 9,height = 7.5)
Hela_plot_up_ppi
dev.off()

pdf(file="./analysis/H1975_vs_R67/20231016_HepG2_up_Clus.pdf",width = 9,height = 7.5)
HepG2_plot_up_ppi
dev.off()
write.csv(A549_up,file="./analysis/hela_a549_l02_HepG2g2/A549_up.csv")
write.csv(Hela_up,file="./analysis/hela_a549_l02_HepG2g2/Hela_up.csv")
write.csv(L02_up,file="./analysis/hela_a549_l02_HepG2g2/L02_up.csv")
write.csv(HepG2G2_up,file="./analysis/hela_a549_l02_HepG2g2/HepG2G2_up.csv")

cell_specific_data<-function(data_express,group_data,p_value=0.01,method="BH",Fc_value=1.5){
  #data_express<-data_naiyao
  #group_data<-filter_data_naiyao
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
  gene <- ge %>% clusterProfilƒer::bitr(fromType = "SYMBOL", 
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
msigdb_GSEA<-function(DEP_list,category = "H"){
 #DEP_list<-DEP_sc
  library("msigdbr")
  m <- msigdbr(species = "Homo sapiens", category = category)
  geneset2 <- m[,c(3,5)]
  geneList <- DEP_list$data%>%
    dplyr::select(Log2FC) %>%
    dplyr::arrange(-Log2FC) %>% 
    dplyr::mutate(gene = rownames(.)) %>%
    dplyr::select(gene,Log2FC) %>%
    filter(gene %in% c(DEP_list$DEG_up,DEP_list$DEG_down))
  DEA_gsea<- bitr(geneList$gene,fromType="SYMBOL", toType="ENTREZID",OrgDb="org.Hs.eg.db")
  DEA_gsea<-distinct(DEA_gsea,ENTREZID,.keep_all = TRUE)
  names(DEA_gsea)[1]<-"gene"
  geneList1<-filter(geneList,gene %in% DEA_gsea$gene)
  me3<-merge(geneList1,DEA_gsea,by="gene",all = T)
  geneList2<-me3[,2]
  names(geneList2) = as.character(me3[,3])
  geneList3 = sort(geneList2, decreasing = TRUE)
 # p1 <- cnetplot(edox, foldChange=geneList3)
  #p2 <- cnetplot(edox, foldChange=geneList3, circular = TRUE, colorEdge = TRUE) 
  #p3 <- heatplot(edox, showCategory=5)
  #p4 <- heatplot(edox, foldChange=geneList3, showCategory=5)
  #edox2 <- pairwise_termsim(edox)
  #p1 <- treeplot(edox2)
  #p2 <- treeplot(edox2, hclust_method = "average")
  #aplot::plot_list(p1, p2, tag_levels='A')
  #cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])
 #   edox <- setReadable(aaa, 'org.Hs.eg.db', 'ENTREZID')
  ABCD <- GSEA(geneList = geneList3,
               TERM2GENE = geneset2,
               minGSSize = 1, maxGSSize = Inf,
               verbose = T,
               eps=0,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               seed = T  )
  return(ABCD)
}

