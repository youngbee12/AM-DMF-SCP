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
library(Seurat)

pbmc<-Read10X('')


g2m_genes<-cc.genes$g2m.genes
g2m_genes<-CaseMatch(search = g2m_genes,match = rownames(data_hela_a549))
s_genes<-cc.genes$s.genes
s_genes<-CaseMatch(search = s_genes,match = rownames(data_hela_a549))
pbmc <- NormalizeData(dd)
cc_info<-CellCycleScoring(pbmc,s.features = s_genes,g2m.features = g2m_genes)
cc_info@meta.data
dd<-CreateSeuratObject(data_hela_a549,project = "test",min.cells = 3,min.features = 300)
library(GSVA)
gsva_out<-gsva(scale(as.matrix(data_hela_a549)), list(S=s_genes,G2M=g2m_genes),method="gsva", min.sz = 1,max.sz=Inf,verbose=T,parallel.sz=1)

RidgePlot(pbmc,features = c())

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
overview_umap<-function(data,data2,group=c("TNM_stage"),eclipse_level=0.95,group_values=c("#999999", "#E69F00","#56B4E9"),title=c("")){
  library(umap)
  library(ggplot2)
  library(ggrepel)
  blca_data<-data.frame(t(data))
  group_all<-data2
  blca.label<-group_all[,group]
  blca.label2<-group_all[,"sample"]
  blca.umap = umap(blca_data, n_components = 2, random_state = 15) 
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
aa<-cell_specific_data(data_express=data_hela_a549,group_data=hela_a549_sample,p_value=0.05,method="BH",Fc_value=1.5)
A549_up_clus<-clustern_kegg_go_list(aa[[1]]$DEG_up)
Hela_up_clus<-clustern_kegg_go_list(aa[[2]]$DEG_up)

HepG2_up_clus<-clustern_kegg_go_list(aa[[3]]$DEG_up)
A549_up<-aa[[1]]$DEG_up
A549_plot_up_ppi<-plot_stringdb(A549_up,layout = "kk",deg_filter=2,fil_nod = 2)



Hela_up<-aa[[2]]$DEG_up
Hela_plot_up_ppi<-plot_stringdb(Hela_up,layout = "kk",deg_filter=5,fil_nod = 20)

HepG2_up<-aa[[3]]$DEG_up

HepG2_plot_up_ppi<-plot_stringdb(HepG2_up,layout = "kk",deg_filter=5,fil_nod = 5)
pdf(file="./analysis/hela_a549_Hepg2/20231011_A549_up_Clus.pdf",width = 6,height = 5)
A549_plot_up_ppi
dev.off()
pdf(file="./analysis/hela_a549_Hepg2/20231011_Hela_up_Clus.pdf",width = 6,height = 10)
Hela_plot_up_ppi
dev.off()

pdf(file="./analysis/hela_a549_Hepg2/20231011_HepG2_up_Clus.pdf",width = 6,height = 5)
HepG2_plot_up_ppi
dev.off()
write.csv(A549_up,file="./analysis/hela_a549_hepg2/A549_up.csv")
write.csv(Hela_up,file="./analysis/hela_a549_hepg2/Hela_up.csv")

write.csv(HepG2_up,file="./analysis/hela_a549_hepg2/HepG2_up.csv")

cell_specific_data<-function(data_express,group_data,p_value=0.01,method="BH",Fc_value=1.5){
#data_express<-data_hela_a549
#group_data<-filter_data_hela_a549
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


