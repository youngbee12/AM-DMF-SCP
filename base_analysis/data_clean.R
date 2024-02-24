########################
######导入数据##########
########################

load(file="./data/20231011_scdata_pre.Rdata")
library(tidyverse)
library(colortools)
cor_res<-data.frame(cor(log2(sc_data),method="pearson",use="complete.obs"))





Group_info<-unique(sc_group$Group)
Group_info[1]
listn <- vector("list",length(Group_info))
for (i in 1:length(Group_info)){
  listn[i]<-list(assign(paste("data",Group_info[i],sep = "_"),cor_res[unlist(sc_group %>% filter(Group==Group_info[i]) %>% .[,"sample"]),unlist(sc_group %>% filter(Group==Group_info[i]) %>% .[,"sample"])]))
}
names(listn)<-Group_info
listn <- lapply(listn, function(df) {
  df[df == 1] <- NA
  return(df)
})

mean_rows_list <- lapply(listn, function(df) {
  row_means <- rowMeans(df, na.rm = T)
  return(row_means)
})
del_sample<-c("A549_16")
#del_sample<-c("A549_2","A549_11","HepG2_17","A549_16")
sc_data<-sc_data[,!colnames(sc_data)%in%del_sample]
sc_group <-sc_group[colnames(sc_data),]
data_30<-filter_na_bygroup(sc_data,group_info = dplyr::select(sc_group,c("sample","Group")),x=0.8)
#data_30<-half_quant(sc_data,x=0.)
#half_data_meadian_impute<-log2(data_30)
#Knn补缺
library(impute)
half_data_meadian_impute<-impute.knn(as.matrix(log2(data_30)) ,k = 3, rowmax = 0.7, colmax = 0.8, maxp = 1500, rng.seed=362436069)
half_data_meadian_impute<-data.frame(half_data_meadian_impute$data)
colnames(half_data_meadian_impute)<-colnames(data_30)


#PCA
library(FactoMineR)
library(factoextra)
res.pca1 = PCA(t(half_data_meadian_impute), scale.unit = TRUE, ncp = 8,  graph = T,axes = 1:3)
pca1<-fviz_pca_ind(res.pca1,
                   label = "none",
                   habillage = factor(sc_group$Group),
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
#UMAP
overview_umap(data = half_data_meadian_impute,data2 =sc_group,group = c("Group"),eclipse_level = 0,group_values = a)
a<-wheel("steelblue", num = 8)
sc_data<-half_data_meadian_impute
save(sc_data,sc_group,file="./data/20231011_data_clean_sc_filterdatabygroup_80.Rdata")
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
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+   #空白背景
    theme(
      
      axis.text.x = element_text(color="black", size=16, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
      axis.text.y = element_text(color="black", size=16, face="plain"),
      axis.title.x = element_text(color="black", size=16, face="bold"),
      axis.title.y = element_text(color="black", size=16, face="bold"),
      legend.text=element_text(color="black", size=16, face="plain"),
      legend.title=element_text(color="black", size=16, face="bold"),
      title = element_text(color="black", size=16, face="bold")
    ) 
    p
  return(p)
  
}
filter_na_bygroup<-function(data,group_info,x=0.5){
  colnames(group_info)<-c("sample","group")
  n<-length(unique(group_info$group))
  listn<-list()
  #data=data_total
  for (i in 1:n){
    t_sm<-group_info[group_info$group==unique(group_info$group)[i],]$sample
    listn[i]<-list(data[,colnames(data) %in% t_sm] %>%
                     .[ rowSums ( is.na ( . ) )  < ncol(.)*x,] %>% rownames(.))
  }
  aa<-Reduce(intersect,listn)
  data_final<-data[aa,]
  return(data_final)
}
