########################
######导入数据##########
########################
library(dplyr)
data_system <- sc_group %>% arrange(Group)%>%filter(Group=="A549"|Group=="Hela"|Group=="H1975"|Group=="R67") 

group_system<-group_system%>% mutate(batch=ifelse(batch=="batch3","Chip2","Chip1"))

data_system_0 <-sc_data_raw[,group_system$sample]
listn<-list()

n<-length(unique(group_system$Group))
for (i in 1:n){
  t_sm<-group_system[group_system$Group==unique(group_system$Group)[i],]$sample
  listn[i]<-list(data_system[,colnames(data_system) %in% t_sm] %>%
                   .[ rowSums ( is.na ( . ) )  < ncol(.)*0.5,] %>% rownames(.))
}
data_50<-filter_na_bygroup(data_system,group_info = dplyr::select(group_system,c("sample","Group")),x=0.5)

#Knn补缺
library(impute)
half_data_meadian_impute<-impute.knn(as.matrix(log2(data_50)) ,k = 5, rowmax = 0.7, colmax = 0.8, maxp = 1500, rng.seed=362436069)
half_data_meadian_impute<-data.frame(half_data_meadian_impute$data)
colnames(half_data_meadian_impute)<-colnames(data_50)
data_system<-half_data_meadian_impute

bks<-seq(-1,1,length.out=10000)
pdf("")
pdf(file="analysis/system/20231031_sys_heat.pdf",width = 20,height = 5)
sys_heat
dev.off()
sys_heat<-pheatmap(data_system,
                      scale = "row",
                      show_rownames=F,
                      show_colnames=F,
                      breaks=bks,
                      clustering_distance_cols = "correlation",
                      clustering_method = "ward.D2",
                      cluster_rows = T,
                      cluster_cols = T,
                      cellwidth =12, 
                      cellheight = 0,
                      display_numbers = FALSE,
                      annotation_col = group_system%>%select(batch,Group),
                      annotation_colors = list(batch=c("batch1"="black","batch2"="grey"),
                                               Group=c("A549"="#4682B4","Hela"="#8E46B4",
                                                        "H1975"="#B44656","R67"="#AEB446")),
                      fontsize =8,
         na_col = "grey",
                      color = colorRampPalette(c("white","white","white"))(10000),
                      #gaps_row =c(500,1000,3000),
      
)
sys_tree<-sys_heat$tree_col
sys_cluster <- data.frame(cutree(sys_tree,8))
sys_cluster$cutree.sys_tree..8.<-as.character(sys_cluster$cutree.sys_tree..8.)
colnames(sys_cluster)<-"col_cluster"
sys_cluster$sample<-rownames(sys_cluster)
ff<-left_join(sys_cluster,group_system)
head(ff)
ff1<-ff %>% filter(Group=="A549")
ff2<-ff %>% filter(Group=="Hela")
ff3<-ff %>% filter(Group=="H1975")
ff4<-ff %>% filter(Group=="R67")
fisher.test(ff2$col_cluster,ff2$batch)$p.value
pdf(file="analysis/system/20231031_A549_stack.pdf",width = 6,height = 5)
stack_bar_plot_twocluster(ff%>%filter(Group=="A549")%>% select(col_cluster,batch))+
  xlab("Cluster")
dev.off()

pdf(file="analysis/system/20231031_Hela_stack.pdf",width = 6,height = 5)
stack_bar_plot_twocluster(ff%>%filter(Group=="Hela")%>% select(col_cluster,batch))+
  xlab("Cluster")
dev.off()

pdf(file="analysis/system/20231031_H1975_stack.pdf",width = 6,height = 5)
stack_bar_plot_twocluster(ff%>%filter(Group=="H1975")%>% select(col_cluster,batch))+
  xlab("Cluster")
dev.off()

pdf(file="analysis/system/20231031_r67_stack.pdf",width = 6,height = 5)
stack_bar_plot_twocluster(ff%>%filter(Group=="R67")%>% select(col_cluster,batch))+
  xlab("Cluster")
dev.off()




#CV
library(rknn)

cv_sys_A549<-apply(data_system_0[listn[[1]],1:18],1,function(x)cv.coef(as.numeric(x)))
cv_sys_Hela<-apply(data_system_0[listn[[2]],19:36],1,function(x)cv.coef(as.numeric(x)))
cv_sys_H1975<-apply(data_system_0[listn[[3]],37:54],1,function(x)cv.coef(as.numeric(x)))
cv_sys_r67<-apply(data_system_0[listn[[4]],55:72],1,function(x)cv.coef(as.numeric(x)))
median(cv_sys_A549)
median(cv_sys_Hela)
median(cv_sys_H1975)
median(cv_sys_r67)

library(ggplot2)
library(hexbin)

# Assuming you have a data frame 'df' with two numeric variables 'x' and 'y'
 df <- data.frame(x = rnorm(10000), y = rnorm(10000))

cv_data_A549<-data.frame(CV=cv_sys_A549,
                         mean=apply(data_system_0[listn[[1]],1:18],1,function(x)mean(x,na.rm=T))
                                    ) 
cv_data_Hela<-data.frame(CV=cv_sys_Hela,
                         mean=apply(data_system_0[listn[[2]],19:36],1,function(x)mean(x,na.rm=T))
) 
cv_data_H1975<-data.frame(CV=cv_sys_H1975,
                         mean=apply(data_system_0[listn[[3]],37:54],1,function(x)mean(x,na.rm=T))
) 
cv_data_r67<-data.frame(CV=cv_sys_r67,
                         mean=apply(data_system_0[listn[[4]],55:72],1,function(x)mean(x,na.rm=T))
) 
length(cv_data_r67$CV)
cv_all_data<-rbind(cv_data_A549,cv_data_Hela,cv_data_H1975,cv_data_r67)
cv_all_data$Group<-c(rep("A549",1957),
                     rep("NCI-H1975",2181),
                     rep("HeLa",2196),
                     rep("67R",2070))
cv_all_data$Group<-factor(cv_all_data$Group,levels = c("NCI-H1975","67R","HeLa","A549"))

pdf(file="analysis/system/20231126_sys_cv.pdf",width = 16,height = 5)
ggplot(cv_all_data, aes(x = CV, y = log(mean))) + 
  facet_wrap(~ Group, scales = "free_x", nrow = 1) +
  
  geom_hex(bins = 30) + # Adjust the number of bins
  scale_fill_gradient(low = "#B0E0E0", high = "navy") + # Change the color gradient
  xlim(0,0.5)+
  xlab("CV")+
  ylab("Log2(Intensity)")+
  #ylim(6,16)+
  geom_smooth(method="gam")+
  theme_bw() +
  
  theme(  #legend.position = "none",
    #axis.line = element_line(size=1.5),
    #axis.text.x = element_text(color="black", size=12, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
    axis.text.x = element_text(color="black", size=16, face="plain") ,
    axis.text.y = element_text(color="black", size=12, face="plain"),
    axis.title.x = element_text(color="black", size=16, face="bold"),
    axis.title.y = element_text(color="black", size=16, face="bold"),
    legend.text=element_text(color="black", size=16, face="bold"),
    legend.title=element_text(color="black", size=16, face="bold"),
    title = element_text(color="black", size=16, face="bold"),
    strip.text = element_text(color="black", size=16, face="plain"))
dev.off()





ggplot(cv_data_A549, aes(x = CV, y = log(mean))) + 
  geom_hex(bins = 30) + # Adjust the number of bins
  scale_fill_gradient(low = "blue", high = "red") + # Change the color gradient
  xlim(0,0.5)+
  ylim(6,16)+
  theme_bw() +
  
  theme(  #legend.position = "none",
          #axis.line = element_line(size=1.5),
          #axis.text.x = element_text(color="black", size=12, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
          axis.text.x = element_text(color="black", size=16, face="plain") ,
          axis.text.y = element_text(color="black", size=12, face="plain"),
          axis.title.x = element_text(color="black", size=16, face="bold"),
          axis.title.y = element_text(color="black", size=16, face="bold"),
          legend.text=element_text(color="black", size=16, face="bold"),
          legend.title=element_text(color="black", size=16, face="bold"),
          title = element_text(color="black", size=16, face="bold"))

ggplot(cv_data_Hela, aes(x = CV, y = log(mean))) + 
  geom_hex(bins = 30) + # Adjust the number of bins
  scale_fill_gradient(low = "blue", high = "red") + # Change the color gradient
  xlim(0,0.5)+
  ylim(6,16)+
  theme_bw() +
  
  theme(  #legend.position = "none",
    #axis.line = element_line(size=1.5),
    #axis.text.x = element_text(color="black", size=12, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
    axis.text.x = element_text(color="black", size=16, face="plain") ,
    axis.text.y = element_text(color="black", size=12, face="plain"),
    axis.title.x = element_text(color="black", size=16, face="bold"),
    axis.title.y = element_text(color="black", size=16, face="bold"),
    legend.text=element_text(color="black", size=16, face="bold"),
    legend.title=element_text(color="black", size=16, face="bold"),
    title = element_text(color="black", size=16, face="bold"))
ggplot(cv_data_H1975, aes(x = CV, y = log(mean))) + 
  geom_hex(bins = 30) + # Adjust the number of bins
  scale_fill_gradient(low = "blue", high = "red") + # Change the color gradient
  xlim(0,0.5)+
  ylim(6,16)+
  theme_bw() +
  
  theme(  #legend.position = "none",
    #axis.line = element_line(size=1.5),
    #axis.text.x = element_text(color="black", size=12, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
    axis.text.x = element_text(color="black", size=16, face="plain") ,
    axis.text.y = element_text(color="black", size=12, face="plain"),
    axis.title.x = element_text(color="black", size=16, face="bold"),
    axis.title.y = element_text(color="black", size=16, face="bold"),
    legend.text=element_text(color="black", size=16, face="bold"),
    legend.title=element_text(color="black", size=16, face="bold"),
    title = element_text(color="black", size=16, face="bold"))
ggplot(cv_data_r67, aes(x = CV, y = log(mean))) + 
  geom_hex(bins = 30) + # Adjust the number of bins
  scale_fill_gradient(low = "blue", high = "red") + # Change the color gradient
  xlim(0,0.5)+
  ylim(6,16)+
  geom_smooth()+
  theme_bw() +
  
  theme(  #legend.position = "none",
    #axis.line = element_line(size=1.5),
    #axis.text.x = element_text(color="black", size=12, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
    axis.text.x = element_text(color="black", size=16, face="plain") ,
    axis.text.y = element_text(color="black", size=12, face="plain"),
    axis.title.x = element_text(color="black", size=16, face="bold"),
    axis.title.y = element_text(color="black", size=16, face="bold"),
    legend.text=element_text(color="black", size=16, face="bold"),
    legend.title=element_text(color="black", size=16, face="bold"),
    title = element_text(color="black", size=16, face="bold"))


#overlap
##A549
sys_A549_data<-data_system_0[,1:18]
group_system_a549<-group_system %>% filter(Group=="A549")
listn_A549<-list()
n<-length(unique(group_system_a549$batch))
for (i in 1:n){
  t_sm<-group_system_a549[group_system_a549$batch==unique(group_system_a549$batch)[i],]$sample
  listn_A549[i]<-list(sys_A549_data[,colnames(sys_A549_data) %in% t_sm] %>%
                   .[ rowSums ( is.na ( . ) )  < ncol(.)*0.3,] %>% rownames(.))
}

listn_A549[[1]]
listn_A549[[2]]
overlap_two_vector(listn_A549[[1]],listn_A549[[2]])

##Hela

sys_hela_data<-data_system_0[,37:54]
group_system_hela<-group_system %>% filter(Group=="Hela")
listn_hela<-list()
n<-length(unique(group_system_hela$batch))
for (i in 1:n){
  t_sm<-group_system_hela[group_system_hela$batch==unique(group_system_hela$batch)[i],]$sample
  listn_hela[i]<-list(sys_hela_data[,colnames(sys_hela_data) %in% t_sm] %>%
                        .[ rowSums ( is.na ( . ) )  < ncol(.)*0.3,] %>% rownames(.))
}

overlap_two_vector(listn_hela[[1]],listn_hela[[2]])

##H1975

sys_h1_data<-data_system_0[,19:36]
group_system_h1<-group_system %>% filter(Group=="H1975")
listn_h1<-list()
n<-length(unique(group_system_h1$batch))
for (i in 1:n){
  t_sm<-group_system_h1[group_system_h1$batch==unique(group_system_h1$batch)[i],]$sample
  listn_h1[i]<-list(sys_h1_data[,colnames(sys_h1_data) %in% t_sm] %>%
                        .[ rowSums ( is.na ( . ) )  < ncol(.)*0.5,] %>% rownames(.))
}

overlap_two_vector(listn_h1[[1]],listn_h1[[2]])
##67r

sys_h11_data<-data_system_0[,55:72]
group_system_h11<-group_system %>% filter(Group=="R67")
listn_h11<-list()
n<-length(unique(group_system_h11$batch))
for (i in 1:n){
  t_sm<-group_system_h11[group_system_h11$batch==unique(group_system_h11$batch)[i],]$sample
  listn_h11[i]<-list(sys_h11_data[,colnames(sys_h11_data) %in% t_sm] %>%
                      .[ rowSums ( is.na ( . ) )  < ncol(.)*0.5,] %>% rownames(.))
}

overlap_two_vector(listn_h11[[1]],listn_h11[[2]])



#Intensity cor
##A549
sys_A549_data$mean_chip1<-apply(sys_A549_data[,1:12],1,function(x)mean(x,na.rm=T))
sys_A549_data$mean_chip2<-apply(sys_A549_data[,13:18],1,function(x)mean(x,na.rm=T))
p_sys_1<-scatter_plot_with_r_and_p_3(log2(sys_A549_data[,19:20]))+
  
  scale_color_distiller(palette= 3, direction=1)+
  xlab("Chip1_A549")+
  ylab("Chip2_A549")+
  theme_bw() +
  xlim(9,22)+
  ylim(9,22)+
  theme(  legend.position = "none",
          #axis.line = element_line(size=1.5),
          #axis.text.x = element_text(color="black", size=12, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
          axis.text.x = element_text(color="black", size=12, face="plain") ,
          axis.text.y = element_text(color="black", size=12, face="plain"),
          axis.title.x = element_text(color="black", size=16, face="bold"),
          axis.title.y = element_text(color="black", size=16, face="bold"),
          legend.text=element_text(color="black", size=16, face="bold"),
          legend.title=element_text(color="black", size=16, face="bold"),
          title = element_text(color="black", size=16, face="bold"))

##Hela
sys_hela_data$mean_chip1<-apply(sys_hela_data[,1:12],1,function(x)mean(x,na.rm=T))
sys_hela_data$mean_chip2<-apply(sys_hela_data[,13:18],1,function(x)mean(x,na.rm=T))
p_sys_2<-scatter_plot_with_r_and_p_3(log2(sys_hela_data[,19:20]))+
  
  scale_color_distiller(palette= 3, direction=1)+
  xlab("Chip1_HeLa")+
  ylab("Chip2_HeLa")+
  theme_bw() +
  xlim(9,22)+
  ylim(9,22)+
  theme(  legend.position = "none",
          #axis.line = element_line(size=1.5),
          #axis.text.x = element_text(color="black", size=12, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
          axis.text.x = element_text(color="black", size=12, face="plain") ,
          axis.text.y = element_text(color="black", size=12, face="plain"),
          axis.title.x = element_text(color="black", size=16, face="bold"),
          axis.title.y = element_text(color="black", size=16, face="bold"),
          legend.text=element_text(color="black", size=16, face="bold"),
          legend.title=element_text(color="black", size=16, face="bold"),
          title = element_text(color="black", size=16, face="bold"))

##H1975
sys_h1_data$mean_chip1<-apply(sys_h1_data[,1:12],1,function(x)mean(x,na.rm=T))
sys_h1_data$mean_chip2<-apply(sys_h1_data[,13:18],1,function(x)mean(x,na.rm=T))
p_sys_3<-scatter_plot_with_r_and_p_3(log2(sys_h1_data[,19:20]))+
  
  scale_color_distiller(palette= 3, direction=1)+
  xlab("Chip1_NCI-H1975")+
  ylab("Chip2_NCI-H1975")+
  theme_bw() +
  xlim(9,22)+
  ylim(9,22)+
  theme(  legend.position = "none",
          #axis.line = element_line(size=1.5),
          #axis.text.x = element_text(color="black", size=12, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
          axis.text.x = element_text(color="black", size=12, face="plain") ,
          axis.text.y = element_text(color="black", size=12, face="plain"),
          axis.title.x = element_text(color="black", size=16, face="bold"),
          axis.title.y = element_text(color="black", size=16, face="bold"),
          legend.text=element_text(color="black", size=16, face="bold"),
          legend.title=element_text(color="black", size=16, face="bold"),
          title = element_text(color="black", size=16, face="bold"))

##67R
sys_h11_data$mean_chip1<-apply(sys_h11_data[,1:12],1,function(x)mean(x,na.rm=T))
sys_h11_data$mean_chip2<-apply(sys_h11_data[,13:18],1,function(x)mean(x,na.rm=T))
p_sys_4<-scatter_plot_with_r_and_p_3(log2(sys_h11_data[,19:20]))+
  
  scale_color_distiller(palette= 3, direction=1)+
  xlab("Chip1_67R")+
  ylab("Chip2_67R")+
  theme_bw() +
  xlim(9,22)+
  ylim(9,22)+
  theme(  legend.position = "none",
          #axis.line = element_line(size=1.5),
          text = element_text(size = 16, face="plain"),
          axis.text = element_text(color="black",size=16, face="plain"),
          
          #axis.text.x = element_text(color="black", size=12, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
          axis.text.x = element_text(color="black", size=12, face="plain") ,
          axis.text.y = element_text(color="black", size=12, face="plain"),
          axis.title.x = element_text(color="black", size=16, face="bold"),
          axis.title.y = element_text(color="black", size=16, face="bold"),
          legend.text=element_text(color="black", size=16, face="bold"),
          legend.title=element_text(color="black", size=16, face="bold"),
          title = element_text(color="black", size=16, face="bold"))
library(patchwork)

p_sys_patch<-(p_sys_1+p_sys_2+p_sys_3+p_sys_4)+plot_layout(ncol = 4)
pdf(file="/Volumes/YZC_X6/experiment/program/single_cell/analysis/sc_analysis/analysis/system/20231110_sys_patch_cor.pdf",width = 20,height =5)
p_sys_patch
dev.off()





#systemUMAP
umap_67r<-overview_umap(data = log2(sys_h11_data[,1:18]),data2 = group_system_h11,group = "batch",eclipse_level = 0.95,group_values = c("#4682B4","#8E46B4"),title = "")+theme(legend.position="none")
umap_h1975<-overview_umap(data = log2(sys_h1_data[,1:18]),data2 = group_system_h1,group = "batch",eclipse_level = 0.95,group_values = c("#4682B4","#8E46B4"),title = "")+theme(legend.position="none")
umap_hela<-overview_umap(data = log2(sys_hela_data[,1:18]),data2 = group_system_hela,group = "batch",eclipse_level = 0.95,group_values = c("#4682B4","#8E46B4"),title = "")+theme(legend.position="none")
umap_a549<-overview_umap(data = log2(sys_A549_data[,1:18]),data2 = group_system_a549,group = "batch",eclipse_level = 0.95,group_values = c("#4682B4","#8E46B4"),title = "")+theme(legend.position="none")
p_sys_ymap_patch<-(umap_a549+umap_hela+umap_h1975+umap_67r)+plot_layout(ncol = 4)
pdf(file="analysis/system/20231110_sys_patch_umap.pdf",width = 20,height =5)
p_sys_ymap_patch
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
overview_umap<-function(data,data2,group=c("TNM_stage"),eclipse_level=0.95,group_values=c("#999999", "#E69F00","#56B4E9"),title=c("")){
  library(umap)
  library(ggplot2)
  library(ggrepel)
  data[is.na(data)]<-0.01
  blca_data<-data.frame(t(data))
  
  complete.cases(blca_data)
  #blca_data<- blca_data %>% filter()
  group_all<-data2
  blca.label<-group_all[,group]
  blca.label2<-group_all[,"sample"]
  blca.umap = umap(na.omit(blca_data), n_components = 2, random_state = 15) 
  layout <- blca.umap[["layout"]] 
  layout <- data.frame(layout) 
  final <- cbind(layout, blca.label,blca.label2)
  
    p<-ggplot(final, aes(x=X1, y=X2,colour=blca.label,fill=blca.label)) + 
      geom_point()+
      #labs(title=title)+
      #geom_text_repel(label=blca.label2)+
      xlab("UMAP1")+
      ylab("UMAP2")+
      theme_bw() +
    scale_color_manual(values=group_values)+
      scale_fill_manual(values=group_values)+
    stat_ellipse(level = eclipse_level,geom = "polygon",alpha=0.1)+
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+   #空白背景
    theme(
      
      axis.text.x = element_text(color="black", size=16, face="plain") ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
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

stack_bar_plot_twocluster<-function(df_test){
  #NMF做分型
  #df_test<-metadata_t1_98 %>% filter(Type=="NAC"&pro_rank5_80_nmf=="3")%>%dplyr::select(pro_rank5_80_nmf,Response_res) 
  # 计算每个组合的计数
  library(reshape2)
  library(dplyr)
  cat("谁为x轴，谁就在左边")
  counts <- table(df_test[,2],df_test[,1])
  print(counts)
  # 将计数数据转换为数据框格式
  counts_df <- as.data.frame.matrix(counts) %>% 
    mutate(cluster=rownames(.)) %>% # 仅对数值列进行归一化
    mutate_if(is.numeric, ~ .*100/sum(.)) %>%
    melt(.)
  
  n=length(unique(counts_df$cluster))
  
  # 使用ggplot2创建堆叠柱状图
  p<-ggplot(counts_df, aes(x = as.character(variable), y = value, fill = cluster)) +
    geom_bar(stat = "identity", position = "stack") +
    labs( y = "Percent(%)") +
    scale_fill_manual(values = c("#4682B4","#8E46B4","#B44656","#AEB446","#46B462")[1:n])+
    ggtitle("") +
    #geom_text(label=)
    theme_bw()+
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
  return(p)
}
cv.coef<-function(x){
  sd(x,na.rm = T)/mean(x,na.rm = T)
}
