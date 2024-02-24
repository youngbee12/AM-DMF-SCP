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

A549_speci<-c("AKR1B1","BSG","PHGDH","SLC16A3","BLVRB","G6PD")
Hela_speci<-c("NSUN2","SET","CBX3","ANXA4","UGDH","B2M","SUB1","AKR1C3","TARS1",
              "PHGDH","ALDH2","SYNGR2","SRI","COX7A2","DDX11","OXSR1","PTMS",
              "HLA.A","H1.5","PDLIM1")
HepG2_speci<-c("PLEC","HYOU1","SIPA1","TAGLN2","S100A11",
               "ANXA1","ANXA2","EML4","TXNRD1","H1.0","KYNU")

library(reshape2)
library(ggplot2)
library(broom)
library(tidyverse)


p1_A549<-multi_group_violinplot(data=data_hela_a549,group=hela_a549_sample,selectgenes = A549_speci)
p2_A549<-depmap_lineplot(Depmap_data=select_depmap_Data,selectgenes = A549_speci)
combined_plot_A549 <-  p2_A549 / p1_A549+ plot_layout(ncol = 1, heights = c(1, 2))
pdf(file="./analysis/hela_a549_Hepg2/20231012_hela_a549_violin_A549.pdf",width = 7.5,height = 7.5)
combined_plot_A549
dev.off()
p1_hela<-multi_group_violinplot(data=data_hela_a549,group=hela_a549_sample,selectgenes = Hela_speci)
p2_hela<-depmap_lineplot(Depmap_data=select_depmap_Data,selectgenes = Hela_speci)

combined_plot_hela <-  p2_hela/ p1_hela+ plot_layout(ncol = 1, heights = c(1, 2))
pdf(file="./analysis/hela_a549_Hepg2/20231012_hela_a549_violin_hela.pdf",width = 21,height = 7.5)

combined_plot_hela
dev.off()



p1_HepG2<-multi_group_violinplot(data=data_hela_a549,group=hela_a549_sample,selectgenes = HepG2_speci)
p2_HepG2<-depmap_lineplot(Depmap_data=select_depmap_Data,selectgenes = HepG2_speci)
combined_plot_HepG2 <-  p2_HepG2/ p1_HepG2+ plot_layout(ncol = 1, heights = c(1, 2))
pdf(file="./analysis/hela_a549_Hepg2/20231012_hela_a549_violin_HepG2.pdf",width = 13.5,height = 7.5)
combined_plot_HepG2
dev.off()



#比较两个趋势的相关性
result <- cor.test(x, y, method = "kendall")
aaa<-depmap_lineplot_data(Depmap_data=select_depmap_Data,selectgenes = HepG2_speci)
bbb<-multi_group_violinplot_data(data=data_hela_a549,group=hela_a549_sample,selectgenes = HepG2_speci)


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

multi_group_violinplot_data<-function(data,group,selectgenes=c("TP53","ALDOA")){
  library(reshape2)
  library(ggplot2)
  library(broom)
  library(tidyverse)
  test_data<-data[selectgenes,]
  # 使用split和sapply来计算均值
  grouped_columns <- split(group$sample, group$Group)
  
  result <- sapply(grouped_columns, function(g) {
    rowMeans(test_data[, g, drop = FALSE],na.rm = T)
  })
  
  # 将结果转为数据框并设置行名
  result <- as.data.frame(result)
  rownames(result) <- rownames(test_data)
  
  # 打印结果
 
  return(result)
}
multi_group_violinplot<-function(data,group,selectgenes=c("TP53","ALDOA")){
  library(reshape2)
  library(ggplot2)
  library(broom)
  library(tidyverse)
  test_data<-data[selectgenes,]
  test_data$pro<-rownames(test_data)
  voli<-melt(test_data)
  colnames(voli) <- c("gene", "sample", "expression")
  voli <- voli %>% left_join(group)
  voli$gene<-factor(voli$gene,levels = selectgenes)
  voli$sample<-as.character(voli$sample)
  
  
  # 使用 tidyverse 和 broom 包计算 p 值
  library(broom)
  p_values <- voli %>%
    group_by(gene)%>%
    do(tidy(anova(lm(.$expression ~ .$Group))))%>%
    mutate(significance = case_when(
      p.value < 0.0001 ~ "****",
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    ))
  
  
  # 将 p 值添加到原始数据中
  voli <- inner_join(voli, p_values, by = "gene")
  
 
  p<-voli %>% ggplot(aes(x=gene, y=expression,color=Group)) + 
    geom_violin(scale="width", adjust=0.5) +
    #geom_point(position = position_jitter(width = 0.2), size = 2)+
    geom_point(position = position_jitterdodge(jitter.width = 0.1,seed = 666,dodge.width = 1)) +
    ggtitle("")+
    ylab("Log2(Intensity)")+
    scale_color_manual(values=c("#4682B4" ,"#AF46B4","#B47846","black"))+
    scale_fill_manual(values=c("#4682B4" ,"#AF46B4","#B47846","black"))+
    geom_text(aes(label = significance, y = 21), vjust = 0) +
    #coord_flip()+
    facet_wrap(.~gene, scales = "free_x",nrow=1)+
    #geom_jitter(width=0.3,alpha=0.2)+
    theme_bw()+
    theme(  legend.position = "none",
      #axis.line = element_line(size=1.5),
      text = element_text(size = 16, face="plain"),
      #axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),# Values for face are one of "plain", "italic", "bold" and "bold.italic"
      axis.text.x=element_blank(),
      axis.text.y = element_text(color="black", size=12, face="plain"),
      #axis.title.x = element_text(color="black", size=16, face="bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(color="black", size=16, face="bold"),
      legend.text=element_text(color="black", size=16, face="bold"),
      legend.title=element_text(color="black", size=16, face="bold"),
      #title = element_text(color="black", size=16, face="bold"),
      strip.text = element_blank(),
      strip.background = element_blank()) 
  return(p)
}
depmap_lineplot_data<-function(Depmap_data,selectgenes=c("TP53","ALDOA")){
  library(reshape2)
  library(ggplot2)
  library(broom)
  library(tidyverse)
  
  test_data2<-Depmap_data[selectgenes,]
  return(test_data2)
}
depmap_lineplot<-function(Depmap_data,selectgenes=c("TP53","ALDOA")){
  library(reshape2)
  library(ggplot2)
  library(broom)
  library(tidyverse)

  test_data2<-Depmap_data[selectgenes,]
  test_data2$pro<-selectgenes
  line<-melt(test_data2)
  colnames(line) <- c("gene", "sample", "expression")
  line$gene<-factor(line$gene,levels = selectgenes)
  
  lineplot<-line %>%
    ggplot(aes(x = sample, y = expression, color = sample)) +
    geom_line(aes(x=sample,y = expression,group = gene), color = "grey", size = 0.5) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, seed = 666, dodge.width = 1)) +
    scale_color_manual(values=c("#4682B4" ,"#AF46B4","#B47846","black"))+
    scale_fill_manual(values=c("#4682B4" ,"#AF46B4","#B47846","black"))+
    ylab("Log2(TPM+1)")+
    facet_wrap(~ gene, scales = "free_x", nrow = 1) +
    theme_bw() +
    theme(
      legend.position = "none",
      #axis.line = element_line(size=1.5),
          text = element_text(size = 16, face="plain"),
          #axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),# Values for face are one of "plain", "italic", "bold" and "bold.italic"
          axis.text.y = element_text(color="black", size=12, face="plain"),
          #axis.title.x = element_text(color="black", size=16, face="bold"),
          axis.title.y = element_text(color="black", size=16, face="bold"),
          legend.text=element_text(color="black", size=16, face="bold"),
          legend.title=element_text(color="black", size=16, face="bold"),
          #title = element_text(color="black", size=16, face="bold"),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          #strip.text = element_blank(),
          #strip.background = element_blank()
    )
  return(lineplot)
}
Kendall_tau_cor_2_df<-function(df1,df2){
  df1<-aaa
  df2<-bbb
  a<-overlap_two_vector(rownames(df1),rownames(df2))$common
  df1<-df1[a,]
  df2<-df2[a,]
  cor.test( as.numeric(df1[5,]),as.numeric(df2[5,]),method = "spearman")
  kendall_tau_results <- apply(1:nrow(df1), 1, function(i) {
    cor.test(as.numeric(df1[i,]), as.numeric(df2[i,]), method = "kendall")$estimate
  })
  
  # 结果转为数据框并设置行名
  result <- data.frame(Kendall_Tau = kendall_tau_results)
  rownames(result) <- rownames(df1)
  return(result)
}
Kendall_tau_cor_2_df(aaa,bbb)
