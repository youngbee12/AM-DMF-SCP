########################
######导入数据##########
########################

naiyao_sample<-filter_data %>% arrange(Group)%>%filter(Group=="H1975"|Group=="R67") %>% .[,c(1,4)]
data_naiyao <-data_order[,naiyao_sample$sample]

res.pca2 = PCA(t(data_naiyao), scale.unit = TRUE, ncp = 8,  graph = T,axes = 1:3)
pca2<-fviz_pca_ind(res.pca2,
                   label = "none",
                   habillage = factor(naiyao_sample$Group),
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

DEP_sc<-wil_cox_nopair_test(data_naiyao,n1=1,n2=18,n3=19,n4=36,Fc_value = 1.5,method = "BH",p_value = 0.05)
vocanol_DEP_sc<-vocanol_plot(DEP_sc$data,p.value = 0.05,Fc_value = 1.5)+
  xlab("Protein mean( Log2 fold-change, 67R vs NCI-H1975)")+
  xlim(-4,5)+theme(legend.position = "bottom",legend.background = element_blank())
pdf(file="./analysis/H1975_vs_R67/20231016_vol_H1975_vs_R67.pdf",height = 9,width = 7 )
vocanol_DEP_sc
dev.off()
up_sc_clus<-clustern_kegg_go_list(DEP_sc$DEG_up)
down_sc_clus<-clustern_kegg_go_list(DEP_sc$DEG_down)
plot1<-plot_after_KEGG_GO(up_sc_clus,down_sc_clus)
########################
########用到函数########
########################
plot_after_KEGG_GO<-function(term_list_up,term_list_down){
  #term_list_up<-up_sc_clus
  #term_list_down<-down_sc_clus
  plotxx<-function(term_df_up,term_df_down){
    
    plot_dat<-term_df_up %>% 
      arrange(p.adjust) %>% 
      slice_head(n = 8) %>%
      mutate(group=rep("up",8))%>%
      bind_rows(term_df_down %>% 
                  arrange(pvalue) %>% 
                  slice_head(n = 8) %>%
                  mutate(group=rep("down",8))%>%
                  mutate(Count=-.[,"Count"])) %>%
      mutate(name=paste(unname(.[,"ID"]),unname(.[,"Description"]), sep = ":"))%>%
      .[!duplicated(.["name"]),] #去除重复的行，仅保留第一个
    
    
    
    p<- ggplot(plot_dat, aes(x = reorder(name,Count), y = Count,fill=-log10(pvalue))) +
      geom_bar(stat = "identity")+ 
      coord_flip()+
      scale_fill_gradient(low = "#999999", high ="#E69F00")+
      ylab("")+
      theme_bw()+
      theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+   #空白背景
      theme(
        
        axis.text.x = element_text(color="black", size=16, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
        axis.text.y = element_text(color="black", size=12, face="plain"),
        axis.title.x = element_text(color="black", size=16, face="bold"),
        axis.title.y = element_text(color="black", size=16, face="bold"),
        legend.text=element_text(color="black", size=16, face="plain"),
        legend.title=element_text(color="black", size=16, face="bold"),
        title = element_text(color="black", size=16, face="bold")
      ) 
    return(p)
  }
  
  #KEGG
  KEGG_up  <-term_list_up$richterm$KEGG@result
  KEGG_down<-term_list_down$richterm$KEGG@result
  kegg_plot<-plotxx(KEGG_up,KEGG_down)
  
  #GOBP 
  GOBP_up<-term_list_up$cut_richterm$cut_BP@result
  GOBP_down<-term_list_down$cut_richterm$cut_BP@result
  GOBP_plot<-plotxx(GOBP_up,GOBP_down)
  
  #GOCC
  GOCC_up<-term_list_up$cut_richterm$cut_CC@result
  GOCC_down<-term_list_down$cut_richterm$cut_CC@result
  GOCC_plot<-plotxx(GOCC_up,GOCC_down)
  
  #GOMF
  GOMF_up<-term_list_up$cut_richterm$cut_MF@result
  GOMF_down<-term_list_down$cut_richterm$cut_MF@result
  GOMF_plot<-plotxx(GOMF_up,GOMF_down)
  
  list_all<-list(kegg_plot,GOBP_plot,GOCC_plot,GOMF_plot)
  names(list_all)<-c("KEGG","GOBP","GOCC","GOMF")
  return(list_all)
}
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
wil_cox_nopair_test<-function(data_log2,n1=1,n2=120,n3=121,n4=240,p_value=0.01,Fc_value=2,method="BH"){
  data<-data_log2
  data_BH <-  p.adjust(apply(data, 1, function(x)
    wilcox.test(x[seq(n1,n2,1)], x[seq(n3,n4,1)],
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
vocanol_plot<-function(data,p.value=0.01,Fc_value=2){
  n<-ncol(data)
  library(ggplot2)
  library(ggrepel)
  data$sign<-"no"
  data[rownames(data[(data$BH_adjust<p.value)&(data$Log2FC>0),]),n+1]<-"up"
  data[rownames(data[(data$BH_adjust<p.value)&(data$Log2FC< 0),]),n+1]<-"down"
  data[rownames(data[(data$BH_adjust<p.value)&(data$Log2FC>log2(Fc_value)),]),n+1]<-paste0(Fc_value,"x up")
  data[rownames(data[(data$BH_adjust<p.value)&(data$Log2FC< -log2(Fc_value)),]),n+1]<-paste0(Fc_value,"x down")
  cat(table(data$sign))
  data_text<-data[data$sign!= "no",]
  data_text$rowname<-rownames(data_text)
  p<-ggplot(data,aes(x=Log2FC,y=-log10(BH_adjust),color=sign))+
    geom_point()+
    geom_text_repel(data=data_text,aes(label=rowname),size=4)+
    #ggtitle("DiconsensusMIBCerentially expressed genes(FC>1.5;P<0.05)")+
    #xlab("protein median( Log2 fold-change, nonResponse vs Response)")+
    ylab("-Log10 ( adjust p value)")+
    scale_x_continuous()+
    #ylim(0,3)+
    #xlim(-6,7)+
    geom_hline(aes(yintercept=-log10(p.value)),linetype="dashed")+
    geom_vline(aes(xintercept=log2(Fc_value)),linetype="dashed")+
    geom_vline(aes(xintercept=-log2(Fc_value)),linetype="dashed")+
    scale_color_manual(values=c("royalblue4","firebrick4","royalblue1","grey70","firebrick1"))+
    theme_bw() +
    
    theme(  #legend.position = "none",
            axis.line = element_line(size=1.5),
            axis.text.x = element_text(color="black", size=12, face="plain",angle=0,hjust=0.5) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
            axis.text.y = element_text(color="black", size=12, face="plain"),
            axis.title.x = element_text(color="black", size=16, face="bold"),
            axis.title.y = element_text(color="black", size=16, face="bold"),
            legend.text=element_text(color="black", size=16, face="bold"),
            legend.title=element_text(color="black", size=16, face="bold"),
            title = element_text(color="black", size=16, face="bold"))
  return(p)
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

