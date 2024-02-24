data_pep_select
hela_hand_sample
data_pep_50<-filter_na_bygroup(data_pep_select,group_info = dplyr::select(hela_hand_sample,c("sample","Group")),x=0.5)
#Knn补缺
library(impute)
half_data_meadian_impute_pep<-impute.knn(as.matrix(log2(data_pep_50)) ,k = 3, rowmax = 0.7, colmax = 0.8, maxp = 1500, rng.seed=362436069)
half_data_meadian_impute_pep<-data.frame(half_data_meadian_impute_pep$data)
colnames(half_data_meadian_impute_pep)<-colnames(data_pep_50)
data_hela_hand_pep<-half_data_meadian_impute_pep  

data_hela_hand_pep

DEP_pep<-wil_cox_nopair_test(data_hela_hand_pep,n1=7,n2=12,n3=1,n4=6,Fc_value = 1.2,method = "none",p_value = 0.05)
outchip_up_G<-unlist(lapply(DEP_pep$DEG_up,function(x)pep_GRVAY_score_calc(input = x)))
inchip_up_G<-unlist(lapply(DEP_pep$DEG_down,function(x)pep_GRVAY_score_calc(input = x)))

outchip_up_misc<-unlist(lapply(DEP_pep$DEG_up,function(x)pep_GRVAY_score_calc(input = x)))
inchip_up_misc<-unlist(lapply(DEP_pep$DEG_down,function(x)pep_GRVAY_score_calc(input = x)))

t.test(outchip_up_G,inchip_up_G)
plot(xx)
listn<-list()
#data=data_total
n<-length(unique(hela_hand_sample$Group))
for (i in 1:n){
  t_sm<-hela_hand_sample[hela_hand_sample$Group==unique(hela_hand_sample$Group)[i],]$sample
  listn[i]<-list(data_pep_select[,colnames(data_pep_select) %in% t_sm] %>%
                   .[ rowSums ( is.na ( . ) )  < ncol(.)*0.5,] %>% rownames(.))
}
overlap_two_vector(listn[[1]],listn[[2]])
com<-overlap_two_vector(listn[[1]],listn[[2]])$common
up_out<-overlap_two_vector(listn[[1]],listn[[2]])$uni_y
up_in<-overlap_two_vector(listn[[1]],listn[[2]])$uni_x
outchip_up_G<-unlist(lapply(c(DEP_pep$DEG_up,up_out),function(x)pep_GRVAY_score_calc(input = x)))
inchip_up_G<-unlist(lapply(c(DEP_pep$DEG_down,up_in),function(x)pep_GRVAY_score_calc(input = x)))
Grvay_df<-data.frame(Group=c(rep("inchip",3404),rep("outchip",707)),GRVAY_score=c(inchip_up_G,outchip_up_G))
pdf(file="analysis/hela_hand/20231020_GRVAY_hela_hand.pdf",width = 5,height = 5)
Grvay_df %>% ggplot(aes(x=Group, y=GRVAY_score,color=Group)) + 
  geom_boxplot(width=0.5) +
  #geom_point(position = position_jitterdodge()) +
  ggtitle("GRVAY_score_out_in")+
  ylab("GRVAY_score")+
  scale_color_manual(values=c("#4682B4" ,"#AF46B4"))+
  stat_compare_means(method = "t.test")+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "inchip")+                  # Pairwise comparison against all
  #geom_jitter(color="black", size=0.1, alpha=0.2)+
  theme_bw() +
  
  theme(  legend.position = "none",
          #axis.line = element_line(size=1.5),
          #axis.text.x = element_text(color="black", size=12, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
          axis.text.x = element_text(color="black", size=16, face="plain") ,
          axis.text.y = element_text(color="black", size=16, face="plain"),
          axis.title.x = element_text(color="black", size=16, face="bold"),
          axis.title.y = element_text(color="black", size=16, face="bold"),
          legend.text=element_text(color="black", size=16, face="bold"),
          legend.title=element_text(color="black", size=16, face="bold"),
          title = element_text(color="black", size=16, face="bold"))
dev.off()

test_pep<-DEP_pep$data
test_pep$G_score<-unlist(lapply(rownames(test_pep),function(x)pep_GRVAY_score_calc(input = x)))
pdf("./analysis/hela_hand/20231011_corr_hela_hand.pdf",height = 5,width = 5)
scatter_plot_with_r_and_p_3(df1)+
  
  scale_color_distiller(palette= 3, direction=1)+
  xlab("Inchip")+
  ylab("Outchip")+
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
dev.off()
df1<-data.frame(Outchip=apply(data_hela_hand[,1:6],1,mean),Inchip=apply(data_hela_hand[,7:12],1,mean))

scatter_plot_with_r_and_p <- function(df) {
  library(ggpointdensity)
  # 确保数据框包含两列
  if (ncol(df) != 2) {
    stop("Data frame must have exactly two columns.")
  }
  colnames(df)<-c("num_1","num_2")
  # 计算相关系数、 R² 值和 p 值
  correlation <- cor.test(df[[1]], df[[2]])
  r_squared <- correlation$estimate
  p_value <- correlation$p.value
  
  # 创建散点图并添加 R² 值和 p 值
  p<-ggplot(df, aes(x = num_1, y = num_2))+
    geom_point(size=0.5) +
    geom_pointdensity() +
   # scale_color_viridis(discrete = F)+
    #scale_color_distiller(palette= 3, direction=1) +
    geom_smooth(method = "lm", se = FALSE, col = "red", linetype = "dashed",linewidth=0.5) +
    labs(title = paste("R = ", round(r_squared, digits = 2), ", p = ", formatC(p_value, format = "e", digits = 3)))+
    
    #xlim(0,1)+
    #ylim(0,1)+
    theme_bw() +
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background =      element_blank(),axis.line =  element_line(colour = "black"))+   #空白背景
    theme(
      #legend.position  ="none",
      axis.text.x = element_text(color="black",size=16, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
      axis.text.y = element_text(color="black", size=16, face="plain"),
      axis.title.x = element_text(color="black", size=24, face="bold"),
      axis.title.y = element_text(color="black", size=24, face="bold"))+
    theme(text = element_text(family = "Helvetica"))
  return(p)
}
scatter_plot_with_r_and_p_1 <- function(df) {
  # 确保数据框包含两列
  if (ncol(df) != 2) {
    stop("Data frame must have exactly two columns.")
  }
  colnames(df)<-c("num_1","num_2")
  # 计算相关系数、 R² 值和 p 值
  correlation <- cor.test(df[[1]], df[[2]])
  r_squared <- correlation$estimate
  p_value <- correlation$p.value
  
  # 创建散点图并添加 R² 值和 p 值
  ggplot(df, aes(x=num_1, y=num_2) ) +
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
    scale_fill_distiller(palette= 4, direction=1) +
    geom_smooth(method = "lm", se = FALSE, col = "red", linetype = "dashed") +
    #scale_x_continuous(expand = c(0, 0)) +
    #scale_y_continuous(expand = c(0, 0)) +
    theme_bw() +
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background =      element_blank(),axis.line =  element_line(colour = "black"))+   #空白背景
    theme(
      #legend.position  ="none",
      axis.text.x = element_text(color="black",size=16, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
      axis.text.y = element_text(color="black", size=16, face="plain"),
      axis.title.x = element_text(color="black", size=24, face="bold"),
      axis.title.y = element_text(color="black", size=24, face="bold"))+
    theme(text = element_text(family = "Helvetica"))
  return(p)
}
df<-log2(sys_A549_data[,19:20])
scatter_plot_with_r_and_p_3 <- function(df) {
  library(ggpointdensity)
  # 确保数据框包含两列
  if (ncol(df) != 2) {
    stop("Data frame must have exactly two columns.")
  }
  colnames(df)<-c("num_1","num_2")
  # 计算相关系数、 R² 值和 p 值
  correlation <- cor.test(df[[1]], df[[2]])
  r <- correlation$estimate
  pvalue <- correlation$p.value
  pvalue<-paste0("r = ", round(r, 2), ", ",ifelse(pvalue<0.0001,"P<0.0001",pvalue))
  # 创建散点图并添加 R² 值和 p 值
  p<-ggplot(df, aes(x = num_1, y = num_2))+
    geom_point(size=0.5) +
    geom_pointdensity() +
    # scale_color_viridis(discrete = F)+
    #scale_color_distiller(palette= 3, direction=1) +
    #geom_smooth(method = "lm", se = FALSE, col = "red", linetype = "dashed",linewidth=0.5) +
    geom_text(aes(x = 10, y = 20, label =pvalue ),
              hjust = 0, vjust = 1,size=8,color="black") +  
    
   # labs(title = paste("R = ", round(r_squared, digits = 2), ", p = ", formatC(p_value, format = "e", digits = 3)))+
    #xlim(0,1)+
    #ylim(0,1)+
    theme_bw() +
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background =      element_blank(),axis.line =  element_line(colour = "black"))+   #空白背景
    theme(
      #legend.position  ="none",
      text = element_text(color="black",size = 16, face="plain"),
      axis.text = element_text(color="black",size=16, face="plain"),
      axis.text.x = element_text(color="black",size=16, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
      axis.text.y = element_text(color="black", size=16, face="plain"),
      axis.title.x = element_text(color="black", size=24, face="bold"),
      axis.title.y = element_text(color="black", size=24, face="bold"))
    
  return(p)
}
