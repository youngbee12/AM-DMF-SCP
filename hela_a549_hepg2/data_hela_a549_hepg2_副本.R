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

hela_a549_sample<-sc_group %>% arrange(Group)%>%filter(Group=="A549"|Group=="Hela"|Group=="HepG2") %>% .[,c(1,3)]
data_hela_a549 <-sc_data[,hela_a549_sample$sample]



load(file="./data/20231011_scdata_pre.Rdata")
sc_data_raw<-sc_data






hela_a549_sample<-sc_group %>% arrange(Group)%>%filter(Group=="A549"|Group=="Hela"|Group=="HepG2") %>% .[,c(1,3)]
data_hela_a549 <-sc_data[,hela_a549_sample$sample]


#Intensity boxplot
#强度
data_hela_a549_int<-data_hela_a549
data_hela_a549_int$Genename<-rownames(data_hela_a549_int)
newprotein_after<-tidyr::gather(data_hela_a549_int,key="group",value="intensity",-Genename)#宽数据转长数据
rownames(newprotein_after)<-NULL
rownames(newprotein_after)<-newprotein_after$Genename

ggplot(newprotein_after,aes(x=group,y=as.numeric(intensity)))+
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







listn<-list()
#data=data_total
#n<-length(unique(hela_hand_sample$Group))
#for (i in 1:n){
# t_sm<-hela_hand_sample[hela_hand_sample$Group==unique(hela_hand_sample$Group)[i],]$sample
# listn[i]<-list(sc_data_raw_hela_hand[,colnames(sc_data_raw_hela_hand) %in% t_sm] %>%
#                  .[ rowSums ( is.na ( . ) )  < ncol(.)*0.5,] %>% rownames(.))
#}
#overlap_two_vector(listn[[1]],listn[[2]])
data_50<-filter_na_bygroup(data_hela_a549,group_info = dplyr::select(hela_a549_sample,c("sample","Group")),x=0.5)
#Knn补缺
library(impute)
half_data_meadian_impute<-impute.knn(as.matrix(log2(data_50)) ,k = 3, rowmax = 0.7, colmax = 0.8, maxp = 1500, rng.seed=362436069)
half_data_meadian_impute<-data.frame(half_data_meadian_impute$data)
colnames(half_data_meadian_impute)<-colnames(data_50)
data_hela_a549<-half_data_meadian_impute

pdf(file="./analysis/hela_a549_Hepg2/20231016_umap_hela_a549.pdf",height = 5,width = 6 )
overview_umap(data = data_hela_a549,data2 =hela_a549_sample,group = c("Group"),eclipse_level = 0.95,group_values = wheel("steelblue", num = 4),n_neighbors = 15)
dev.off()
res.pca4 = PCA(t(data_hela_a549), scale.unit = TRUE, ncp = 8,  graph = T,axes = 1:3)
pca4<-fviz_pca_ind(res.pca4,
                   label = "none",
                   habillage = factor(hela_a549_sample$Group),
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
wide_df_hela_a549<-wide_df %>% filter(Group=="A549"|Group=="Hela"|Group=="HepG2")
filter_data_hela_a549<-sc_group %>%filter(Group=="A549"|Group=="Hela"|Group=="HepG2")
pdf(file="./analysis/hela_a549_Hepg2/20231009_identi_hela_a549.pdf",height = 5,width = 5 )
ggplot(wide_df_hela_a549, aes(x=Group, y=mean,color=Group,fill=Group)) +
  geom_bar(position=position_dodge(), stat="identity",colour='white',width=.6,alpha=0.8) +
  # scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.4)+
  scale_color_manual(values=c("#4682B4","#AF46B4","#B47846"))+
  scale_fill_manual(values=c("#4682B4","#AF46B4","#B47846"))+
  geom_jitter(data = filter_data_hela_a549,aes(x=Group,y=p_num),width=.1)+
  scale_y_continuous(limits=c(0,3100),breaks=seq(0,3000,500))+
  geom_text(aes(label = ceiling(mean)), vjust = -9, colour = "black")+
  ylab("Proteins")+
  xlab("")+
  theme_bw() +
  
  theme(  legend.position = "none",
          axis.line = element_line(size=1.5),
          axis.text.x = element_text(color="black", size=12, face="plain",angle=0,hjust=0.5) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
          axis.text.y = element_text(color="black", size=12, face="plain"),
          axis.title.x = element_text(color="black", size=16, face="bold"),
          axis.title.y = element_text(color="black", size=16, face="bold"),
          legend.text=element_text(color="black", size=16, face="bold"),
          legend.title=element_text(color="black", size=16, face="bold"),
          title = element_text(color="black", size=16, face="bold"))

dev.off()

pdf(file = "./analysis/hela_a549_Hepg2/20231024_corr_hela_a549.pdf",height = 5,width = 5)

library(corrplot)
corrplot(cor(log2(sc_data[,hela_a549_sample$sample]),method="pearson",use="complete.obs"),
         method="color",tl.col = "black",
         tl.pos="n",
         col=c(rep("white", 500),rep('grey', 100) ,COL2('RdBu', 400)),
         col.lim=c(0.5,1),
         #addCoef.col = "white",
         #type="upper")
)
dev.off()

group<-factor(paste0("group",arrange(filter_data_hela_a549,Group)$Group))
design<-model.matrix(~0+group)
colnames(design)<-levels(group)
library(limma)
#创建对比矩阵
group_names <-levels(group)
contrasts_names<-combn(group_names,2,paste,collapse="-")
contrast_formulas<-combn(group_names,2,function(x)(paste0(x[1],"-",x[2])))
contrast_matrix<-makeContrasts(contrasts = setNames(contrast_formulas,contrasts_names),levels = design)
fit <-lmFit(data_hela_a549,design)
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
#data_heat<-data_hela_a549[c(aa[[1]]$DEG_up,aa[[2]]$DEG_up,aa[[3]]$DEG_up),]
data_heat<-data_hela_a549[c(aa[[1]]$DEG_up,
                            aa[[2]]$DEG_up,aa[[1]]$DEG_down,aa[[2]]$DEG_down,
                            aa[[3]]$DEG_up,aa[[3]]$DEG_down),]
pdf(file = "./analysis/hela_a549_Hepg2/20231012_HEATMAP_hela_a549.pdf",height = 5,width = 6)
A549_up_clus<-clustern_kegg_go_list(aa[[1]]$DEG_up)
A549_up_clus_1<-rbind(A549_up_clus$richterm$KEGG@result,
          A549_up_clus$cut_richterm$cut_BP@result,
          A549_up_clus$cut_richterm$cut_CC@result,
          A549_up_clus$cut_richterm$cut_MF@result)%>% arrange(p.adjust,Count)
Hela_up_clus_1<-rbind(Hela_up_clus$richterm$KEGG@result,
                      Hela_up_clus$cut_richterm$cut_BP@result,
                      Hela_up_clus$cut_richterm$cut_CC@result,
                      Hela_up_clus$cut_richterm$cut_MF@result)%>% arrange(p.adjust,Count)
A549_down_clus_1<-rbind(A549_down_clus$richterm$KEGG@result,
                        A549_down_clus$cut_richterm$cut_BP@result,
                      A549_down_clus$cut_richterm$cut_CC@result,
                      A549_down_clus$cut_richterm$cut_MF@result)%>% arrange(p.adjust,Count)
HepG2_down_clus_1<-rbind(HepG2_down_clus$richterm$KEGG@result,
                         HepG2_down_clus$cut_richterm$cut_BP@result,
                        HepG2_down_clus$cut_richterm$cut_CC@result,
                        HepG2_down_clus$cut_richterm$cut_MF@result)%>% arrange(p.adjust,Count)
HepG2_up_clus_1<-rbind(HepG2_up_clus$richterm$KEGG@result,
                       HepG2_up_clus$cut_richterm$cut_BP@result,
                       HepG2_up_clus$cut_richterm$cut_CC@result,
                       HepG2_up_clus$cut_richterm$cut_MF@result) %>% arrange(p.adjust,Count)
p1<-ggplot(A549_up_clus_1[c(2,3),], aes(x =reorder(Description,Count), y = Count,fill=-log10(p.adjust))) +
  geom_bar(stat = "identity")+ 
  coord_flip()+
  scale_y_reverse() +
  #scale_x_discrete(position = "top")+
  scale_fill_gradient(low = "blue", high ="black", limits = c(0, 15))+
  ylab("")+
  xlab("")+
  geom_text(aes(label = ceiling(Count)), hjust =-2, colour = "white",size=8)+
  #geom_text(aes(label = Description),  colour = "white",nudge_x = 0,nudge_y = 2.5)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+   #空白背景
  theme(
    #legend.position = "none",
    axis.text.x=element_blank(),
    axis.ticks = element_blank(),
    #axis.text.x = element_text(color="black", size=16, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
    axis.text.y = element_text(color="black", size=12, face="plain"),
    #axis.text.y = element_blank(),
    axis.title.x = element_text(color="black", size=16, face="bold"),
    axis.title.y = element_text(color="black", size=16, face="bold"),
    legend.text=element_text(color="black", size=16, face="plain"),
    legend.title=element_text(color="black", size=16, face="bold"),
    title = element_text(color="black", size=16, face="bold")
  ) 
p2<-ggplot(Hela_up_clus_1[c(1,3,11),], aes(x = reorder(Description,Count), y = Count,fill=-log10(p.adjust))) +
  geom_bar(stat = "identity")+ 
  coord_flip()+
  scale_y_reverse() +
  #scale_x_discrete(position = "top")+
  scale_fill_gradient(low = "blue", high ="black", limits = c(0, 15))+
  ylab("")+
  xlab("")+
  geom_text(aes(label = ceiling(Count)), hjust =-2, colour = "white",size=8)+
  #geom_text(aes(label = Description),  colour = "white",nudge_x = 0,nudge_y = 2.5)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+   #空白背景
  theme(
    #legend.position = "none",
    axis.text.x=element_blank(),
    axis.ticks = element_blank(),
    #axis.text.x = element_text(color="black", size=16, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
    axis.text.y = element_text(color="black", size=12, face="plain"),
    #axis.text.y = element_blank(),
    axis.title.x = element_text(color="black", size=16, face="bold"),
    axis.title.y = element_text(color="black", size=16, face="bold"),
    legend.text=element_text(color="black", size=16, face="plain"),
    legend.title=element_text(color="black", size=16, face="bold"),
    title = element_text(color="black", size=16, face="bold")
  ) 
p3<-ggplot(A549_down_clus_1[c(1),], aes(x = reorder(Description,Count), y = Count,fill=-log10(p.adjust))) +
  geom_bar(stat = "identity")+ 
  coord_flip()+
  scale_y_reverse() +
  #scale_x_discrete(position = "top")+
  scale_fill_gradient(low = "blue", high ="black", limits = c(0, 15))+
  ylab("")+
  xlab("")+
  geom_text(aes(label = ceiling(Count)), hjust =-2, colour = "white",size=8)+
  #geom_text(aes(label = Description),  colour = "white",nudge_x = 0,nudge_y = 2.5)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+   #空白背景
  theme(
    #legend.position = "none",
    axis.text.x=element_blank(),
    axis.ticks = element_blank(),
    #axis.text.x = element_text(color="black", size=16, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
    axis.text.y = element_text(color="black", size=12, face="plain"),
    #axis.text.y = element_blank(),
    axis.title.x = element_text(color="black", size=16, face="bold"),
    axis.title.y = element_text(color="black", size=16, face="bold"),
    legend.text=element_text(color="black", size=16, face="plain"),
    legend.title=element_text(color="black", size=16, face="bold"),
    title = element_text(color="black", size=16, face="bold")
  ) 
p4<-ggplot(HepG2_up_clus_1[c(1,3,19,22),], aes(x = reorder(Description,Count), y = Count,fill=-log10(p.adjust))) +
  geom_bar(stat = "identity")+ 
  coord_flip()+
  scale_y_reverse() +
 # scale_x_discrete(position = "top")+
  scale_fill_gradient(low = "blue", high ="black", limits = c(0, 15))+
  ylab("")+
  xlab("")+
  geom_text(aes(label = ceiling(Count)), hjust =-2, colour = "white",size=8)+
  #geom_text(aes(label = Description),  colour = "white",nudge_x = 0,nudge_y = 2.5)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+   #空白背景
  theme(
    #legend.position = "none",
    axis.text.x=element_blank(),
    axis.ticks = element_blank(),
    #axis.text.x = element_text(color="black", size=16, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
    axis.text.y = element_text(color="black", size=12, face="plain"),
    #axis.text.y = element_blank(),
    axis.title.x = element_text(color="black", size=16, face="bold"),
    axis.title.y = element_text(color="black", size=16, face="bold"),
    legend.text=element_text(color="black", size=16, face="plain"),
    legend.title=element_text(color="black", size=16, face="bold"),
    title = element_text(color="black", size=16, face="bold")
  ) 
p5<-ggplot(HepG2_down_clus_1[c(2,6,17,20),], aes(x = reorder(Description,Count), y = Count,fill=-log10(p.adjust))) +
  geom_bar(stat = "identity")+ 
  coord_flip()+
  scale_y_reverse() +
  #scale_x_discrete(position = "top")+
  scale_fill_gradient(low = "blue", high ="black", limits = c(0, 15))+
  ylab("")+
  xlab("")+
  geom_text(aes(label = ceiling(Count)), hjust =-1, colour = "white",size=8)+
  #geom_text(aes(label = Description),  colour = "white",nudge_x = 0,nudge_y = 2.5)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+   #空白背景
  theme(
    #legend.position = "none",
    axis.text.x=element_blank(),
    axis.ticks = element_blank(),
    #axis.text.x = element_text(color="black", size=16, face="plain",angle=45,hjust=1) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
    axis.text.y = element_text(color="black", size=12, face="plain"),
    #axis.text.y = element_blank(),
    axis.title.x = element_text(color="black", size=16, face="bold"),
    axis.title.y = element_text(color="black", size=16, face="bold"),
    legend.text=element_text(color="black", size=16, face="plain"),
    legend.title=element_text(color="black", size=16, face="bold"),
    title = element_text(color="black", size=16, face="bold")
  ) 
library(patchwork)
right_panel<-p1/p2/p3/p4/p5+plot_layout(ncol = 1, heights = c(0.3, 0.5,0.2,1,1))



# Grab the HEL plot (assuming you've already defined a function to draw it)
left_grid <- grid.grabExpr(draw(HEL))

# Combine the left grid and right plot
final_grid <- plot_grid(right_panel, left_grid, ncol = 2, rel_widths = c(1, 1.5), labels = c("A", "B"))
pdf(file = "./analysis/hela_a549_Hepg2/20231024_HEATMAP_hela_a549_cob.pdf",height = 8,width = 20)

final_grid
dev.off()
plot_grid(
 
  plot_grid(plotlist = list(p1,p2,p3,p4,p5), nrow = 5),
  grid.grabExpr(draw(HEL)),
  ncol = 2,
  rel_widths = c(1, 2),
 rel_heights = c(c(4,2,1,1,1),1),
  labels = c("A", "B", "C")  # 可选：添加图形标签
)
class(A549_up_clus$cut_richterm$cut_BP@result)

Hela_up_clus<-clustern_kegg_go_list(aa[[2]]$DEG_up)

HepG2_up_clus<-clustern_kegg_go_list(aa[[3]]$DEG_up,aa[[2]]$DEG_down)

A549_down_clus<-clustern_kegg_go_list(aa[[1]]$DEG_down)

HepG2_down_clus<-clustern_kegg_go_list(aa[[3]]$DEG_down)
anno_df = data.frame(
  
  mod = c(rep("A549_up",length(c(aa[[1]]$DEG_up))),
              rep("Hela_up",length(c(aa[[2]]$DEG_up))),
          rep("A549_down",length(c(aa[[1]]$DEG_down))),
          rep("HepG2_up",length(c(aa[[2]]$DEG_down,aa[[3]]$DEG_up))),
              rep("HepG2_down",length(c(aa[[3]]$DEG_down)))
))
ha = rowAnnotation(df = anno_df,
                       col = list(Group=c("A549_up"=wheel("steelblue", num = 5)[1],
                                          "Hela_up"=wheel("steelblue", num = 5)[2],
                                          "A549_down"=wheel("steelblue", num = 5)[3],
                                          "HepG2_up"=wheel("steelblue", num = 5)[4],
                                          "HepG2_down"=wheel("steelblue", num = 5)[5])))
HEL<-Heatmap(
  t(scale_neg1_1(t(data_heat))), 
  name = "mat",
  #col = col_fun,
  cluster_columns = F,
  cluster_rows=F,
  show_row_names = F,
  show_column_names = F,
  show_row_dend =F,
  show_column_dend = F,
  #row_km = 5,
  column_split = arrange(hela_a549_sample,Group)$Group,
  row_names_gp = gpar(fontsize = 10),
#right_annotation = ha,
top_annotation = HeatmapAnnotation(df=hela_a549_sample %>% arrange(Group)%>%
                                     dplyr::select(Group),
                                   col = list(Group=c("A549"=wheel("steelblue", num = 4)[1],
                                                      "Hela"=wheel("steelblue", num = 4)[2],
                                                      
                                                      "HepG2"=wheel("steelblue", num = 4)[3]))
)
)
dev.off()

write.table(protein_data,file="20231016_sc_Data.txt",sep = "\t")

library(cowplot)

library(ggplot2)

p1 <- ggplot(mtcars, aes(x=mpg, y=disp)) + geom_point() + ggtitle("MPG vs Displacement")
p2 <- ggplot(mtcars, aes(x=mpg, y=hp)) + geom_point() + ggtitle("MPG vs Horsepower")

plot_grid(p1, p2, labels = "AUTO")
plot_grid(p1, p2, labels = "AUTO", rel_widths = c(1, 1.5))



library(e1071)  # 加载e1071包,它包含SVM算法
library(pROC)
# 准备数据
# 假设你有一个包含特征蛋白数据的数据框，其中包括特征列和目标列
# 特征列是一些数值型特征的值，目标列是类别标签（例如0或1）
# 假设数据框名为"protein_data"，特征列名为"feature1"、"feature2"等，目标列名为"target"
data_raw[is.na(data_raw)]<-0.0001

protein_data<-data.frame(scale(t(data_hela_a549)))

protein_data$label<-as.factor(hela_a549_sample$Group)
# 划分数据集为训练集和测试集
set.seed(123)  # 设置随机种子以确保结果可复现

train_indices <- sample(nrow(protein_data), nrow(protein_data) * 0.7)  # 随机选择80%的数据作为训练集
train_data <- protein_data[train_indices, ]  # 训练集数据
test_data <- protein_data[-train_indices, ]  # 测试集数据

tune.out=tune(svm ,label~.,data=train_data ,kernel ="linear",
              ranges =list(cost=c(0.001,0.01,0.1, 1,5,10,100)))
summary(tune.out)
bestmod =tune.out$best.model
support_vector_indices <- bestmod$index  # 获取支持向量的索引
support_vectors <- train_data[support_vector_indices, ]  # 提取支持向量
data_med[is.na(data_med)]==0.01
prediction <- predict(bestmod, newdata = data.frame(t(log2(data_med)))[1:5,])
plot(bestmod)

bestmod
# 在测试集上进行预测
predictions <- predict(bestmod, newdata = test_data)
predictions <- predict(bestmod, newdata = test_data, probability = TRUE)
# 创建混淆矩阵
confusion_mat <- table(Actual = test_data$label, Predicted = predictions)
print(confusion_mat)
# 计算准确率
accuracy <- sum(diag(confusion_mat)) / sum(confusion_mat)
print(paste("Accuracy:", accuracy))
# 加载ROCR包
library(ROCR)
data(ROCR.simple)
pred <- prediction(ROCR.simple$predictions,ROCR.simple$labels)
# 创建一个prediction对象
pred <- prediction(predictions, labels)

# 计算真正例率和假正例率
perf <- performance(pred, "tpr", "fpr")


# 绘制ROC曲线
plot(perf, main = "ROC Curve", xlab = "False Positive Rate", ylab = "True Positive Rate")

# 计算AUC
auc_obj <- auc(roc_obj)
auc_value <- as.numeric(auc_obj)
print(paste("AUC:", auc_value))
# 构建SVM分类器模型
svm_model <- svm(label ~ ., data = train_data, kernel = "linear", cost = 1)  # 使用线性核函数和成本参数1构建SVM模型

# 在测试集上进行预测
predictions <- predict(bestmod, newdata = test_data)

# 评估分类器性能
accuracy <- sum(predictions == final_data$label) / nrow(final_data)  # 计算准确率

# 在38T1上进行预测
predictions <- data.frame(predict(svm_model, newdata =final_data))

fisher.test(predictions$predict.svm_model..newdata...final_data.,metadata_T1_38_1$Response_res)


OS_38<-km_OS_PLOT_cluster(cluster =predictions,
                          meatadata = metadata_T1_38_1,
                          title=c("OS"))











###hallmark,KEGG,reactome
library(msigdbr)
hall <- msigdbr(species = "Homo sapiens", category = "H")
kegg <-msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:KEGG")
reactome <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:REACTOME")

merge_term<-rbind(hall,kegg,reactome)

merge_term_list <- merge_term %>%
  split(.$gs_name) %>%
  lapply(function(df) {
    df$gene_symbol
  })

names(merge_term_list) <- merge_term$gs_name[!duplicated(merge_term$gs_name)]








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
aa<-cell_specific_data(data_express=data_hela_a549,group_data=hela_a549_sample,p_value=0.01,method="BH",Fc_value=2)
A549_up_clus<-clustern_kegg_go_list(aa[[1]]$DEG_up)
Hela_up_clus<-clustern_kegg_go_list(aa[[2]]$DEG_up)

HepG2_up_clus<-clustern_kegg_go_list(aa[[3]]$DEG_up)

A549_down_clus<-clustern_kegg_go_list(aa[[1]]$DEG_down)
Hela_down_clus<-clustern_kegg_go_list(aa[[2]]$DEG_down)
HepG2_down_clus<-clustern_kegg_go_list(aa[[3]]$DEG_down)

A549_up<-aa[[1]]$DEG_up
A549_plot_up_ppi<-plot_stringdb(A549_up,layout = "kk",deg_filter=0,fil_nod = 0)



Hela_up<-aa[[2]]$DEG_up
Hela_plot_up_ppi<-plot_stringdb(Hela_up,layout = "kk",deg_filter=0,fil_nod = 0)

HepG2_up<-aa[[3]]$DEG_up

HepG2_plot_up_ppi<-plot_stringdb(HepG2_up,layout = "kk",deg_filter=1,fil_nod = 1)
pdf(file="./analysis/hela_a549_Hepg2/20231016_A549_up_Clus.pdf",width = 9,height = 7.5)
A549_plot_up_ppi
dev.off()
pdf(file="./analysis/hela_a549_Hepg2/20231016_Hela_up_Clus.pdf",width = 9,height = 7.5)
Hela_plot_up_ppi
dev.off()

pdf(file="./analysis/hela_a549_Hepg2/20231016_HepG2_up_Clus.pdf",width = 9,height = 7.5)
HepG2_plot_up_ppi
dev.off()
write.csv(A549_up,file="./analysis/hela_a549_l02_HepG2g2/A549_up.csv")
write.csv(Hela_up,file="./analysis/hela_a549_l02_HepG2g2/Hela_up.csv")
write.csv(L02_up,file="./analysis/hela_a549_l02_HepG2g2/L02_up.csv")
write.csv(HepG2G2_up,file="./analysis/hela_a549_l02_HepG2g2/HepG2G2_up.csv")

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
cell_specific_data=phos_NMF_sample
extrat_up_Genes<-function(cell_specific_data){
  n<-length(cell_specific_data)
  DEG<-c()
  for(i in 1:n){
    
    DEG<-c(DEG,cell_specific_data[[i]]$DEG_up)
  }
  DEG<- unique(DEG)
  return(DEG)
  
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


