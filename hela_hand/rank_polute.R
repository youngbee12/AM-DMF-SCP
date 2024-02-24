

data_pro_select_hela_hand<-data_genes[,hela_hand_sample$sample]



n<-length(unique(hela_hand_sample$Group))
for (i in 1:n){
  t_sm<-hela_hand_sample[hela_hand_sample$Group==unique(hela_hand_sample$Group)[i],]$sample
  listn[i]<-list(data_pro_select_hela_hand[,colnames(data_pro_select_hela_hand) %in% t_sm] %>%
                   .[ rowSums ( is.na ( . ) )  < ncol(.)*0.5,] %>% rownames(.))
}

inchip_sel<-data_pro_select_hela_hand[listn[[1]],1:6]
inchip_sel$mean<-apply(inchip_sel,1,function(x)mean(x,na.rm=T))
inchip_sel <- inchip_sel %>% arrange(desc(mean))
inchip_sel$rank <-seq(1,nrow(inchip_sel),1)
inchip_sel$inte<-log10(inchip_sel$mean)
inchip_sel$Name<-rownames(inchip_sel)
inchip_sel <- inchip_sel %>% mutate(Polu=ifelse(grepl("CON_",Name),"Con","not_Con"))
length(inchip_sel[inchip_sel$Polu=="CON_",])
length(outchip_sel[outchip_sel$Polu=="CON_",])
inchip_sel <- inchip_sel %>% mutate(Polu=ifelse(Name %in% polute_pro_uniprot,"Con",Polu))

rank_inchip<-ggplot(inchip_sel, aes(x = rank, y = inte,color=Polu,alpha=Polu)) + 
  geom_point( size = 1) +
  #geom_boxplot(data = tp53_df, aes(x = rank, y = values),
  #width = 500, color = "red", fill = NA) +
  # geom_jitter(data = tp53_df, aes(x = rank, y = values), width = 0.2, size = 1, color = "black") + # 添加散点
  scale_alpha_manual(values = c(1, 0.1)) +
  ylim(2,8)+
  xlab("")+
  ylab("Log10(Intensity)")+
  geom_text(aes(x = 1000, y = 6, label = paste0("Inchip")),
            hjust = 0, vjust = 1, color = "black", size = 8) + 
  #scale_color_manual(values = c("Con" = "red", "not_Con" = "white")) +
  #labs(title = "Scatter plot of Rank vs Log2 Transformed Pro Tumor Median with Boxplot for TP53",
  #x = "Rank", y = "Log10 Transformed Pro Tumor Median", color = "Legend") +
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


outchip_sel<-data_pro_select_hela_hand[listn[[2]],7:12]
outchip_sel$mean<-apply(outchip_sel,1,function(x)mean(x,na.rm=T))
outchip_sel <- outchip_sel %>% arrange(desc(mean))
outchip_sel$rank <-seq(1,nrow(outchip_sel),1)
outchip_sel$inte<-log10(outchip_sel$mean)
outchip_sel$Name<-rownames(outchip_sel)
outchip_sel <- outchip_sel %>% mutate(Polu=ifelse(grepl("CON_",Name),"Con","not_Con"))
outchip_sel <- outchip_sel %>% mutate(Polu=ifelse(Name %in% polute_pro_uniprot,"Con",Polu))
rank_outchip<-ggplot(outchip_sel, aes(x = rank, y = inte,color=Polu,alpha=Polu)) + 
  geom_point( size = 1) +
  #geom_boxplot(data = tp53_df, aes(x = rank, y = values),
  #width = 500, color = "red", fill = NA) +
  # geom_jitter(data = tp53_df, aes(x = rank, y = values), width = 0.2, size = 1, color = "black") + # 添加散点
  scale_alpha_manual(values = c(1, 0.1)) +
  ylim(2,8)+
  xlab("")+
  ylab("Log10(Intensity)")+
  geom_text(aes(x = 750, y = 6, label = paste0("Outchip")),
            hjust = 0, vjust = 1, color = "black", size = 8) + 
  #scale_color_manual(values = c("Con" = "red", "not_Con" = "white")) +
  #labs(title = "Scatter plot of Rank vs Log2 Transformed Pro Tumor Median with Boxplot for TP53",
       #x = "Rank", y = "Log10 Transformed Pro Tumor Median", color = "Legend") +
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


library(patchwork)
pdf(file="analysis/hela_hand/20231020_polute_rank_hela_hand_inchip_outchip.pdf",width=6,height = 5)
rank_inchip/rank_outchip
dev.off()
pdf(file="analysis/hela_hand/20231020_polute_rank_hela_hand_inchip.pdf",width=5,height = 5)
rank_inchip
dev.off()
pdf(file="analysis/hela_hand/20231020_polute_rank_hela_hand_outchip.pdf",width=5,height = 5)
rank_outchip
dev.off()
# 获取TP53的数据
tp53_rank <- data_test_log$rank[data_test_log$protein == "TP53"]
tp53_values <- unlist(data_test_log[data_test_log$protein == "TP53", 5:ncol(data_test_log)])
tp53_df <- data.frame(rank = rep(tp53_rank, times = length(tp53_values)), values = tp53_values)
# 绘制散点图
# 绘制散点图
ggplot(data_test_log, aes(x = rank, y = pro_tumor_median)) + 
  geom_point(aes(color = "Pro Tumor Median"), size = 1.5, alpha = 0.6) +
  geom_boxplot(data = tp53_df, aes(x = rank, y = values),
               width = 500, color = "red", fill = NA) +
  geom_jitter(data = tp53_df, aes(x = rank, y = values), width = 0.2, size = 1, color = "black") + # 添加散点
  
  scale_color_manual(values = c("Pro Tumor Median" = "blue", "TP53" = "red")) +
  labs(title = "Scatter plot of Rank vs Log2 Transformed Pro Tumor Median with Boxplot for TP53",
       x = "Rank", y = "Log2 Transformed Pro Tumor Median", color = "Legend") +
  theme_minimal()





library(Biostrings)
library(dplyr)
sequences <- readBStringSet("analysis/external_data/contaminants.fasta")
sq<-data.frame(sequences)
sq$name<-rownames(sq)
se1<-sq
extract_gene_name <- function(str) {
  return(sub(".*GN=([^ ]+).*", "\\1", str))
}

gene_names <- sapply(sq$name, extract_gene_name)
se1$gene_names<-gene_names
polute_pro_uniprot<-substr(rownames(se1),1,6)
