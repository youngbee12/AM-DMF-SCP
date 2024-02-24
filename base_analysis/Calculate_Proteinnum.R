########################
######导入数据##########
########################

load(file="./data/20231008_scdata_pre.Rdata")
library(tidyverse)


sc_group
nacol<-data.frame(nrow(sc_data) - colSums ( is.na ( sc_data ) ) )
colnames(nacol)<-"p_num"
sc_group$p_num <- nacol$p_num

filter_data <- sc_group %>% group_by(Group) %>%
  arrange(-p_num) %>%
  filter(n()<10|row_number() <= n() - 2) %>%
  ungroup()

median_data<-filter_data %>%group_by(Group) %>%
  summarize(mean_value = median(p_num))
mean_data<-filter_data %>%group_by(Group) %>%
  summarize(mean_value = mean(p_num))

wide_df <- filter_data %>% 
  pivot_wider(names_from = Group, values_from = p_num) %>% select(-Rep) %>% t(.) %>% data.frame(.)
mean_num<-apply(wide_df,1,function(x)mean(x,na.rm=T))
sd_num<-apply(wide_df,1,function(x)sd(x,na.rm=T))
wide_df$Group<-rownames(wide_df)
wide_df$mean<-mean_num
wide_df$sd<-sd_num
ggplot(wide_df, aes(x=Group, y=mean,color=Group,fill=Group)) +
  geom_bar(position=position_dodge(), stat="identity",colour='white',width=.6,alpha=0.2) +
  # scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.4)+
  #scale_color_manual(values=c("#24C9D4","#324B4D","#94B0B3","#B3AFF7","#7C7BBE"))+
  #scale_fill_manual(values=c("#24C9D4","#324B4D","#94B0B3","#B3AFF7","#7C7BBE"))+
  geom_jitter(data = filter_data,aes(x=Group,y=p_num),width=.1)+
  scale_y_continuous(limits=c(0,3000),breaks=seq(0,3000,500))+
  geom_text(aes(label = ceiling(mean)), vjust = 5, colour = "black")+
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




