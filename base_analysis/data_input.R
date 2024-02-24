########################
######导入数据##########
########################
library(readxl)
library(dplyr)
data <- data.frame(read_excel("./data/diann_result1.xlsx",sheet = "Sheet1",col_names = TRUE))

data1 <-data %>% `rownames<-`(.[,1]) %>% dplyr::select(-Genes)

split_string <- unlist(lapply(colnames(data1),function(x)unlist(strsplit(x, "yzc_"))[2]))
split_string1 <- unlist(lapply(split_string,function(x)unlist(strsplit(x, "\\."))[1]))
colnames(data1) <- split_string1



group_info <- data.frame(t(data.frame(lapply(split_string1, function(x)strsplit(sub("^(_*[^_]*)_", "\\1#", x), "#")))))

rownames(group_info) <-split_string1
group_info$sample<-rownames(group_info)
nacol<-data.frame(nrow(data1) - colSums ( is.na ( data1 ) ) )
colnames(nacol)<-"p_num"
group_info$p_num <- nacol$p_num
A<-c(rep("batch1",6),rep("batch2",6),rep("batch3",6),rep("batch1",6),rep("batch2",6),rep("batch3",6),rep("batch1",6),rep("batch3",6),rep("batch2",6),rep("batch1",6),rep("batch1",6),rep("batch2",6),rep("batch3",6),rep("batch1",6),rep("batch2",5),rep("batch3",6),rep("batch1",6),rep("batch2",6),rep("batch3",6),rep("QC",6))
group_info$batch <-A
colnames(group_info) <- c("Group","Rep","sample","p_num","batch") 
group_info$batch_sample<-paste(group_info$Group,group_info$batch,sep="_")

data1<-median_correct(data1)
sc_data <- data1[,1:113]
qc_data <- data1[,114:119]
sc_group <- group_info[1:113,]
qc_group <- group_info[114:119,]
wide_df <- group_info %>% 
  pivot_wider(names_from = Group, values_from = p_num) %>% dplyr::select(-Rep,-sample,-batch,-batch_sample) %>% t(.) %>% data.frame(.)
colnames(wide_df)<-rownames(group_info)
mean_num<-apply(wide_df,1,function(x)mean(x,na.rm=T))
sd_num<-apply(wide_df,1,function(x)sd(x,na.rm=T))
wide_df$Group<-rownames(wide_df)
wide_df$mean<-mean_num
wide_df$sd<-sd_num
save(sc_data,sc_group,qc_data,qc_group,wide_df,file="./data/20231011_scdata_pre.Rdata")
# 输入字符串
#s <- "aa/sd/asd/2131"

# 将第二个斜杠替换为一个特殊字符
#s <- sub("^([^/]*/[^/]*)/", "\\1#", s)

# 使用特殊字符进行分割
#result <- strsplit(s, "#")

# 打印结果
#print(unlist(result))
#在这个示例中，sub函数的第一个参数是一个正则表达式，用于匹配字符串开头的两个斜杠之间的部分
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
