########################
######导入数据##########
########################

library(PerformanceAnalytics)
library(corrplot)
library(RColorBrewer)
pdf(file = "./result/cor.pdf",width=15)
corrplot(cor(sc_data,method="pearson"),
         method="color",tl.col = "black"
        
         #col=c(rep("grey", 200), COL2('BrBG', 200)),
         #col.lim=c(0,1),
         #addCoef.col = "white",
         #type="upper")
)



dev.off()
ax<-data.frame(cor(log2(sc_data),method="pearson",use="complete.obs"))
corrplot(cor(log2(data_med),method="pearson",use="complete.obs"),
         method="color",tl.col = "black",
         
         col=c(rep("white", 700), COL2('RdYlBu', 300)),
         col.lim=c(0.5,1),
         #addCoef.col = "white",
         #type="upper")
)
  ########################
########用到函数########
########################
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

