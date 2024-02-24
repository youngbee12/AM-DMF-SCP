########################
######导入数据##########
########################
data_depmap<-read.csv("./analysis/external_data/DepMap_data/OmicsExpressionProteinCodingGenesTPMLogp1.csv")
Depmap_anno<-read.csv("./analysis/external_data/DepMap_data/Model.csv")
colna<-colnames(data_depmap)
colna1<-unlist(lapply(colna,function(x)unlist(strsplit(x,"\\.."))[1]))
colnames(data_depmap)<-colna1
colnames(data_depmap)[1]<-"ModelID"
rownames(data_depmap) <-data_depmap[,"ModelID"]
raw_depmap<-data_depmap %>% dplyr::select("ModelID") %>% left_join(Depmap_anno)


#A549=ACH-000681
#Hela=ACH-001086
#HepG2=ACH-000739
select_depmap_Data<-data_depmap[c("ACH-000681","ACH-001086","ACH-000739"),]
rownames(select_depmap_Data)<-c("A549","Hela","HepG2")
select_depmap_Data<-data.frame(t(select_depmap_Data[,-1]))
select_depmap_Data1<-select_depmap_Data
select_depmap_Data1$a<-rownames(select_depmap_Data1)
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

