########################
######导入数据##########
########################

### Histogram 


### Set colors
color_all <- "black"
color_both <- "blue"
transparency <- 0.1

### Generate plot

  proteins_all <- read.csv("analysis/external_data/Protein_abs_copynum/HeLaHF14k.csv")
  rownames(proteins_all)<-proteins_all$Accession
  gene<-data.frame(proteins_all$Accession)
  print(dim(gene)[1])
  colnames(gene)<-"gene"
  DEA<- clusterProfiler::bitr(gene$gene,fromType="UNIPROT", toType="SYMBOL", OrgDb="org.Hs.eg.db")
  DEA<-distinct(DEA,SYMBOL,.keep_all = TRUE)
  
  proteins_all_final<-proteins_all[DEA$UNIPROT,]
  proteins_all_final$Symbol<-DEA$SYMBOL
  row.names(DEA)<-DEA$SYMBOL
  csv2 <- listn[[1]]
  csv3 <- listn[[2]]
  ### Get the overlap
  proteins_both <- proteins_all_final %>%
    dplyr::filter(Symbol %in% csv2)
  proteins_both2 <- proteins_all_final %>%
    filter(Symbol %in% csv3)
  ### Histogram
  pdf(file="analysis/hela_hand/20231111_Copynum.pdf",width = 5,height = 5)
  ggplot(proteins_all_final) +
    xlab("Log")+
    ylab("Count")+
    geom_histogram(
      data = proteins_all_final,
      aes_string(x = "Log"), 
      color = color_all,
      fill = color_all,
      alpha = transparency,
      bins = 20
    ) +geom_vline(xintercept=mean(proteins_all_final$Log),linetype = "dashed")+
    geom_histogram(
      data = proteins_both,
      aes_string(x = "Log"), 
      color = "#4682B4",
      fill = "#4682B4",
      alpha = transparency,
      bins = 20
    )+geom_vline(xintercept=mean(proteins_both$Log),color="#4682B4",linetype = "dashed")+
    geom_histogram(
      data = proteins_both2,
      aes_string(x = "Log"), 
      color = "#AF46B4",
      fill = color_both,
      alpha = transparency,
      bins = 20
    ) +geom_vline(xintercept=mean(proteins_both2$Log),color="#AF46B4",linetype = "dashed")+
    
    
    geom_text(aes(x =4, y = 1000, label = paste0("Mean(All) = ", round(mean(proteins_all_final$Log), 2))),
                    hjust = 0, vjust = 1, color = "black", size = 6) +  
    geom_text(aes(x =4, y = 1000, label = paste0("\nMean(Inchip)= ",round(mean(proteins_both$Log), 2))),
                    hjust = 0, vjust = 1, color = "#4682B4", size =6) +  
    geom_text(aes(x =4, y = 1000, label = paste0("\n\nMean(Outchip)= ",round(mean(proteins_both2$Log), 2))),
                    hjust = 0, vjust = 1, color = "#AF46B4", size = 6) +
    theme_bw() +
    
    theme(  #legend.position = "none",
            #axis.line = element_line(size=1.5),
            axis.text.x = element_text(color="black", size=16, face="plain",angle=0,hjust=0.5) ,# Values for face are one of "plain", "italic", "bold" and "bold.italic"
            axis.text.y = element_text(color="black", size=16, face="plain"),
            axis.title.x = element_text(color="black", size=16, face="bold"),
            axis.title.y = element_text(color="black", size=16, face="bold"),
            legend.text=element_text(color="black", size=16, face="bold"),
            legend.title=element_text(color="black", size=16, face="bold"),
            title = element_text(color="black", size=16, face="bold"))

    
dev.off()
dev.new()
### Density 

library(ggrepel)
### Set colors
color_all <- "black"
color_both <- "blue"
transparency <- 0.1

### Generate plot

  ### Load files
  file <- input$file1
  ext <- tools::file_ext(file$datapath)
  req(file)
  validate(need(ext == "csv", "Please upload a csv file"))
  proteins_all <- read.csv(file$datapath)
  file <- input$file2
  ext <- tools::file_ext(file$datapath)
  req(file)
  validate(need(ext == "csv", "Please upload a csv file"))
  csv2 <- read.csv(file$datapath)
  
  ### Get the overlap
  proteins_both <- proteins_all %>%
    filter(Accession %in% csv2$Accession)
  
  ### Histogram
  ggplot(proteins_all) +
    geom_density(
      data = proteins_all,
      aes_string(x = input$variable), 
      color = color_all,
      fill = color_all,
      alpha = transparency,
      bins = input$bins
    ) +
    geom_density(
      data = proteins_both,
      aes_string(x = input$variable), 
      color = color_both,
      fill = color_both,
      alpha = transparency,
      bins = input$bins
    ) +
    labs(
      title = paste0("Histogram of ", input$variable)
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  
  ggplot(proteins_all_final) +
    geom_density(
      data = proteins_all_final,
      aes_string(x = "Log"), 
      color = color_all,
      fill = color_all,
      alpha = transparency,
      bins = 20
    ) +
    geom_density(
      data = proteins_both,
      aes_string(x = "Log"), 
      color = color_both,
      fill = color_both,
      alpha = transparency,
      bins = 20
    )+
    geom_density(
      data = proteins_both2,
      aes_string(x = "Log"), 
      color = "red",
      fill = color_both,
      alpha = transparency,
      bins = 20
    ) +
    labs(
      title = paste0("Histogram of ")
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  
  
