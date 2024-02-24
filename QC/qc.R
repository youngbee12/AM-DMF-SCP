



#相关性

panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
   {
       usr <- par("usr"); on.exit(par(usr)) 
       par(usr = c(0, 1, 0, 1))
      r <- abs(cor(x, y,use="pairwise.complete.obs")) 
       txt <- format(c(r, 0.123456789), digits=digits)[1] 
       txt <- paste(prefix, txt, sep="") 
       if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
       
           test <- cor.test(x,y) 
           # borrowed from printCoefmat
             Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                                                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                                                  symbols = c("***", "**", "*", ".", " ")) 
             
               text(0.5, 0.5, txt, cex = cex * r) 
             text(.8, .8, Signif, cex=cex, col=2) 
           }
 pdf('./analysis/QC/cor.pdf', width=5,height=5)
 
 pairs(log2(median_correct(qc_data)), upper.panel=panel.cor, gap=0)
 dev.off()
cor(log2(median_correct(qc_data)),use="pairwise.complete.obs")   #temp为相关性系数

#CV曲线
library(rknn)
cv_qc<-apply(qc_data,1,function(x)cv.coef(as.numeric(x)))
pdf('./analysis/QC/cv.pdf', width=5,height=5)
plot(density(na.omit(cv_qc)),main = NA, xlab = "CV", ylab = "Density",cex=2) 
abline(v = 0.135, col = "red", lty = 2)# 添加标
# 添加文本标签
text(x = 0.3, y = 4, labels = "Median(CV) = 0.1536757", col = "blue",cex=1)
dev.off()
median(cv_qc,na.rm = T)


#BOXPLOT


