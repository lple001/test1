myDescriptStat <- function(x){
  #x<-pheno_stem[,2]
  #x<-as.numeric(x)
  #x<-impute(x,mean)
  x<-na.omit(x)
  n <- length(x)                    #样本数据个数
  m <- mean(x,na.rm=T)                      #均值
  me <- median(x,na.rm=T)                   #中位数
  mo <- names(table(x))[which.max(table(x))]  #众数
  sd <- sd(x,na.rm=T)                       #标准差
  v <- var(x,na.rm=T)                       #方差
  maxx<-max(x,na.rm=T)
  minn<- min(x,na.rm=T)
  r <- max(x,na.rm=T) - min(x,na.rm=T)              #极差
  cv <- 100 * sd/m                  #变异系数
  css <- sum(x - m)^2               #样本校正平方和
  uss <- sum(x^2)                   #样本未校正平方和
  R1 <- quantile(x,0.75,na.rm=T) - quantile(x,0.25,na.rm=T)     #四分位差
  se <- sd/sqrt(n)                              #标准误
  
  g1 <- n/((n-1)*(n-2)*sd^3)*sum((x-m)^3)/sd^3  #偏度系数
  g2 <- ((n*(n+1))/((n-1)*(n-2)*(n-3))*sum((x-m)^4)/sd^4 -(3*(n-1)^2)/((n-2)*(n-3))) #峰度系数
  
  data.frame(N=n,Mean=m,Median=me,Mode=mo,
             Std_dev=sd,Variance=v,max=maxx,min=minn,Range=r,
             CV=cv,CSS=css,USS=uss,
             R1=R1,SE=se,Skewness=g1,Kurtosis=g2,
             row.names=1)
}


panel.lm <- function (x, y, col, bg = NA, pch = par("pch"),
                      cex =3, ...) {
  p<- as.numeric(cor.test(x, y, method="pearson")[3])
  
  if(p>0.05){points(x, y, pch = pch, col ="gray", bg = bg, cex = cex)
    abline(stats::lm(y ~ x), col= "gray", ...)}
  if(p<0.05){points(x, y, pch = pch, col =col, bg = bg, cex = cex)
    abline(stats::lm(y ~ x), col= "red", ...)}
}
panel.nlm <- function (x, y, col , bg = NA, pch = par("pch"),
                       cex =3, ...) {
  p<- as.numeric(cor.test(x, y, method="pearson")[3])
  if(p>0.05){points(x, y, pch = pch, col ="gray", bg = bg, cex = cex)
    lines(stats::lowess(y ~ x, f = 5, iter = 2), col= "gray", ...)}
  if(p<0.05){points(x, y, pch = pch, col =col, bg = bg, cex = cex)
    lines(stats::lowess(y ~ x, f = 5, iter = 2), col= "red", ...)}
}
panel.cor <- function(x, y, digits=2, prefix="", cex.cor,col, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r1 <- as.numeric(cor.test(x, y, method="pearson")[4]);r<-abs(r1)
  p<- as.numeric(cor.test(x, y, method="pearson")[3])
  txt <- format(c(r1, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  if(p>0.01 & p<=0.05){text(0.5, 0.6,  c("*"),col=col,  cex =  1.5);text(0.5, 0.4,  txt,col=col,  cex =  1.5)}
  if(p>0.001 & p<=0.01){text(0.5, 0.6,  c("**"),col=col,  cex =  1.5);text(0.5, 0.4,  txt,col=col,  cex =  1.5)}
  if(p>0.0001 & p<=0.001){text(0.5, 0.6,  c("***"),col=col,  cex =  1.5);text(0.5, 0.4,  txt,col=col,  cex =  1.5)}
  if(p<=0.0001){text(0.5, 0.6,  c("****"),col=col,  cex =  1.5);text(0.5, 0.4,  txt,col=col,  cex =  1.5)}
  if(p>0.05){text(0.5, 0.6, c("NS"),col="gray", cex =  1.3);text(0.5, 0.4,txt,col="gray", cex =  1.5)}
}

panel.hist <- function(x, col, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5)) #
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks
  nB <- length(breaks)
  y <- h$counts
  y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col=col,font=3, ...)#col="cyan", ...)
}


phenocheck.func<-function(pheno.file){

  if (!requireNamespace("Hmisc", quietly = TRUE)) {
      install.packages("Hmisc")
  }

  library(Hmisc)
  
  if(!exists("./datacheck.result")) dir.create("./datacheck.result")
  
  pheno0<-read.csv(pheno.file,head=T)
  file.copy(pheno.file, "./datacheck.result/pheno.csv", overwrite = TRUE)
  file.copy(pheno.file, "./datacheck.result/phenofortxt.csv", overwrite = TRUE)
  file.rename("./datacheck.result/phenofortxt.csv","./datacheck.result/pheno.txt")

  pheno1<-pheno0[,-1]
  pheno1<-as.data.frame(lapply(pheno1,as.numeric))
  col_name<-colnames(pheno1)
  pheno_sum0<-apply(pheno1, 2, myDescriptStat)
  row_name<-names(pheno_sum0[[1]])
  pheno_sum1 <- data.frame(t(matrix(unlist(pheno_sum0), nrow=dim(pheno1)[2], byrow=T)))
  rownames(pheno_sum1)<-row_name
  colnames(pheno_sum1)<-col_name
  write.csv(pheno_sum1,"./datacheck.result/pheno_sum.csv")
  
  pairsdf<-pheno1
  
  labels<-colnames(pairsdf)
  
  
  
  pdffile<-"./datacheck.result/Pheno.pairs.plot.pdf"
  pdf(pdffile, length(col_name)*2,length(col_name)*2);
  #op <- par(xaxt="n", yaxt="n")
  pairs(pairsdf, pch=".",
        upper.panel=panel.cor,
        diag.panel = panel.hist,
        lower.panel = panel.lm,
        cex.labels=2, 
        labels=labels,
        col="#4169E1")
  #par(op);
  dev.off()
}

