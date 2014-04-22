nahu <-
function(data){
 
  
  a <- nrow(data)
  b <- ncol(data)
  data.m <- as.matrix(data)
  l.data <- length (data.m)
  data.r <- matrix(NA,a,b) 
  S2.1 <- matrix(NA,a,b)
  S3.1 <- matrix(NA,a,b)
  S6.1 <- matrix(NA,a,b)
  S2 <- numeric()
  S3 <- numeric()
  S6 <-  numeric()
  S.1 <- matrix(NA,a,b)
  S1 <- numeric()
  k <-  2
  
  for (i in 1:nrow(data.m)){
    for(j in 1:ncol(data.m)){
      data.r[i,j] <- (data.m[i,j]) - (mean(data.m[i,]))+(mean(data.m))
    }
  }
  ranks <- matrix(NA,a,b)
  
  for (j in 1:ncol(data.r)){
    
    ranks[,j] <- rank(data.r[,j])
  }
  
  ranks.y <- matrix(NA,a,b)
  
  for (j in 1:ncol(data.m)){
    
    ranks.y[,j] <- rank(data.m[,j])
  }
  
  r.means <- rowMeans(ranks)
  r.means.y <- rowMeans(ranks.y)
  
  for (i in 1:nrow(data)){
    for (j in 1:ncol(data)){
      S2.1[i,j] <- (ranks[i,j]-r.means[i])^2
      S3.1[i,j] <- (ranks.y[i,j]-r.means.y[i])^2 
      S6.1[i,j] <- (abs(ranks.y[i,j]-r.means.y[i]))
    }
    S2[i]<-as.data.frame((sum(S2.1[i,])) / (b-1))
    S3[i]<-as.data.frame(sum(S3.1[i,]) / (r.means.y[i]))
    S6[i]<-as.data.frame(sum(S6.1[i,]) / (r.means.y[i]))
  }
  for (i in 1:nrow(data)){
    for (j in 1:(b-1)){
      
      S.1[i,j] <- abs(ranks[i,j] - ranks[i,k]) 
      while(k < b)
        k <- k + 1
    }
    S1[i] <- (2*(sum(S.1[i,j]))) / (b*(b-1))
    S1[i] <- as.data.frame(S1[i])
  }
  hist(data.m, prob=T,breaks=sqrt(l.data), xlab="Data",main="Diagnostic Graph", border=c("green"))
  densi.d <- density(data.m)
  lines(density(data.m,bw="SJ-ste"))
  set.seed(l.data)
  z <- rnorm(l.data, mean=mean(data.m), sd=sd(data.m))
  densi.s <- density(z)
  lines(density(z) , col=2)
  legend(summary(densi.d$x)[5],max(densi.d$y),  c('data', 'normal.dist'), lty=1,col=c(1,2))
  
  means <- as.vector(rowMeans(data))
  result <- cbind(rownames(data),means,S1,S2,S3,S6)
  colnames(result) <- c("Gen","Mean","S1","S2","S3","S6")
  
  s1 <- as.numeric(result[,3])
  s2 <- as.numeric(result[,4])
  s3 <- as.numeric(result[,5])
  s6 <- as.numeric(result[,6])
  
  plot(means,s1,pch=19,cex=0.5,main="Means x S1",xlab=expression(Mean[Phenotypic]),
       ylab=expression(S[1]),xlim=c(min(means),max(means)),
       ylim=c(min(s1),max(s1)))
  m <- apply(cbind(means,s1),2,mean)
  textxy(means,s1,1:a,m=m,cex=1,col="blue")
  origin(m)
  plot(means,s2,pch=19,cex=0.5,main="Means x S2",xlab=expression(Mean[Phenotypic]),
       ylab=expression(S[2]),xlim=c(min(means),max(means)),
       ylim=c(min(s2),max(s2)))
  m <- apply(cbind(means,s2),2,mean)
  textxy(means,s2,1:a,m=m,cex=1,col="blue")
  origin(m)
  plot(means,s3,pch=19,cex=0.5,main="Means x S3",xlab=expression(Mean[Phenotypic]),
       ylab=expression(S[3]),xlim=c(min(means),max(means)),
       ylim=c(min(s3),max(s3)))
  m <- apply(cbind(means,s3),2,mean)
  textxy(means,s3,1:a,m=m,cex=1,col="blue")
  origin(m)
  plot(means,s6,pch=19,cex=0.5,main="Means x S6",xlab=expression(Mean[Phenotypic]),
       ylab=expression(S[6]),xlim=c(min(means),max(means)),
       ylim=c(min(s6),max(s6)))
  m <- apply(cbind(means,s6),2,mean)
  textxy(means,s6,1:a,m=m,cex=1,col="blue")
  origin(m)
  
  return(result)
}
