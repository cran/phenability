fox <- function(data){
   
  a <- nrow(data)
  b <- ncol(data)
  data.m <- as.matrix(data)
  l.data <- length (data.m)  
  ranks <- matrix(NA,a,b)
  y <- NULL
  k <- matrix(NA)
  k2 <- matrix(NA) 
  for(i in 1:nrow(data.m)){
    for (j in 1:ncol(data.m)){
      
      ranks[,j] <- rank(-data.m[,j])
    }
    
    y <- which(ranks[i,] <= 3)
    k[i] <- length(y)
    k[i] <- as.data.frame(k[i])
    k2[i] <- as.matrix(k[i])
    k2 <- as.numeric(k2[i])
    
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
  result <- cbind(rownames(data),means,k)
  colnames(result) <- c("Gen","Mean", "TOP")
  
  result.g <- cbind(rownames(data),means,k2)
  
  top <- as.numeric(result[,3])
  plot(means,top,pch=19,cex=0.5,main="Means x TOP",xlab=expression(Mean[Phenotypic]),
       ylab=expression(TOP[third]),xlim=c(min(means),max(means)),
       ylim=c(min(top),max(top)))
  m <- apply(cbind(means,top),2,mean)
  textxy(means,top,1:a,m=m,cex=1, col="blue")
  origin(m)
  
  return(result)
  
}
