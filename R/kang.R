kang <- 
function (data){
  
  a <- nrow(data)
  b <- ncol(data)
  data.m <- as.matrix(data)
  l.data <- length (data.m)
  yield <- numeric()
  vars <- numeric()
  
  for (i in 1:nrow(data.m)){
    
    yield[i] <- sum(data.m[i,])  
    vars[i] <- var(data.m[i,])
  }
  
  rank.y <- rank(-yield)
  rank.v <- rank(vars)
  rank.sum <- rank.y + rank.v
  rank.sum <- as.data.frame(rank.sum)
  
  hist(data.m, prob=T,breaks=sqrt(l.data), xlab="Data",main="Diagnostic Graph", border=c("green"))
  densi.d <- density(data.m)
  lines(density(data.m,bw="SJ-ste"))
  set.seed(l.data)
  z <- rnorm(l.data, mean=mean(data.m), sd=sd(data.m))
  densi.s <- density(z)
  lines(density(z) , col=2)
  legend(summary(densi.d$x)[5],max(densi.d$y),  c('data', 'normal.dist'), lty=1,col=c(1,2))
  
  means <- as.vector(rowMeans(data))
  result <- cbind(rownames(data),means,rank.y,rank.v,rank.sum)
  colnames(result) <- c("Gen","Mean","(Rank PROD)","(Rank VAR)","(RS-rank sum)")
  
  rs <- as.numeric(result[,5])
  plot(result[,2],result[,5],pch=19,cex=0.5,main="Means x Rank Sum",
       xlab=expression(Mean[phenotypic]),ylab=expression(rs["(rank sum)"]),
       xlim=c(min(result[,2]),max(result[,2])),
       ylim=c(min(result[,5]),max(result[,5])))
  m <- apply(cbind(result[,2],result[,5]),2,mean)
  textxy(means,rs,1:a,m=m,cex=1,col="blue")
  origin(m)
  
  
  return(result)
}
