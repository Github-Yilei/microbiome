get_quant <- function(xx, yy, q) { xx[which(cumsum(yy)/sum(yy) >=q)][[1]] }
means <- data.frame(rowsum(qs$Score*qs$Count, qs$Cycle)/rowsum(qs$Count, qs$Cycle))
means$Cycle<-rownames(means)
means$Cycle<-as.numeric(means$Cycle)
q25s <- by(qs, qs$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.25), simplify=TRUE)
q50s <- by(qs, qs$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.5), simplify=TRUE)
q75s <- by(qs, qs$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.75), simplify=TRUE)
cums <- by(qs, qs$Cycle, function(foo) sum(foo$Count), simplify=TRUE)
if(!all(sapply(list(names(q25s), names(q50s), names(q75s), names(cums)), identical, rownames(means)))) {
  stop("Calculated quantiles/means weren't compatible.")
}

Q25=as.vector(q25s)
Q75=as.vector(q75s)
Q50=as.vector(q50s)

sum_reads<-cbind(means, Q25,Q50,Q75)

qs_plot<-qs%>%ggplot(aes(x=Cycle, y=Score)) + 
  geom_tile(aes(fill=log2(Count),alpha=0.1))+ 
  scale_fill_gradient(low="white", high="orange")+
  geom_line(data=sum_reads,aes(Cycle,sum_reads[,1]))+theme(legend.position="none")+
  scale_x_continuous(breaks=c(0, 50, 100,150,200,250,300))
