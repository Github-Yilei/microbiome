plotErrors <- function(dq, nti=c("A","C","G","T"), ntj=c("A","C","G","T"), obs=TRUE, err_out=TRUE, err_in=FALSE, nominalQ=FALSE) {
  ACGT <- c("A", "C", "G", "T")
  if(!(all(nti %in% ACGT) && all(ntj %in% ACGT)) || any(duplicated(nti)) || any(duplicated(ntj))) {
    stop("nti and ntj must be nucleotide(s): A/C/G/T.")
  }
  
  dq <- err
  
  if(!is.null(dq$trans)) {
    if(ncol(dq$trans) <= 1) {
      stop("plotErrors only supported when using quality scores in the error model (i.e. USE_QUALS=TRUE).")
    }
    transdf<-dq$trans
    transdf<-data.frame(transdf)
    transdf$Var1<-rownames(transdf)
    transdf<-pivot_longer(transdf, starts_with("X"), values_to="values", names_to="Var2")
    transdf$Var2<-str_replace(transdf$Var2, "X","")
    transdf$Var1<-factor(transdf$Var1)
    transdf$Var2<-as.integer(transdf$Var2)
    transdf<-data.frame(transdf)
    colnames(transdf) <- c("Transition", "Qual", "count")
  } else if(!is.null(dq$err_out)) {
    if(ncol(dq$err_out) <= 1) {
      stop("plotErrors only supported when using quality scores in the error model (i.e. USE_QUALS=TRUE).")
    }
    transdf<-dqF$err_out
    transdf<-data.frame(transdf)
    transdf$Var1<-rownames(transdf)
    transdf<-pivot_longer(transdf, starts_with("X"), values_to="values", names_to="Var2")
    transdf$Var2<-str_replace(transdf$Var2, "X","")
    transdf$Var1<-factor(transdf$Var1)
    transdf$Var2<-as.integer(transdf$Var2)
    transdf<-data.frame(transdf)
    colnames(transdf) <- c("Transition", "Qual", "est")
  } else {
    stop("Non-null observed and/or estimated error rates (dq$trans or dq$err_out) must be provided.")
  }
  transdf$from <- substr(transdf$Transition, 1, 1)
  transdf$to <- substr(transdf$Transition, 3, 3)
  
  if(!is.null(dq$trans)) {
    tot.count <- tapply(transdf$count, list(transdf$from, transdf$Qual), sum)
    transdf$tot <- mapply(function(x,y) tot.count[x,y], transdf$from, as.character(transdf$Qual))
    transdf$Observed <- transdf$count/transdf$tot
  } else {
    transdf$Observed <- NA
    obs <- FALSE
  }
  if(!is.null(dq$err_out)) {
    transdf$Estimated <- mapply(function(x,y) dq$err_out[x,y], transdf$Transition, as.character(transdf$Qual))
  } else {
    transdf$Estimated <- NA
    err_out <- FALSE
  }
  if(!is.null(dq$err_in)) {
    # If selfConsist, then err_in is a list. Use the first err_in, the initial error rates provided.
    ei <- dq$err_in; if(is.list(ei)) ei <- ei[[1]]
    transdf$Input <- mapply(function(x,y) ei[x,y], transdf$Transition, as.character(transdf$Qual))
  } else {
    transdf$Input <- NA
    err_in <- FALSE
  }
  transdf$Nominal <- (1/3)*10^-(transdf$Qual/10)
  transdf$Nominal[transdf$Transition %in% c("A2A", "C2C", "G2G", "T2T")] <- 1 - 10^-(transdf$Qual[transdf$Transition %in% c("A2A", "C2C", "G2G", "T2T")]/10)
  transition_data<<-transdf
  p <- ggplot(data=transdf[transdf$from %in% nti & transdf$to %in% ntj,], aes(x=Qual))
  if(obs) p <- p + geom_point(aes(y=Observed), color="gray40", na.rm=TRUE)
  if(err_in)   p <- p + geom_line(aes(y=Input), linetype="dashed")
  if(err_out)  p <- p + geom_line(aes(y=Estimated))
  if(nominalQ) p <- p + geom_line(aes(y=Nominal), color="red")
  p <- p + scale_y_log10()
  p <- p + facet_wrap(~Transition, nrow=length(nti))
  p <- p + xlab("Consensus quality score") + ylab("Error frequency (log10)")
  p <- p + theme_bw()
  p

}

