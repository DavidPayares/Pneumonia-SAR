assign("moranbi.testA",
  function(X,Y,W,N,test=c("positive","negative","two-sided"),graph=F,print.results=T){
  observed<-moran.bi(X,Y,W)
  DF <- data.frame(1:length(X),X,Y)
  names(DF) <- c("Obs","X","Y")
  if(length(X)<8){
    require(combinat)
    X1<-unique(permn(DF$Obs))
  }
  else{
    X1<-randomize_vector(DF$Obs,N)
  }
  store<-rep(NA,length(X1))
  for(i in 1:length(store)){
    store[i]<-moran.bi(X[X1[[i]]],Y[X1[[i]]],W)
  }

  if(test=="NULL")
  if(observed>=0){
    p.value<-(sum(ifelse(store>observed,1,0))+1)/(length(store)+1)
    expected=(-1/(length(X)-1))
  }
  else if(observed<0){
    p.value<-(sum(ifelse(store<observed,1,0))+1)/(length(store)+1)
    expected=(-1/(length(X)-1))
  }


#  if(length(test)>1){test=test[1]}
#  if(test=="positive"){
#    p.value<-(sum(ifelse(store>observed,1,0))+1)/(length(store)+1)
#    expected=(-1/(length(X)-1))
#  }
#  else if(test=="negative"){
#    p.value<-(sum(ifelse(store<observed,1,0))+1)/(length(store)+1)
#    expected=(-1/(length(X)-1))
#  }
#  else if(test=="two-sided"){
#    store<-abs(store-mean(store))
#    observed<-observed-mean(store)
#    expected=abs(-1/(length(X)-1))
#    p.value<-(sum(ifelse(store>observed,1,0))+1)/(length(store)+1)
#  }
#  else{stop("test must be at 'positive', 'negative' or 'two-sided'")}
  if(graph==T){
    require(ggplot2)
    tmp.dat<-data.frame(store=store,observed=observed)
    export.graph<-ggplot(tmp.dat,aes(x=store))+
      scale_y_sqrt()+geom_density()+
      geom_vline(aes(xintercept=observed),color="red",size=1)+
      xlab("Bivariate Moran's I Coefficient") +
      ylab("Empirical Density")+theme_bw()
    print(export.graph)
  }
  if(print.results==T){
    print(list(Observed=observed,Expected=expected,p.value=p.value))
  }
  z<-list(Observed=observed,Expected=expected,p.value=p.value,Values=store)
})