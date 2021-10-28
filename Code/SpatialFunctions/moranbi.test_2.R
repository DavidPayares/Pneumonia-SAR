assign("moranbi.test",
  function(X,Y,listw,zero.policy=NULL,adjust.n=TRUE,N,graph=F,print.results=T,...){
  observed<-moran.bi(X,Y,listw=listw,zero.policy=zero.policy,adjust.n = adjust.n,...)
  DF <- data.frame(1:length(X),X,Y)
  names(DF) <- c("Obs","X","Y")
  if(length(X)<8){
    require(combinat)
    X1<-unique(permn(DF$Obs))
  }
  else{
    X1<-randomize_vector(DF$Obs,N)
  }
  if(length(Y)<8){
    require(combinat)
    Y1<-unique(permn(DF$Obs))
  }
  else{
    Y1<-randomize_vector(DF$Obs,N)
  }
  rxy <- cor(X,Y)
  store<-rep(NA,length(X1))
  for(i in 1:length(store)){
    store[i]<-moran.bi(X[X1[[i]]],Y[Y1[[i]]],listw,zero.policy = zero.policy, adjust.n = adjust.n)
  }
  if(observed>=0){
    p.value<-(sum(ifelse(store>observed,1,0))+1)/(length(store)+1)
    expected <- -rxy/(length(X)-1)
  }
  else if(observed<0){
    p.value<-(sum(ifelse(store<observed,1,0))+1)/(length(store)+1)
    expected <- -rxy/(length(X)-1)
  }
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