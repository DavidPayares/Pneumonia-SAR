assign("randomize_vector",
   function(X,N){
   lst<-list()
   for(i in 1:N){
     lst[[i]]<-sample(X,length(X))
   }
   return(lst)
})