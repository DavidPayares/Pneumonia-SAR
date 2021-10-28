assign("moran.bi",
  function(X,Y,listw,zero.policy = NULL,adjust.n = TRUE){
    if (!inherits(listw, "listw")) 
        stop(paste(deparse(substitute(listw)), "is not a listw object"))
    if (!is.numeric(X)) 
        stop(paste(deparse(substitute(X)), "is not a numeric vector"))
    if (!is.numeric(Y)) 
        stop(paste(deparse(substitute(Y)), "is not a numeric vector"))  
    if (is.null(zero.policy)) 
        zero.policy <- get("zeroPolicy", envir = .spdepOptions)
    stopifnot(is.logical(zero.policy))  
   wc <- spweights.constants(listw, zero.policy = zero.policy, adjust.n = adjust.n)
   n <- wc$n
   morans<-(n/(wc$S0))%*%((t(scale(X))%*%as.matrix(as_dgRMatrix_listw(listw))%*%scale(Y))/(t(scale(X))%*%scale(X)))
   return(as.vector(morans))
})
