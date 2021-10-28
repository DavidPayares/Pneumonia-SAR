assign("moran.bi1",
  function(x,y,listw,zero.policy = NULL,adjust.n = TRUE, NAOK=FALSE,...){
    if (!inherits(listw, "listw")) 
        stop(paste(deparse(substitute(listw)), "is not a listw object"))
    if (!is.numeric(x)) 
        stop(paste(deparse(substitute(x)), "is not a numeric vector"))
    if (!is.numeric(y)) 
        stop(paste(deparse(substitute(y)), "is not a numeric vector"))  
    if (is.null(zero.policy)) 
        zero.policy <- get("zeroPolicy", envir = .spdepOptions)
    stopifnot(is.logical(zero.policy))  
   wc <- spweights.constants(listw, zero.policy = zero.policy, adjust.n = adjust.n)
   n <- wc$n
   morans<-(n/(wc$S0))%*%((t(scale(x))%*%as.matrix(as_dgRMatrix_listw(listw))%*%scale(y))/(n-1))
 	xx <- mean(x, na.rm=NAOK)
	zx <- x - xx
	yy <- mean(y, na.rm=NAOK)
	zy <- y - yy
	zz <- sum((zx^2)*(zy^2), na.rm=NAOK)
	K <- (length(x)*sum(zz, na.rm=NAOK))/((sum(zx^2, na.rm=NAOK))*sum(zy^2, na.rm=NAOK))
 	res <- list(I=as.vector(morans), K=K)
	res
})
