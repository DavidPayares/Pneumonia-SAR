assign("moranbi1.test",
  function(x, y, listw, randomisation=TRUE, zero.policy=NULL,
	alternative="greater", rank = FALSE, na.action=na.fail, spChk=NULL, 
	adjust.n=TRUE) {
	alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
	if (!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if (!is.numeric(x)) stop(paste(deparse(substitute(x)),
		"is not a numeric vector"))
	if (!is.numeric(y)) stop(paste(deparse(substitute(y)),
		"is not a numeric vector"))   
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
	if (is.null(spChk)) spChk <- get.spChkOption()
	if (spChk && !chkIDs(x, listw))
		stop("Check of data and weights ID integrity failed")
 	if (spChk && !chkIDs(y, listw))
		stop("Check of data and weights ID integrity failed")

#	if (any(is.na(x))) stop("NA in X")
	xname <- deparse(substitute(x))
 	yname <- deparse(substitute(y))
	wname <- deparse(substitute(listw))
	NAOK <- deparse(substitute(na.action)) == "na.pass"
	x <- na.action(x)
	xna.act <- attr(x, "na.action")
	if (!is.null(xna.act)) {
	    subset <- !(1:length(listw$neighbours) %in% xna.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
	}
 	y <- na.action(y)
	yna.act <- attr(y, "na.action")
	if (!is.null(yna.act)) {
	    subset <- !(1:length(listw$neighbours) %in% yna.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
	}

	n <- length(listw$neighbours)
	if (n != length(x)) stop("objects of different length")
	if (n != length(y)) stop("objects of different length")	
	wc <- spweights.constants(listw, zero.policy=zero.policy, 
		adjust.n=adjust.n)
#	S02 <- wc$S0*wc$S0
	res <- moran.bi1(x, y, listw, wc$n, wc$S0, zero.policy=zero.policy,NAOK=NAOK)
	I <- res$I
	K <- res$K
	if (rank) K <- (3*(3*wc$n^2 -7))/(5*(wc$n^2 - 1))
	
	W <- as.matrix(as_dgRMatrix_listw(listw))   # wc$n=49
	S0 <- t(rep(1,n))%*%(W)%*%rep(1,n)
	S02 <- S0^2
	S3 <- t(rep(1,n))%*%(W*t(W))%*%rep(1,n)
	S4 <- t(rep(1,n))%*%(W*W)%*%rep(1,n)
	S5 <- t(rep(1,n))%*%(W%*%W)%*%rep(1,n)
	S6 <- t(rep(1,n))%*%(t(W)%*%W+W%*%t(W))%*%rep(1,n)
	S1 <- S3+S4 # t(rep(1,wc$n))%*%(W*W+W*t(W))%*%rep(1,wc$n)
	S2 <- 2*S5+S6
	rxy <- cor(x,y)
#	K <- sum((as.numeric(scale(x,center=TRUE, scale=F))^2)*as.numeric(scale(y,center=TRUE, scale=F))^2)/(((var(x)*wc$n1)/wc$n)*((var(y)*wc$n1)/wc$n))
	EI <- -rxy/(wc$n1)
#	EI2 <- ((rxy^2)*wc$n*(2*(S02-S2+S1)+(2*S3-2*S5)*wc$n3+S3*wc$n2*wc$n3)-K*(6*(S02-S2+S1)+(4*S1-2*S2)*wc$n3+S1*wc$n2*wc$n3)+
#	          wc$n*((S02-S2+S1)+(2*S4-S6)*wc$n3+S4*wc$n2*wc$n3))/(wc$n1*wc$n2*wc$n3*S02)
#	VI <- EI2-EI^2
#	Z(Ixy) <- (I-EI)/sqrt(VI)
	
	if(randomisation) {
		VI <- ((rxy^2)*wc$n*(2*(S02-S2+S1)+(2*S3-2*S5)*wc$n3+S3*wc$n2*wc$n3)+wc$n*((S02-S2+S1)+(2*S4-S6)*wc$n3+S4*wc$n2*wc$n3))
		tmp <- K*(6*(S02-S2+S1)+(4*S1-2*S2)*wc$n3+S1*wc$n2*wc$n3)
                if (tmp > VI) warning("Kurtosis overflow,\ndistribution of variable does not meet test assumptions")
		VI <- (VI - tmp) / (wc$n1*wc$n2*wc$n3*S02)
                tmp <- (VI - EI^2)
                if (tmp < 0) warning("Negative variance,\ndistribution of variable does not meet test assumptions")
		VI <- tmp
	} else {
		VI <- (wc$nn*wc$S1 - wc$n*wc$S2 + 3*S02) / (S02*(wc$nn - 1))
                tmp <- (VI - EI^2)
                if (tmp < 0) warning("Negative variance,\ndistribution of variable does not meet test assumptions")
		VI <- tmp
	}
	ZI <- (I - EI) / sqrt(VI)
	statistic <- ZI
	names(statistic) <- "Bivariate Moran Z(I) statistic"
        if (alternative == "two.sided") 
		PrI <- 2 * pnorm(abs(ZI), lower.tail=FALSE)
        else if (alternative == "greater")
            PrI <- pnorm(ZI, lower.tail=FALSE)
        else PrI <- pnorm(ZI)
	if (!is.finite(PrI) || PrI < 0 || PrI > 1) 
		warning("Out-of-range p-value: reconsider test arguments")
	vec <- c(I, EI, VI)
	names(vec) <- c("Bivariate Moran I statistic", "Expectation", "Variance")
	method <- paste("Bivariate Moran I test under", ifelse(randomisation,
	    "randomisation", "normality"))
	data.name <- paste(xname, ifelse(rank,
		"using rank correction",""), "\nweights:",
		wname, ifelse(is.null(xna.act), "", paste("\nomitted:", 
	    paste(xna.act, collapse=", "))), ifelse(is.null(yna.act), "", paste("\nomitted:", 
	    paste(yna.act, collapse=", "))),"\n")
	res <- list(statistic=statistic, p.value=PrI, estimate=vec, 
	    alternative=alternative, method=method, data.name=data.name)
	if (!is.null(xna.act)) attr(res, "na.action") <- xna.act
 	if (!is.null(yna.act)) attr(res, "na.action") <- yna.act
	class(res) <- "htest"
	res
})

