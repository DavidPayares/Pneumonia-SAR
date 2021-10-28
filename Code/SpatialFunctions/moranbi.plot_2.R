assign("moranbi.plot",
  function(x, y, listw, spChk = NULL, labels = NULL, xlab = NULL, ylab = NULL,
           zero.policy = NULL, quiet = NULL, ...){ # zero.policy = NULL, quiet = NULL, ...)
  if (!inherits(listw, "listw"))
      stop(paste(deparse(substitute(listw)), "is not a listw object"))
  if (is.null(quiet))
      quiet <- !get("verbose", envir = .spdepOptions)
  stopifnot(is.vector(x))
  stopifnot(is.vector(y))
  stopifnot(is.logical(quiet))
  if (is.null(zero.policy))
      zero.policy <- get("zeroPolicy", envir = .spdepOptions)
  stopifnot(is.logical(zero.policy))
  xname <- deparse(substitute(x))
  yname <- deparse(substitute(y))
  if (!is.numeric(x))
    stop(paste(xname, "is not a numeric vector"))
  if (any(is.na(x)))
    stop("NA in X")
  if (!is.numeric(y))
    stop(paste(yname, "is not a numeric vector"))
  if (any(is.na(y)))
    stop("NA in Y")
  n <- length(listw$neighbours)
  if (n != length(x))
    stop("objects of different length")
  if (is.null(spChk))
    spChk <- get.spChkOption()
  if (spChk && !chkIDs(x, y, listw))
    stop("Check of data and weights ID integrity failed")
  labs <- TRUE
  if (is.logical(labels) && !labels)
    labs <- FALSE
  if (is.null(labels) || length(labels) != n)
    labels <- as.character(attr(listw, "region.id"))
  wy <- lag.listw(listw, y, zero.policy = zero.policy)
  if (is.null(xlab))
    xlab <- xname
  if (is.null(ylab))
    ylab <- paste("spatially lagged", yname)
  plot(x, wy, xlab = xlab, ylab = ylab, ...)
  if (zero.policy) {
  n0 <- wy == 0
  if (any(n0)) {
  symbols(x[n0], wy[n0], inches = FALSE, circles = rep(diff(range(x))/50,
                                                           length(which(n0))), bg = "grey", add = TRUE)
    }
  }
  xwy.lm <- lm(wy ~ x)
  abline(xwy.lm)
  abline(h = mean(wy), lty = 2)
  abline(v = mean(x), lty = 2)
  infl.xwy <- influence.measures(xwy.lm)
  is.inf <- which(apply(infl.xwy$is.inf, 1, any))
  points(x[is.inf], wy[is.inf], pch = 9, cex = 1.2)
  if (labs)
    text(x[is.inf], wy[is.inf], labels = labels[is.inf],
         pos = 2, cex = 0.7)
  rownames(infl.xwy$infmat) <- labels
  if (!quiet)
  summary(infl.xwy)
  invisible(infl.xwy)
})
