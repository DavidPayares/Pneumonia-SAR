##################################################
#                 Libraries                      #
##################################################
# Sys.setenv(R_INSTALL_STAGED = FALSE)
# options(warn = -1)


packages <- c('lmtest','RColorBrewer','classInt','spdep','TeachingDemos','shapefiles','sp','maptools',
             'scatterplot3d','geoR','spatial','fBasics','car','aplpack','RODBC','ggplot2','spgrass6',
             'adespatial','RANN','ade4','olsrr','rgeos','rgdal','spdep','spgwr','GWmodel','nnet','olsrr',
             'stats','classInt','gridExtra','lmtest','spatialreg','car','MASS','caret','glmnet')

for(p in packages){
  if(!require(p,character.only = TRUE)) install.packages(p, dependencies = TRUE)
  library(p,character.only = TRUE)
}



# Set working directory
dir <- '/mnt/d/M.Sc. Gesopatial Tecnologies/GeoMundus/GeoMundus 2019/Neumonia'
setwd(file.path(dir))

# Load some useful functions
funDir <- file.path(dir, 'SpatialFunctions')
spFun <- list.files(funDir, pattern = '.r', full.names = T)
invisible(lapply(spFun, source))

# Reset par for sequential plotting
resetPar <- function() { dev.new(); op <- par(no.readonly = TRUE); dev.off(); op}


##################################################
#                     Data                       #
##################################################

# Load Data
years <- list('04','07','11','14')
pneuData <- paste0(dir, '/SHPFinal/Neumonia', years, '.shp')

# read Data
pneuShp <- lapply(pneuData, st_read)
pneuNames <- as.character(paste0('pneu',years))
names(pneuShp) <- pneuNames

##################################################
#             Contiguity Matrices                #
##################################################

# Contiguity criteria

# Queen
queenW <- lapply(pneuShp, poly2nb)
# Root
rootW <- lapply(pneuShp, poly2nb, queen = F)

# Distance criteria

# Coords
coords <- lapply(pneuShp, function(x) {st_sfc(st_centroid(x$geometry))})

# Delanuay Triangulation
triangW <- lapply(coords, tri2nb)
# Sphere of Influence
sphereW <- lapply(1:length(pneuNames), function(x) {graph2nb(soi.graph(triangW[[x]], coords[[x]]))})
names(sphereW) <- pneuNames
# Gabriel
gabrielW <- lapply(1:length(pneuNames), function(x) {graph2nb(gabrielneigh(coords[[x]]), sym = T)})
names(gabrielW) <- pneuNames
# Relative Neighbors
relativeW <- lapply(1:length(pneuNames), function(x) {graph2nb(relativeneigh(coords[[x]]), sym = T)})
names(relativeW) <- pneuNames

# K-neighbors critreria

# Polygons IDs
IDs <- row.names(as.data.frame(pneuShp$pneu04))

# K-1 neighbor
kn1W <- lapply(1:length(pneuNames), function(x) {knn2nb(knearneigh(coords[[x]], k = 1), row.names = IDs)})
names(kn1W) <- pneuNames

# K-2 neighbors
kn2W <- lapply(1:length(pneuNames), function(x) {knn2nb(knearneigh(coords[[x]], k = 2), row.names = IDs)})
names(kn2W) <- pneuNames

# K-4 neighbors
kn4W <- lapply(1:length(pneuNames), function(x) {knn2nb(knearneigh(coords[[x]], k = 4), row.names = IDs)})
names(kn4W) <- pneuNames


# Plot weights matrices
op=par(mfrow=c(3,3), oma=c(0,0,0,0), mar=c(1,1,1,1))
matrices <- c('queen','root', 'triang', 'sphere', 'gabriel', 'relative', 'kn1', 'kn2', 'kn4')

for (w in matrices){
  plot(st_geometry(pneuShp$pneu04), border = 'gray')
  plot(get(paste0(w, 'W'))$pneu04, coords$pneu04, add=T,  pch=20, cex.main=1.5, frame.plot=FALSE , col="gray25")
  title(w)
}


## Spatial weights matrices PCNM

# Function to get AIC for each matrix
getMatrixAIC <- function(x, w){
  pcnm <- test.W(pneuShp[[x]]$SMR, w[[x]])
  aic <- pcnm$all$AICc
  return(aic)
}

# empty list
pcnm <- list()
# Get AICs
for (w in matrices){
  aic <- lapply(1:length(pneuNames), function(x){getMatrixAIC(x,get(paste0(w, 'W')))})
  pcnm[[w]] <- aic
}

# AICs per year and matrix
pcnm <- as.data.frame(do.call(rbind, pcnm))
names(pcnm) <- pneuNames
pcnm

# Best weights matrices per year (based on AIC)
bestMatrices <- lapply(1:length(pneuNames),function(x){rownames(pcnm[which.min(pcnm[[x]]),])})
names(bestMatrices) <- pneuNames # 2004: queen, 2007: kn2, 2011: kn4, 2014; kn2
bestMatrices

## Compute Moran's I to asses matrices

#set seed
set.seed(123)

# Moran's I test
getMoranPvalue <- function(x, w){
  moran <- moran.mc(pneuShp[[x]]$SMR, nb2listw(w[[x]]), nsim = 1000, zero.policy = TRUE)
  p <- moran$p.value
  return(p)
}

# empty list
moranp <- list()
# Get p values
for (w in matrices){
  pvalue <- lapply(1:length(pneuNames), function(x){getMoranPvalue(x,get(paste0(w, 'W')))})
  moranp[[w]] <- pvalue
}

# p-values per year and matrix
moranp <- as.data.frame(do.call(rbind, moranp))
names(moranp) <- pneuNames # 2004: queen, 2007: kn2, 2011: kn4, 2014: queen
moranp

# Select matrices based on both PCNM and Moran's I
spatialW <- c('queenW','kn2W', 'kn4W', 'queenW')


# Plot Moran's I
par(resetPar())
op=par(mfrow=c(2,2))
for (i in 1:length(spatialW)){
  moran.plot(pneuShp[[i]]$SMR, nb2listw(get(spatialW[[i]])[[i]]), xlab = 'SMR', ylab = 'Spatially lagged SMR', main = paste0('20', years[[i]], '\n ', spatialW[[i]] ))
}

# Plot Moran's I correlograms
par(resetPar())
op=par(mfrow=c(2,2))
for (i in 1:length(spatialW)){
  correlogram <- sp.correlogram(get(spatialW[[i]])[[i]], pneuShp[[i]]$SMR, order=5, method= "I",style="W", zero.policy=T)
  print(correlogram)
  plot(correlogram, main = paste0('20', years[[i]], '\n ', spatialW[[i]] ))
} # Besides the first lag, the third lag is significative for 2004 and 2011

# Define spatial lag variables
lagQueenW3 <- nblag(queenW$pneu04, 3)[[3]]
lagKn4W3 <- nblag(kn4W$pneu11, 3)[[3]]

# Add spatial lags to original data
pneuShp$pneu04$lag3SMR <- lag.listw(nb2listw(lagQueenW3, style = 'W'), pneuShp$pneu04$SMR)
pneuShp$pneu11$lag3SMR <- lag.listw(nb2listw(lagKn4W3, style = 'W'), pneuShp$pneu11$SMR)


##################################################
#             Exploratory Analysis               #
##################################################

## Global Analysis

# Moran's I
moranI <- lapply(1:length(pneuNames), function(x) {moran.test(pneuShp[[x]]$SMR, nb2listw(get(spatialW[[x]])[[x]]))})
names(moranI) <- pneuNames
moranI

# Getid's Ord
getisO <- lapply(1:length(pneuNames), function(x) {globalG.test(pneuShp[[x]]$SMR, nb2listw(get(spatialW[[x]])[[x]], style = 'B'))})
names(getisO) <- pneuNames
getisO

# Bivariate Moran's I

# Change column names (data integrity)
for (y in 1:length(pneuNames)){
  for(col in 77:91){
    colnames(pneuShp[[y]])[col] <-  sub(paste0(years[[y]] ,'.*'), "", colnames(pneuShp[[y]])[col])
  }
}

# Bivariate Moran's I (p-value and plot)
par(resetPar())
op=par(mfrow=c(4,4), mar=c(4,4,1,1),oma=c(1,1,1,1))
for (col in c(78:87, 89:91)){
  mi <- moranbi.test(pneuShp$pneu11$SMR, pneuShp$pneu11[[col]], nb2listw(get(spatialW[[3]])[[3]]), N= 999)
  moranbi.plot(pneuShp$pneu11$SMR, pneuShp$pneu11[[col]], nb2listw(get(spatialW[[3]])[[3]]), N= 999,graph=T, quiet = T, main = paste0('I = ', round(mi$Observed,3), ', p-value = ', mi$p.value), xlab = 'SMR', ylab = colnames(pneuShp$pneu11[col])[1], cex.main=0.9)
} # For 2014 (For other years the variables need to be adressed for each dataset, that is, pneu'year'.  e.g., pneu07)

# 2004: TEM, CPM, CVV, ESC.
# 2007: NUT, IPSE
# 2011: CPM, CVV, CVD, ESC, IPSE
# 2014: CPM, NUT, DEP, NBI, CVV, VAC

## Local Analysis

# LISA maps
invisible(lapply(1:length(pneuNames), function (x) {moran.cluster(pneuShp[[x]]$SMR, nb2listw(get(spatialW[[x]])[[x]]), zero.policy = T, pneuShp$pneu04$geometry, significant=T, main = paste0('20', years[[x]]))}))

# Geary's C
gearyC <- lapply(1:length(pneuNames), function(x) {geary.test(pneuShp[[x]]$SMR, nb2listw(get(spatialW[[x]])[[x]]), zero.policy=F)})
names(gearyC) <- pneuNames
gearyC

##################################################
#                OLS Regression                  #
#################################################

selectLinearModel <- function(dataset, variablesNames, formula){
  
  # Add spatial order lag for 2004, 2011
  if (substr(formula, nchar(formula)-6, nchar(formula)) == 'lag3SMR'){
    variablesNames <- append(variablesNames, 'lag3SMR')
  }
  
  # Data frame without geometry
  noGeoData <<- as.data.frame(dataset[variablesNames])
  noGeoData <<- noGeoData[,-ncol(noGeoData)]
  
  ### Variables criteria
  
  ### Stepwise Method
  m <- stats::step(lm(as.formula(formula), data = noGeoData), trace=0)
  # Get selected variables and compute regression
  formulaStep <- paste0(variablesNames[1], ' ~ ', paste(names(m$model)[-1], collapse = ' + '))
  steplm <- lm(as.formula(formulaStep) , data = noGeoData)
  
  ### Lasso method
  lasso <- cv.glmnet(x=as.matrix(noGeoData[,-1]), y=as.matrix(noGeoData$SMR), alpha=1, nfolds = nrow(noGeoData), type.measure="mse", family="gaussian", grouped = F)
  c <- coef(lasso, s = lasso$lambda.min)
  lassoPredictors <- row.names(c)[which(c!=0)][-1]
  # Get selected variables
  formulaLasso <- paste0(variablesNames[1], ' ~ ', paste(lassoPredictors, collapse = ' + '))
  lassolm <- lm(as.formula(formulaLasso) , data = noGeoData)
  
  ###  Elastic Net Method
  # Validation method
  train.control <- trainControl(method = "LOOCV")
  # Select optimal parameters (alpha and lambda)
  elastic_net <- train(as.formula(formula), data = noGeoData, method = "glmnet", trControl = train.control)
  alpha <- elastic_net$bestTune$alpha
  elastic <- cv.glmnet(x=as.matrix(noGeoData[,-1]), y=as.matrix(noGeoData$SMR), alpha= alpha, nfolds = nrow(noGeoData), type.measure="mse", family="gaussian", grouped = F)
  c <- coef(elastic, s = elastic$lambda.min)
  elasticPredictors <- row.names(c)[which(c!=0)][-1]
  # Get selected variables
  formulaElastic <- paste0(variablesNames[1], ' ~ ', paste(elasticPredictors, collapse = ' + '))
  elasticlm <- lm(as.formula(formulaElastic) , data = noGeoData)
  
  ### Cross validation
  
  # stepwise (AIC)
  stepcv <- train(as.formula(formulaStep), data = noGeoData, method = "lm", trControl = train.control)
  # Lasso
  lassocv <- train(as.formula(formulaLasso), data = noGeoData, method = "lm", trControl = train.control)
  # Elastic Net
  elasticcv <- train(as.formula(formulaElastic), data = noGeoData, method = "lm", trControl = train.control)
  
  RMSE <- c(AIC(stepcv$finalModel), AIC(lassocv$finalModel), AIC(elasticcv$finalModel))
  
  results <- cbind(RMSE)
  rownames(results) <- c('classic', 'lasso', 'elasticNet')
  
  bestCriteria <- rownames(results)[which.min(results)]
  
  if (bestCriteria == 'classic'){
    bestModel <- formulaStep
  } else if (bestCriteria == 'lasso'){
    bestModel <- formulaLasso
  } else{
    bestModel <- formulaElastic    
  }
  
  message(paste("The best criteria for selecting covariates is:", bestCriteria, "\nFormula selected model:\n", bestModel, '\n'))
  return(bestModel)
}

# Variables Names
variablesNames <- c('SMR','IDD','TEM','CPM','DEP','NUT','CVV','CVD','ESC','IPSE','VAC','NBI','ACU')
# Set linear model
baseFormula <- paste0(variablesNames[1], ' ~ ', paste(variablesNames[2:length(variablesNames)], collapse = ' + '))
# Add higher order spatial lags (2004, 2011)
formulas <- c(paste0(baseFormula, ' + lag3SMR'), baseFormula, paste0(baseFormula, ' + lag3SMR'), baseFormula )

# ols Models
olsModels <- lapply(1:length(pneuNames), function (x){ cat(paste0('***20', years[[x]], '***\n')); selectLinearModel(dataset = pneuShp[[x]], variablesNames = variablesNames , formula = formulas[[x]])})

# Final OLS models
names(olsModels) <- pneuNames
olslm <- lapply(1:length(pneuNames), function (x){lm(as.formula(olsModels[[x]]), data = pneuShp[[x]])})
names(olslm) <- pneuNames

# Residuals analysis

# Spatial autocorrelation in residuals
lapply(1:length(pneuNames), function (x){
  lm.morantest(olslm[[x]], nb2listw(get(spatialW[[x]])[[x]]))
})   # only 2011 has residual spatial autocorrelation

#Plot moran's I
par(resetPar())
par(mfrow = c(2,2))
invisible(lapply(1:length(pneuNames), function (x){
  moran.plot(olslm[[x]]$residuals, nb2listw(get(spatialW[[x]])[[x]]), xlab = 'SMR', ylab = 'spatially lagged residuals' , main = paste0('20', years[[x]]))
}))


# OLS regression assumptions
par(resetPar())
par(mfrow=c(2,2))
invisible(lapply(1:length(pneuNames), function (x){
  lmyear <- olslm[[x]];
  cat(paste0('\n***** 20', years[[x]], ' *****\n'))
  # Linearity:
  plot(lmyear, which=c(1), main = paste0('20', years[[x]]))
  # Residuals normality: Shapiro Wilk test
  normality <- shapiro.test(residuals(lmyear))
  if (normality$p.value >= 0.05) {cat('Normal\n')} else {cat('non-normal\n')}
  # Heteroscedasticity: Breusch Pagan test
  homos <- bptest(lmyear)
  if (homos$p.value >= 0.05) {cat('Homocedastic\n')} else {cat('Heteroscedastic\n')}
  # Specificity: Ramsey's RESET test
  specif <- resettest(lmyear)
  if (specif$p.value >= 0.05) {cat('Good specified\n')} else {cat('Poorly specified\n')}
  # Multicolineality: VIF
  vif <- vif(lmyear)
  if (length(vif[vif > 10]) > 0) {cat('Multicolinearity\n')} else {cat('Non-multicolinearity\n')}
}))  # 2007 suffers of non-normality

# Correct non-normality 2007

# log transformation in dependent variables
logFormula07 <- gsub("DEP", "log(DEP)", olsModels$pneu07)
logFormula07 <- gsub("ESC", "log(ESC)", logFormula07)
logFormula07 <- gsub("IPSE", "log(IPSE)", logFormula07)

olsModels$pneu07 <- logFormula07
olslm$pneu07 <- lm(as.formula(logFormula07), data = pneuShp$pneu07)

# Normality test
shapiro.test(residuals(olslm$pneu07)) #normal


##################################################
#                  SAR Models                    #
##################################################

## Langrange Multipliers
lagrange <- lapply(1:length(pneuNames), function(x) {lagMul <- summary(lm.LMtests(lm(as.formula(olsModels[[x]]), data = pneuShp[[x]]), nb2listw(get(spatialW[[x]])[[x]]), test = 'all')); print(lagMul);bestML <- names(lagMul[which.min(lagMul$results$p.value)]) ; return(bestML)}) 
lagrange # Lagrange Multipliers suggests a Lag model for every year

## Spatial autorregresive models
sarModels <- function(lm, data, listw){
  
  # remove spatial lag for models with Wx
  cov <- unlist(strsplit(lm, split = " ~ "))[2]
  if (substr(cov, nchar(cov)-6, nchar(cov)) == 'lag3SMR'){
    wx <- paste0(' ~ ', substr(cov,1, nchar(cov)-10))
  }else{
    wx <- paste0(' ~ ', cov)
  }
  
  lm <- as.formula(lm)
  
  lagSar <- lagsarlm(lm, data ,listw)                                   # Spatial lag model (WY)
  errSar <- errorsarlm(lm, data, listw)                                 # Error model (We)
  durSar <- lagsarlm(lm, data ,listw, Durbin = as.formula(wx))          # Durbin model (WY. WX)
  sacSar <- sacsarlm(lm, data ,listw, type="sac")                       # SARAR model (WY, We)
  slxSar <- lmSLX(lm, data, listw, Durbin = as.formula(wx))             # SLX model (WX)
  durErrSar <- errorsarlm(lm, data, listw, Durbin = as.formula(wx))     # SDEM model (WX, We)
  gnSar <-  sacsarlm(lm, data ,listw, Durbin = as.formula(wx))          # General Nesting (WY,WX,We)
  ols <- lm(lm, data)                                                   # OLS regression
  
  sarModels <- list(lagSar, errSar, slxSar, durSar, sacSar, durErrSar, gnSar, ols) # All models
  
  statLag <- lapply(sarModels, function(x){val <- cbind(aic = AIC(x), loglik = logLik(x)[[1]]); return(val)})
  statLag <- as.data.frame(do.call(rbind, statLag), row.names = c('lag','error','slx', 'durbin','sac','erdur','gns','ols'))
  return(statLag)
}

# Test models
sarResults <- lapply(1:length(pneuNames), function (x) {sarModels(olsModels[[x]], pneuShp[[x]], nb2listw(get(spatialW[[x]])[[x]]))})
names(sarResults) <- pneuNames
sarResults # GNS are overfitted models (see. https://onlinelibrary.wiley.com/doi/10.1111/jors.12188)
# 2004: durbin model, 2007: lag model, 2011: slx model, 2014: error model

# Get best models (based on AIC and Log Likekihood)
bestSar <- lapply(1:length(pneuNames), function(x) {rownames(sarResults[[x]])[c(1:6,8)][which.min(sarResults[[x]]$aic[c(1:5,8)])]}) # No GNS or sDEM (Overfitting)
bestSar

# Selected Models
erd04 <- errorsarlm(as.formula(olsModels$pneu04), pneuShp$pneu04, nb2listw(queenW$pneu04), Durbin = ~ IDD + DEP + CVV + NBI) 
ols07 <- lm(as.formula(olsModels$pneu07), pneuShp$pneu07)
lag11 <- lagsarlm(as.formula(olsModels$pneu11), pneuShp$pneu11, nb2listw(kn4W$pneu11)) 
erd14 <- errorsarlm(as.formula(olsModels$pneu14), pneuShp$pneu14, nb2listw(queenW$pneu14), Durbin = ~ TEM + CPM + NUT + CVV + IPSE + VAC + NBI)
sarReg <- list(erd04, ols07, lag11, erd14)


LR.sarlm(lm(as.formula(olsModels[[1]]), data = pneuShp[[1]]), erd04)


# Models summary
lapply(sarReg, function(x) {summary(x, Nagelkerke=T, Hausman = T)})

# SAR models assumptions

# Residual plot function for SarLm objects
residual.plot <- function(model, year) {
  plot(residuals(model) ~ model$fitted.values, xlab = "fitted values", ylab = "residuals", main = year)
  abline(h=0, lty="dotted")
  lines(lowess(model$fitted.values, residuals(model)), col="red")
}

# Model check
par(resetPar())
par(mfrow=c(2,2))
invisible(lapply(1:length(pneuNames), function (x){
  saryear <- sarReg[[x]];
  cat(paste0('\n***** 20', years[[x]], ' *****\n'))
  # Linearity:
  residual.plot(saryear, paste0('20', years[[x]]))
  # Residuals normality: Shapiro Wilk test
  normality <- shapiro.test(residuals(saryear))
  if (normality$p.value >= 0.05) {cat('Normal\n')} else {cat('non-normal\n')}
  # Heteroscedasticity: Breusch Pagan test
  if (class(saryear)[1] %in% c('SlX','lm')) {homos <- bptest(saryear)} else { homos <- bptest.Sarlm(saryear)}
  if (homos$p.value >= 0.05) {cat('Homocedastic\n')} else {cat('Heteroscedastic\n')}
  # autocorrelation
  moran <- moran.test(residuals(saryear), nb2listw(get(spatialW[[x]])[[x]]))
  if (moran$p.value >= 0.05) {cat('Non-spatially autocorrelated\n')} else {cat('Spatially autocorrelated\n')}
})) 



##################################################
#                   Inference                    #
##################################################

# Impacts
impactModels <- lapply(c(1,3:4), function (x) {summary(impacts(sarReg[[x]], listw = nb2listw(get(spatialW[[x]])[[x]]), R = 1000), zstats = T, short=TRUE, density = T)}) # 2014: error models doesn't have impacts
names(impactModels) <- pneuNames[c(1,3:4)]
impactModels


# Save variables
save.image("~/Pneumonia/Pnuemonia/Peunomia.RData")
