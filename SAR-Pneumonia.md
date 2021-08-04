---
title: "Penumonia SAR Models"
output:
  html_document:
    theme: readable
    highlight: tango
    keep_md: yes
    toc: true
    toc_float:
      collapsed: true
      smooth_scroll: false
    number_sections: true
---

This is a reproducible example of the paper [Pneumonia Sar](). In this research, we perform a spatial autorregresive model analysis to measure the spatial influence of socio-economic, enviromental, behavioural and healthcare factors in pneumonia mortality rates. 

# Getting Started

## Packages

First, we load the packages needed for processing and analyzing our data.


```r
# Omitting warmings (be carefull)
options(warn = -1)

# Packages
packages <- c('lmtest','RColorBrewer','classInt','spdep','TeachingDemos','shapefiles','sp','maptools',
             'scatterplot3d','geoR','spatial','fBasics','car','aplpack','RODBC','ggplot2','spgrass6',
             'adespatial','RANN','ade4','olsrr','rgeos','rgdal','spdep','spgwr','GWmodel','nnet','olsrr',
             'stats','classInt','gridExtra','lmtest','car','MASS','caret','glmnet')
```


For loading and/or installing the packages, you can run the following command. The packages that are installed in your machine be automatically loaded. Those that are not installed be installed and loaded to your R session. Remember to provide permits if linux is used as operative system.


```r
# Permits
#Sys.setenv(R_INSTALL_STAGED = FALSE)

# Load packages
for(p in packages){
  if(!require(p,character.only = TRUE)) install.packages(p, dependencies = TRUE)
  library(p,character.only = TRUE)
}
```


## Working directory

It is easier to work in R using a working directory. That way you ensure your inputs and outputs be loaded and stored in a global directory.


```r
# Set working directory
dir <- '/mnt/d/M.Sc. Gesopatial Tecnologies/GeoMundus/GeoMundus 2019/Neumonia'
setwd(file.path(dir))
```


We also use some useful functions in the exploratory analysis. Load them all from the *SpatialFunctions* folder. Remember this folder has to be in your current working directory.


```r
# Load some useful functions
funDir <- file.path(dir, 'SpatialFunctions')
spFun <- list.files(funDir, pattern = '.r', full.names = T)
invisible(lapply(spFun, source))
```


And we add an extra function for resting the plot options everytime we plot.


```r
# Reset par for sequential plotting
resetPar <- function() { dev.new(); op <- par(no.readonly = TRUE); dev.off(); op}
```


# Data

Our data contains the [standardized mortality ratio (SMR)](https://ibis.health.state.nm.us/resource/SMR_ISR.html#:~:text=Standardized%20Mortality%20Ratio%20(SMR)%20is,the%20same%20age%2Fsex%20groups.) of Pneumonia (our study variable) in 19 districts of Bogotá, Colombia for the years 2004, 2007, 2011 and 2014. The data also contains socio-economic, enviromental and healthcare covariates.

As our data is stored in a *shapefile (.shp)* format, we load it like so.


```r
# Load Data
years <- list('04','07','11','14')
pneuData <- paste0(dir, '/SHPFinal/Neumonia', years, '.shp')
```

And then we read it using the `sp` package.


```r
# read Data
pneuShp <- lapply(pneuData, st_read)
pneuNames <- as.character(paste0('pneu',years))
names(pneuShp) <- pneuNames
```

# Exploratory Analysis

## Spatial weights matrix

One of the main aspects of any spatial autorregresive model is the contiguity matrix, also known as the spatial weights matrix ($\mathbf{W}$). This matrix encodes the spatial dependence and influence of one region with its neighbors. There are many ways to define $\mathbf{W}$. Usually, an expert proposes a potential spatial weights matrix based on their knowledge of the phenomenon (study variable). For our research, we consider most of the parametric matrices configurations as we do not make any assumptions about the underlying spatial structure of the SMR variable. 


```r
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
```


Now we plot the matrices to check the relationships among regions.


```r
# Plot weights matrices
op=par(mfrow=c(3,3), oma=c(0,0,0,0), mar=c(1,1,1,1))
matrices <- c('queen','root', 'triang', 'sphere', 'gabriel', 'relative', 'kn1', 'kn2', 'kn4')

for (w in matrices){
  plot(st_geometry(pneuShp$pneu04), border = 'gray')
  plot(get(paste0(w, 'W'))$pneu04, coords$pneu04, add=T,  pch=20, cex.main=1.5 , col="gray25")
  title(w)
}
```

![](SAR-Pneumonia_files/figure-html/unnamed-chunk-9-1.png)<!-- -->


### Principal coordinates of neighbour matrices (PCNM) {-}

Now we evaluate the best matrix for our data. We use the [Principal coordinates of neighbour matrices (PCNM)](https://www.sciencedirect.com/science/article/abs/pii/S0304380006000925) as our selection criteria. For further details, please refer to the PCNM method.

First we create a function to extract the AIC of the PCNM method to compare and extract the best matrix per each year.


```r
## Spatial weights matrices PCNM

# Function to get AIC for each matrix
getMatrixAIC <- function(x, w){
  pcnm <- test.W(pneuShp[[x]]$SMR, w[[x]])
  aic <- pcnm$all$AICc
  return(aic)
}
```
 

And then we get the AIC values.


```r
# empty list
pcnm <- list()
# Get AICs
for (w in matrices){
  aic <- lapply(1:length(pneuNames), function(x){getMatrixAIC(x,get(paste0(w, 'W')))})
  pcnm[[w]] <- aic
}
```


These are the best matrices per year based on the AIC value.


```r
# AICs per year and matrix
pcnm <- as.data.frame(do.call(rbind, pcnm))
names(pcnm) <- pneuNames
pcnm
```

```
##             pneu04    pneu07    pneu11    pneu14
## queen    -68.77277 -58.97563 -48.16067 -31.85952
## root     -61.15346 -61.13476  -48.8688 -25.39557
## triang   -61.19358 -66.05321 -62.92984 -23.53263
## sphere   -53.70998  -71.8989 -55.78013 -24.97818
## gabriel  -55.85301 -63.07068 -57.84765 -27.43804
## relative -55.71152 -71.29338 -59.65928 -35.22539
## kn1      -59.24907 -66.80696 -56.95576  -20.2536
## kn2      -67.48367 -83.38051 -58.48591  -42.9033
## kn4      -53.45341  -80.2701 -70.07799 -29.22466
```

```r
# Best weights matrices per year (based on AIC)
bestMatrices <- lapply(1:length(pneuNames),function(x){rownames(pcnm[which.min(pcnm[[x]]),])})
names(bestMatrices) <- pneuNames # 2004: queen, 2007: kn2, 2011: kn4, 2014; kn2
bestMatrices
```

```
## $pneu04
## [1] "queen"
## 
## $pneu07
## [1] "kn2"
## 
## $pneu11
## [1] "kn4"
## 
## $pneu14
## [1] "kn2"
```

### Moran's Index {-}

We confirm the PCNM finding using the Moran's Index of our study variable and the spatial weights matrices. We expect to find spatial autocorrelation between the SMR and $\mathbf{W}$. We use the Moran's I statistically significant $p-value$ to assess the matrices. 


```r
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
```

```
##              pneu04     pneu07      pneu11      pneu14
## queen    0.03396603  0.1188811  0.04995005  0.04395604
## root     0.01798202  0.1658342  0.02297702  0.01798202
## triang    0.1048951  0.1258741 0.000999001 0.000999001
## sphere   0.02597403  0.1228771 0.007992008  0.02997003
## gabriel  0.06393606  0.1468531 0.000999001 0.008991009
## relative 0.02097902  0.1378621  0.01598402   0.1198801
## kn1      0.08591409 0.01398601  0.00999001   0.1868132
## kn2      0.06193806 0.07392607  0.00999001   0.2127872
## kn4      0.05294705 0.09090909 0.000999001  0.01998002
```


And now having the matrices selected by both the PCNM method and the Moran's I, we chose our final $\mathbf{W}$'s for our study period. We selected the matrices whose PCNM's *AIC* were the lowest and Moran's I $p-values$ were statistically significant.


```r
# Selectec matrices based on both PCNM and Moran's I
spatialW <- c('queenW','kn2W', 'kn4W', 'queenW')
```


Now, let us plot the Moran's I of our selected matrices.


```r
# Plot Moran's I
par(resetPar())
op=par(mfrow=c(2,2))
for (i in 1:length(spatialW)){
  moran.plot(pneuShp[[i]]$SMR, nb2listw(get(spatialW[[i]])[[i]]), xlab = 'SMR', ylab = 'Spatially lagged SMR', main = paste0('20', years[[i]], '\n ', spatialW[[i]] ))
}
```

![](SAR-Pneumonia_files/figure-html/unnamed-chunk-15-1.png)<!-- -->


### Higher order matrices {-}

We also want to know if our matrices  have statistically significant spatial lags or greater orders, that is, finding autocorrelation throughout longer distances or neighbors orders. If we find statistically significant orders of the matrices, we have to include them in the spatial autorregresive models to guarantee a proper inclusion of the spatial dependence among the regions.


```r
# Plot Moran's I correlograms
par(resetPar())
op=par(mfrow=c(2,2))
for (i in 1:length(spatialW)){
  correlogram <- sp.correlogram(get(spatialW[[i]])[[i]], pneuShp[[i]]$SMR, order=5, method= "I",style="W", zero.policy=T)
  print(correlogram)
  plot(correlogram, main = paste0('20', years[[i]], '\n ', spatialW[[i]] ))
} # Besides the first lag, the third lag is significative for 2004 and 2011
```

![](SAR-Pneumonia_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

We see that in 2004 and 2011 the 3rd order matrices are statistically significant. We add the lagged variable as covariates to include them in the spatial autorregresive (SAR) models.


```r
# Define spatial lag variables
lagQueenW3 <- nblag(queenW$pneu04, 3)[[3]]
lagKn4W3 <- nblag(kn4W$pneu11, 3)[[3]]

# Add spatial lags to original data
pneuShp$pneu04$lag3SMR <- lag.listw(nb2listw(lagQueenW3), pneuShp$pneu04$SMR)
pneuShp$pneu11$lag3SMR <- lag.listw(nb2listw(lagKn4W3), pneuShp$pneu11$SMR)
```

## Global Analysis

As an essential condition of the SAR models, the dependent variable (SMR), independent variables (covariates) or distrubances (error term) must exhibit spatial autocorrelation. This spatial autocorrelation can be present across study area (Global) or in specific regions (Local). The former can be found using both the Moran's I and Geary's C statistics. 

### Global Moran's I {-}

We assess the presence of global spatial autocorrelation in our dependent variable using the global Moran's I.


```r
# Moran's I
moranI <- lapply(1:length(pneuNames), function(x) {moran.test(pneuShp[[x]]$SMR, nb2listw(get(spatialW[[x]])[[x]]))})
names(moranI) <- pneuNames
moranI
```

```
## $pneu04
## 
## 	Moran I test under randomisation
## 
## data:  pneuShp[[x]]$SMR  
## weights: nb2listw(get(spatialW[[x]])[[x]])    
## 
## Moran I statistic standard deviate = 1.9898, p-value = 0.02331
## alternative hypothesis: greater
## sample estimates:
## Moran I statistic       Expectation          Variance 
##        0.23413855       -0.05555556        0.02119646 
## 
## 
## $pneu07
## 
## 	Moran I test under randomisation
## 
## data:  pneuShp[[x]]$SMR  
## weights: nb2listw(get(spatialW[[x]])[[x]])    
## 
## Moran I statistic standard deviate = 1.5402, p-value = 0.06176
## alternative hypothesis: greater
## sample estimates:
## Moran I statistic       Expectation          Variance 
##        0.22110039       -0.05555556        0.03226555 
## 
## 
## $pneu11
## 
## 	Moran I test under randomisation
## 
## data:  pneuShp[[x]]$SMR  
## weights: nb2listw(get(spatialW[[x]])[[x]])    
## 
## Moran I statistic standard deviate = 3.8994, p-value = 4.821e-05
## alternative hypothesis: greater
## sample estimates:
## Moran I statistic       Expectation          Variance 
##        0.42364583       -0.05555556        0.01510198 
## 
## 
## $pneu14
## 
## 	Moran I test under randomisation
## 
## data:  pneuShp[[x]]$SMR  
## weights: nb2listw(get(spatialW[[x]])[[x]])    
## 
## Moran I statistic standard deviate = 1.8748, p-value = 0.03041
## alternative hypothesis: greater
## sample estimates:
## Moran I statistic       Expectation          Variance 
##        0.19958513       -0.05555556        0.01851956
```

As we can see, all years except 2007 exhibit global spatial autocorrelation using a significance level of $p-value$ < 0.05.

### Geary's C {-}

Geary's C is an attempt to determine if adjacent observations of the same phenomenon are correlated locally. We use this statistics to check for local clusters in our study area. The results are very similar to those made by the Global Moran's I. 2007 does not exhibit spatial autocorrelation. 



```r
gearyC <- lapply(1:length(pneuNames), function(x) {geary.test(pneuShp[[x]]$SMR, nb2listw(get(spatialW[[x]])[[x]]), zero.policy=F)})
names(gearyC) <- pneuNames
gearyC
```

```
## $pneu04
## 
## 	Geary C test under randomisation
## 
## data:  pneuShp[[x]]$SMR 
## weights: nb2listw(get(spatialW[[x]])[[x]]) 
## 
## Geary C statistic standard deviate = 1.9454, p-value = 0.02586
## alternative hypothesis: Expectation greater than statistic
## sample estimates:
## Geary C statistic       Expectation          Variance 
##        0.71801475        1.00000000        0.02101025 
## 
## 
## $pneu07
## 
## 	Geary C test under randomisation
## 
## data:  pneuShp[[x]]$SMR 
## weights: nb2listw(get(spatialW[[x]])[[x]]) 
## 
## Geary C statistic standard deviate = 0.50977, p-value = 0.3051
## alternative hypothesis: Expectation greater than statistic
## sample estimates:
## Geary C statistic       Expectation          Variance 
##        0.89449902        1.00000000        0.04283139 
## 
## 
## $pneu11
## 
## 	Geary C test under randomisation
## 
## data:  pneuShp[[x]]$SMR 
## weights: nb2listw(get(spatialW[[x]])[[x]]) 
## 
## Geary C statistic standard deviate = 3.2585, p-value = 0.0005601
## alternative hypothesis: Expectation greater than statistic
## sample estimates:
## Geary C statistic       Expectation          Variance 
##        0.55868473        1.00000000        0.01834302 
## 
## 
## $pneu14
## 
## 	Geary C test under randomisation
## 
## data:  pneuShp[[x]]$SMR 
## weights: nb2listw(get(spatialW[[x]])[[x]]) 
## 
## Geary C statistic standard deviate = 2.3687, p-value = 0.008926
## alternative hypothesis: Expectation greater than statistic
## sample estimates:
## Geary C statistic       Expectation          Variance 
##        0.63791091        1.00000000        0.02336808
```


### Bivariate Moran's I {-}

To find spatial association between the dependent variable and the covariates, we compute the bivariate Moran's Index. Covariates with statistically significant degree of spatial correlation with the SMR are initial potential candidates for the spatial autorregresive models. We calculate the Moran's I $p-value$ and plot to scrutinized the covariates that spatially correlate with the SMR. 

The analysis is perform for the entier study period. However, for the sake of brevity, we display only the results for 2014.


```r
# Change column names (data integrity)
for (y in 1:length(pneuNames)){
  for(col in 77:91){
    colnames(pneuShp[[y]])[col] <-  sub(paste0(years[[y]] ,'.*'), "", colnames(pneuShp[[y]])[col])
  }
}

# Bivariate Moran's I (p-value and plot)
par(resetPar())
op=par(mfrow=c(4,4), mar=c(4,4,1,1),oma=c(1,1,1,1))
for (col in c(78:88, 89:90)){
  mi <- moranbi.test(pneuShp$pneu14$SMR, pneuShp$pneu14[[col]], nb2listw(get(spatialW[[4]])[[4]]), N= 999)
  moranbi.plot(pneuShp$pneu14$SMR, pneuShp$pneu14[[col]], nb2listw(get(spatialW[[4]])[[4]]), N= 999,graph=T, quiet = T, main = paste0('I = ', round(mi$Observed,3), ', p-value = ', mi$p.value), xlab = 'SMR', ylab = colnames(pneuShp$pneu14[col])[1], cex.main=0.9)
} # For 2014 (For other years the variables need to be adressed for each dataset, that is, pneu'year'.  e.g., pneu07)


# 2004: TEM, CPM, CVV, ESC.
# 2007: NUT, IPSE
# 2011: CPM, CVV, CVD, ESC, IPSE
# 2014: CPM, NUT, DEP, NBI, CVV, VAC
```

![](SAR-Pneumonia_files/figure-html/unnamed-chunk-20-1.png)<!-- -->

CPM, NUT, DEP, NBI and CVV show bivariate spatial autocorrelation with SMR. We consider these variables in the SAR modeling as they can contribute to explaining our dependent variable.

## Local Analysis

As we see in the global analysis, such statistics (Moran’s I, Geary's C) are designed to identify global spatial autocorrelation (clustering). Such clustering is a characteristic of the complete spatial pattern and does not provide an indication of the location of the clusters.

Local analysis allows to assess the spatial correlation of each location with its local neighboorhood. For this analysis, we use the Local Spatial Autocorrelation via the Moran's I and the Getis Ord statistic.

The global Moran's I informed us of local spatial autocorrelation in the study period. However, it does not determine the existence of clusters of high/low values.

### LISA maps (Moran's I) {-}

The Local Moran statistic identifies local clusters and local spatial outliers. Local Moran's I allows for a classification of the significant locations as High-High and Low-Low spatial clusters, and High-Low and Low-High spatial outliers.

We see that our data shows local clustering in the center of the city for almost every year except in 2014 wher eonly one region is statistically significant.


```r
invisible(lapply(1:length(pneuNames), function (x) {moran.cluster(pneuShp[[x]]$SMR, nb2listw(get(spatialW[[x]])[[x]]), zero.policy = T, pneuShp$pneu04$geometry, significant=T, main = paste0('20', years[[x]]))}))
```

![](SAR-Pneumonia_files/figure-html/unnamed-chunk-21-1.png)<!-- -->![](SAR-Pneumonia_files/figure-html/unnamed-chunk-21-2.png)<!-- -->![](SAR-Pneumonia_files/figure-html/unnamed-chunk-21-3.png)<!-- -->![](SAR-Pneumonia_files/figure-html/unnamed-chunk-21-4.png)<!-- -->


### Getis Ord {-}

We corroborate the LISA results using the local statistic Getis Ord statistic.

Every year, except 2014, presents positive and statistically significant spatial clustering, that is, there are high values concentration in our study area. The results for 2014 usign Getis Ord are similar to those obtained throught the LISA Moran's I statistic, in which, spatial clustering was not found.


```r
# Getid's Ord
getisO <- lapply(1:length(pneuNames), function(x) {globalG.test(pneuShp[[x]]$SMR, nb2listw(get(spatialW[[x]])[[x]], style = 'B'))})
names(getisO) <- pneuNames
getisO
```

```
## $pneu04
## 
## 	Getis-Ord global G statistic
## 
## data:  pneuShp[[x]]$SMR 
## weights: nb2listw(get(spatialW[[x]])[[x]], style = "B") 
## 
## standard deviate = 2.4654, p-value = 0.006843
## alternative hypothesis: greater
## sample estimates:
## Global G statistic        Expectation           Variance 
##       2.690322e-01       2.456140e-01       9.022585e-05 
## 
## 
## $pneu07
## 
## 	Getis-Ord global G statistic
## 
## data:  pneuShp[[x]]$SMR 
## weights: nb2listw(get(spatialW[[x]])[[x]], style = "B") 
## 
## standard deviate = 2.3359, p-value = 0.009749
## alternative hypothesis: greater
## sample estimates:
## Global G statistic        Expectation           Variance 
##       0.1191726891       0.1111111111       0.0000119107 
## 
## 
## $pneu11
## 
## 	Getis-Ord global G statistic
## 
## data:  pneuShp[[x]]$SMR 
## weights: nb2listw(get(spatialW[[x]])[[x]], style = "B") 
## 
## standard deviate = 2.6555, p-value = 0.003959
## alternative hypothesis: greater
## sample estimates:
## Global G statistic        Expectation           Variance 
##       2.480737e-01       2.222222e-01       9.476995e-05 
## 
## 
## $pneu14
## 
## 	Getis-Ord global G statistic
## 
## data:  pneuShp[[x]]$SMR 
## weights: nb2listw(get(spatialW[[x]])[[x]], style = "B") 
## 
## standard deviate = 0.14643, p-value = 0.4418
## alternative hypothesis: greater
## sample estimates:
## Global G statistic        Expectation           Variance 
##       0.2484907062       0.2456140351       0.0003859437
```
# Spatial Autorregresive Models

Albeit the exploratory analysis insinuates the existence of spatial autocorrelation in our data, it does not inform about the level of influence of the independent variables or the underlying spatial dependence structure. The former can be addressed by performing a linear regression in which the independent variables' influence over the dependet variable can be assesed.

## Linear regression model

We perform three different models to find the best combination of covariates that explain the SMR. We not assume any dependence (except the ones found in the bivariate analysis) between the independet variables and the outcome variable (SMR). We select the classic Ordinary least squeare (OLS regression) and two regularized variations, the Lasso and the Elastic Net regression. The formers reduce the risk of overfitting, the variance, the correlation effect between variables and the influence of irrelevant predictors derived in a normal OLS regression. 

First, we create a function that estimates, assesses and selects the best linear model among the three different regressions. For this, we use a leave one out cross-validation approach and we retain the model with the lowest Root-Mean-Square Error (RMSE) and AIC.


```r
selectLinearModel <- function(dataset, variablesNames, formula){
  
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
  
  RMSE <- c(stepcv$results$RMSE, lassocv$results$RMSE, elasticcv$results$RMSE)
  
  results <- cbind(RMSE)
  rownames(results) <- c('stepAIC', 'lasso', 'elasticNet')
  
  bestCriteria <- rownames(results)[which.min(results)]
  
  if (bestCriteria == 'stepAIC'){
    bestModel <- formulaStep
  } else if (bestCriteria == 'lasso'){
    bestModel <- formulaLasso
  } else{
    bestModel <- formulaElastic    
  }
  
  message(paste("The best criteria for selecting covariates is:", bestCriteria, "\nFormula selected model:\n", bestModel, '\n'))
  return(bestModel)
}
```

We define the covariates and the linear formula.


```r
# Variables Names
variablesNames <- c('SMR','IDD','TEM','CPM','DEP','NUT','CVV','CVD','ESC','IPSE','VAC','NBI','ACU')
# Set linear model
formula <- paste0(variablesNames[1], ' ~ ', paste(variablesNames[2:length(variablesNames)], collapse = ' + '))
```

And now we use the function we just created ```selectLinearModel```. For 2004, the Lasso and Elastic Net approach could not reach convergence, so we selected the OLS model as the best model. For 2007, 2011 and 2014 the best models were selected by the Lasso approach.


```r
#ols Model (2004)
model04 <- step(lm(as.formula(formula), data = pneuShp$pneu04), trace = 0)
# Get selected variables
formula2004 <- paste0(variablesNames[1], ' ~ ', paste(names(model04$model)[-1], collapse = ' + '))
cat('***2004***\n')
```

```
## ***2004***
```

```r
message(paste0('The best criteria for selecting covariates is: stepAIC 
Formula selected model:\n', formula2004,'\n'))
```

```
## The best criteria for selecting covariates is: stepAIC 
## Formula selected model:
## SMR ~ TEM + NUT + CVD + ESC + VAC
```

```r
# ols Models for 2007, 2011, 2014 (method failed with 2004)
olsModels <- lapply(2:length(pneuNames), function (x){ cat(paste0('***20', years[[x]], '***\n')); selectLinearModel(dataset = pneuShp[[x]], variablesNames = variablesNames , formula = formula)})
```

```
## ***2007***
```

```
## The best criteria for selecting covariates is: lasso 
## Formula selected model:
##  SMR ~ DEP + ESC + IPSE
```

```
## ***2011***
```

```
## The best criteria for selecting covariates is: stepAIC 
## Formula selected model:
##  SMR ~ CPM + CVV + IPSE + VAC + ACU
```

```
## ***2014***
```

```
## The best criteria for selecting covariates is: lasso 
## Formula selected model:
##  SMR ~ IDD + TEM + CPM + NUT + CVV + IPSE + VAC
```
The outcome linear models contain potential covariates that could explain the SMR. We save the models for further analysis.


```r
# Final OLS models
olsModels <- c(formula2004, olsModels)
names(olsModels) <- pneuNames
olslm <- lapply(1:length(pneuNames), function (x){lm(as.formula(olsModels[[x]]), data = pneuShp[[x]])})
names(olslm) <- pneuNames
```

### Residuals Analysis {-}

Autocorrelation in the model residuals, regardless of its nature, violates the OLS assumptions. If spatial autocorrelation is identified, a model that accounts for spatial dependence has to used instead.

We assess the linear models' residuals to find residual spatial autocorrelation and to confirm the use of spatial autorregresive models for our data.


```r
# Spatial autocorrelation in residuals
lapply(1:length(pneuNames), function (x){
  lm.morantest(olslm[[x]], nb2listw(get(spatialW[[x]])[[x]]))
})   # only 2011 has residual spatial autocorrelation
```

```
## [[1]]
## 
## 	Global Moran I for regression residuals
## 
## data:  
## model: lm(formula = as.formula(olsModels[[x]]), data = pneuShp[[x]])
## weights: nb2listw(get(spatialW[[x]])[[x]])
## 
## Moran I statistic standard deviate = 1.2402, p-value = 0.1074
## alternative hypothesis: greater
## sample estimates:
## Observed Moran I      Expectation         Variance 
##       0.02231240      -0.12940269       0.01496444 
## 
## 
## [[2]]
## 
## 	Global Moran I for regression residuals
## 
## data:  
## model: lm(formula = as.formula(olsModels[[x]]), data = pneuShp[[x]])
## weights: nb2listw(get(spatialW[[x]])[[x]])
## 
## Moran I statistic standard deviate = 0.11585, p-value = 0.4539
## alternative hypothesis: greater
## sample estimates:
## Observed Moran I      Expectation         Variance 
##      -0.07056865      -0.09025602       0.02888145 
## 
## 
## [[3]]
## 
## 	Global Moran I for regression residuals
## 
## data:  
## model: lm(formula = as.formula(olsModels[[x]]), data = pneuShp[[x]])
## weights: nb2listw(get(spatialW[[x]])[[x]])
## 
## Moran I statistic standard deviate = 2.6324, p-value = 0.00424
## alternative hypothesis: greater
## sample estimates:
## Observed Moran I      Expectation         Variance 
##       0.19844785      -0.09649925       0.01255430 
## 
## 
## [[4]]
## 
## 	Global Moran I for regression residuals
## 
## data:  
## model: lm(formula = as.formula(olsModels[[x]]), data = pneuShp[[x]])
## weights: nb2listw(get(spatialW[[x]])[[x]])
## 
## Moran I statistic standard deviate = -0.42552, p-value = 0.6648
## alternative hypothesis: greater
## sample estimates:
## Observed Moran I      Expectation         Variance 
##      -0.21759796      -0.17145190       0.01176081
```
The Moran's I test in the linear regression residuals concluded that only 2011 present residual spatial autocorrelation. Although 2004, 2007 and 2014 show no evidence of residual spatial autocorrelaton, we detected spatial dependence in the outcome variable. The results from the Moran's I in both the dependent variable and the OLS residuals suggest initial spatial autorregresive model configurations that considers the nature of this dependence. For instance, in 2011 we might need a model that accounts for spatial structure in the error term.

Now let us plot the Moran's I to visualize the results.


```r
#Plot moran's I
par(resetPar())
par(mfrow = c(2,2))
invisible(lapply(1:length(pneuNames), function (x){
  moran.plot(olslm[[x]]$residuals, nb2listw(get(spatialW[[x]])[[x]]), xlab = 'SMR', ylab = 'spatially lagged residuals' , main = paste0('20', years[[x]]))
}))
```

![](SAR-Pneumonia_files/figure-html/unnamed-chunk-28-1.png)<!-- -->

### Linear regression assumptions {-}

We apply also formal diagnostic tests to inspect model assumptions such as linearity, multicollinearity, and homoscedasticity. The explanatory variables that survive the linear regression and its assumptions will be considered potential candidates for a spatial regression model. Normality and spatial dependence in the residuals are also tested.


```r
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
```

```
## 
## ***** 2004 *****
```

```
## Normal
## Homocedastic
## Good specified
## Non-multicolinearity
## 
## ***** 2007 *****
```

```
## non-normal
## Homocedastic
## Good specified
## Non-multicolinearity
## 
## ***** 2011 *****
```

```
## Normal
## Homocedastic
## Good specified
## Non-multicolinearity
## 
## ***** 2014 *****
```

![](SAR-Pneumonia_files/figure-html/unnamed-chunk-29-1.png)<!-- -->

```
## Normal
## Homocedastic
## Good specified
## Non-multicolinearity
```

2007 suffers from non-normality. For tackling this violation of the lienar regression assumptions, we perform a logaritmic transformation to the independent covariates such as our residual will become normal distributed.


```r
# Correct non-normality 2007

# log transformation in dependent variables
logFormula07 <- gsub("DEP", "log(DEP)", olsModels$pneu07)
logFormula07 <- gsub("ESC", "log(ESC)", logFormula07)
logFormula07 <- gsub("IPSE", "log(IPSE)", logFormula07)

olsModels$pneu07 <- logFormula07
olslm$pneu07 <- lm(as.formula(logFormula07), data = pneuShp$pneu07)

# Normality test
shapiro.test(residuals(olslm$pneu07)) #normal
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  residuals(olslm$pneu07)
## W = 0.91001, p-value = 0.07405
```
## Spatial autorregresive models

Now that we count with the lagged variables (see [Higher order matrices](#Higher-order-matrices)), the covariates associated in space with the SMR (see [Bivariate Moran's I](Bivariate-Moran's-I)) and the ones found by the linear model, we can define the SAR models independent variables.

First, we add the lagged variables. Then we add the spatially autocorrelated covariates, if not present already in the linear models.


```r
# Redefine formulas to include spatial lags (2004 and 2011)
olsModels$pneu04 <- paste0(olsModels$pneu04, ' + lag3SMR')
olsModels$pneu11 <- paste0(olsModels$pneu11, ' + lag3SMR')
```


As we identified spatial autocorrelation, a model that accounts for spatial dependence has to be used. Spatial econometrics literature have developed models to incorporate the spatial dependence in different forms: (i) a spatially lagged dependent variable, (ii) spatially lagged independent variables, and (iii) a spatial structure in the error term. The simultaneous interaction of the three forms produce a spatial autorregresive model known as the General Nesting Spatial (GNS) Model. The GNS is expressed as,

$$\mathbf{y}=\rho \mathbf{W}_{\mathbf{y}}+\mathbf{X} \boldsymbol{\beta}+\mathbf{W} \mathbf{X} \theta+\varepsilon \\
\varepsilon=\lambda \mathbf{W} \varepsilon+\mathbf{u} \\
\mathbf{u}\sim \mathcal{N}(0,\,\sigma^{2})\ $$

where $\mathbf{y}$ represents a vector consisting of one observation on the dependent variable for every spatial unit, $\mathbf{X}$ the matrix of independent variables, $\mathbf{W}$ is the spatial weights matrix that describes the structure of dependence between units, $\mathbf{W}_{\mathbf{y}}$ denotes spatially lagged dependent variable, $\mathbf{W}_{\mathbf{X}}$ the spatially lagged independent variable, and $\mathbf{W}\varepsilon$ the spatial interaction effects in the error term. The scalar parameters $\rho$ and $\lambda$ measure the strength of dependence between units, while $\theta$, like $\boldsymbol{\beta}$, is a vector of response parameters. $\mathbf{u}$ is a vector of independently and identically distributed disturbance terms with zero mean and variance $\sigma$.

Other spatial autoregressive models can be obtained by restricting the GNS model spatial interactions, that is, omitting a form of spatial dependence. For example, the spatial lag model is a particular specification in which only the endogenous interactions are considered (spatially lagged dependent variable \mathbf{W}_{\mathbf{y}}).

\begin{equation}
\mathbf{y}=\rho \mathbf{W}_{\mathbf{y}}+\mathbf{X} \boldsymbol{\beta} + \varepsilon \\
\varepsilon\sim \mathcal{N}(0,\,\sigma^{2})\
\end{equation}

We employ the Lagrange Multiplier test to identify the appropriate spatial autoregressive models and their form or forms of spatial dependence.

The Lagrange Multiplier tests suggest a spatial error model for 2004, and a spatial lag model for the remaining years.


```r
## Langrange Multipliers
lagrange <- lapply(1:length(pneuNames), function(x) {lagMul <- summary(lm.LMtests(lm(as.formula(olsModels[[x]]), data = pneuShp[[x]]), nb2listw(get(spatialW[[x]])[[x]]), test = 'all')); bestML <- names(lagMul[which.min(lagMul$results$p.value)]) ; return(bestML)}) 
lagrange # Lagrange Multipliers suggests a Lag model for every year
```

```
## [[1]]
## [1] "RLMlag"
## 
## [[2]]
## [1] "LMlag"
## 
## [[3]]
## [1] "LMlag"
## 
## [[4]]
## [1] "LMlag"
```

However, Lagrange Multiplier test ignores models with exogenous spatial interaction (spatially lagged independent variables). We support its results with the Akaike information criterion (AIC) for the seven possible spatial models fitted using the Maximum likelihood estimation.

First, we define a function to build all the SAR models. For further details about the models refer to [The SLX Model (Halleck et al.)](https://onlinelibrary.wiley.com/doi/10.1111/jors.12188).


```r
## Spatial autorregresive models
sarModels <- function(lm, data, listw){
  lagSar <- lagsarlm(lm, data ,listw)               # Spatial lag model (WY)
  errSar <- errorsarlm(lm, data, listw)             # Error model (We)
  durSar <- lagsarlm(lm, data ,listw, Durbin = T)   # Durbin model (WY. WX)
  sacSar <- sacsarlm(lm, data ,listw)               # SARAR model (WY, We)
  slxSar <- lmSLX(lm, data, listw)                  # SLX model (WX)
  gnSar <-  sacsarlm(lm, data ,listw, Durbin = T)   # General Nesting (WY,WX,We)
  ols <- lm(lm, data)                               # OLS regression
  
  sarModels <- list(lagSar, errSar, durSar, sacSar, slxSar, gnSar, ols) # All models
  
  statLag <- lapply(sarModels, function(x){val <- cbind(aic = AIC(x), loglik = logLik(x)[[1]]); return(val)})
  statLag <- as.data.frame(do.call(rbind, statLag), row.names = c('lag','error','durbin','sac','slx','gns','ols'))
  return(statLag)
}
```


And then we test our models.


```r
# Test models
sarResults <- lapply(1:length(pneuNames), function (x) {sarModels(as.formula(olsModels[[x]]), pneuShp[[x]], nb2listw(get(spatialW[[x]])[[x]]))})
names(sarResults) <- pneuNames
sarResults
```

```
## $pneu04
##               aic   loglik
## lag     -2.914170 10.45709
## error   -2.896077 10.44804
## durbin -16.352729 23.17636
## sac     -7.097018 13.54851
## slx    -12.710649 20.35532
## gns    -30.911939 31.45597
## ols     -4.268996 10.13450
## 
## $pneu07
##              aic    loglik
## lag    -6.527506  9.263753
## error  -6.512847  9.256423
## durbin -1.155354  9.577677
## sac    -4.897860  9.448930
## slx    -2.901140  9.450570
## gns    -1.620479 10.810240
## ols    -8.436453  9.218227
## 
## $pneu11
##              aic     loglik
## lag    12.986909  2.5065457
## error  13.444602  2.2776992
## durbin  7.698594 11.1507029
## sac    13.808901  3.0955494
## slx     7.996370 10.0018150
## gns     4.195890 13.9020549
## ols    16.520817 -0.2604087
## 
## $pneu14
##              aic    loglik
## lag     4.101220  7.949390
## error   2.416914  8.791543
## durbin 11.487987 11.256006
## sac     3.518699  9.240650
## slx    13.836719  9.081641
## gns    -4.427679 20.213839
## ols     5.888563  6.055719
```

Let us select the models with the lowest AIC and highest Log-likelihood. These will be our final SAR models. We will ommit from the selection the OLS model, as there is no spatial interactions in it, and the GNS which has been found to be highly overfitted.


```r
# Get best models (based on AIC and Log Likekihood)
bestSar <- lapply(1:length(pneuNames), function(x) {rownames(sarResults[[x]])[which.min(sarResults[[x]]$aic[1:5])]}) # No GNS or OLS
bestSar
```

```
## [[1]]
## [1] "durbin"
## 
## [[2]]
## [1] "lag"
## 
## [[3]]
## [1] "durbin"
## 
## [[4]]
## [1] "error"
```

And now we compute our definitive model for each year.


```r
# Selected Models
dur04 <- lagsarlm(as.formula(olsModels$pneu04), pneuShp$pneu04, nb2listw(queenW$pneu04), Durbin = ~ CVD + NUT + VAC + TEM + ESC) 
lag07 <- lagsarlm(as.formula(olsModels$pneu07), pneuShp$pneu07 ,nb2listw(kn2W$pneu07))
dur11 <- lagsarlm(as.formula(olsModels$pneu11), pneuShp$pneu11, nb2listw(kn4W$pneu11), Durbin = ~  CPM + CVV + IPSE + VAC) 
err14 <- errorsarlm(as.formula(olsModels$pneu14), pneuShp$pneu14, nb2listw(queenW$pneu14))
sarReg <- list(dur04, lag07, dur11, err14)

# Summary
lapply(sarReg, function(x) {summary(x, Nagelkerke=T)})
```

```
## [[1]]
## 
## Call:lagsarlm(formula = as.formula(olsModels$pneu04), data = pneuShp$pneu04, 
##     listw = nb2listw(queenW$pneu04), Durbin = ~CVD + NUT + VAC + 
##         TEM + ESC)
## 
## Residuals:
##       Min        1Q    Median        3Q       Max 
## -0.177894 -0.049231  0.010979  0.048104  0.101514 
## 
## Type: mixed 
## Coefficients: (asymptotic standard errors) 
##                Estimate  Std. Error z value  Pr(>|z|)
## (Intercept) -20.6375101   3.8824309 -5.3156 1.063e-07
## TEM           0.0704103   0.0195309  3.6051 0.0003121
## NUT           0.0445005   0.0113838  3.9091 9.264e-05
## CVD           0.0538048   0.0131305  4.0977 4.173e-05
## ESC          -0.0525523   0.0093131 -5.6428 1.673e-08
## VAC           0.0077153   0.0013242  5.8263 5.667e-09
## lag3SMR      -1.0143229   0.1826959 -5.5520 2.825e-08
## lag.CVD       0.1191404   0.0439447  2.7111 0.0067052
## lag.NUT       0.0998653   0.0356886  2.7982 0.0051382
## lag.VAC       0.0248744   0.0042936  5.7934 6.899e-09
## lag.TEM       0.2911598   0.0434159  6.7063 1.996e-11
## lag.ESC       0.0169494   0.0354422  0.4782 0.6324891
## 
## Rho: -0.72793, LR test value: 6.3377, p-value: 0.01182
## Asymptotic standard error: 0.21967
##     z-value: -3.3138, p-value: 0.0009204
## Wald statistic: 10.981, p-value: 0.0009204
## 
## Log likelihood: 23.17634 for mixed model
## ML residual variance (sigma squared): 0.0045685, (sigma: 0.067591)
## Nagelkerke pseudo-R-squared: 0.93488 
## Number of observations: 19 
## Number of parameters estimated: 14 
## AIC: -18.353, (AIC for lm: -14.015)
## LM test for residual autocorrelation
## test value: 3.5095, p-value: 0.061019
## 
## 
## [[2]]
## 
## Call:lagsarlm(formula = as.formula(olsModels$pneu07), data = pneuShp$pneu07, 
##     listw = nb2listw(kn2W$pneu07))
## 
## Residuals:
##       Min        1Q    Median        3Q       Max 
## -0.320406 -0.092983  0.016820  0.129424  0.180856 
## 
## Type: lag 
## Coefficients: (asymptotic standard errors) 
##              Estimate Std. Error z value  Pr(>|z|)
## (Intercept) 10.021197   2.567442  3.9032 9.494e-05
## log(DEP)     0.420959   0.211853  1.9870   0.04692
## log(ESC)    -2.991608   0.709487 -4.2166 2.480e-05
## log(IPSE)   -0.028294   0.033675 -0.8402   0.40079
## 
## Rho: -0.062638, LR test value: 0.091053, p-value: 0.76284
## Asymptotic standard error: 0.23116
##     z-value: -0.27098, p-value: 0.78641
## Wald statistic: 0.073428, p-value: 0.78641
## 
## Log likelihood: 9.263753 for lag model
## ML residual variance (sigma squared): 0.022059, (sigma: 0.14852)
## Nagelkerke pseudo-R-squared: 0.53041 
## Number of observations: 19 
## Number of parameters estimated: 6 
## AIC: -6.5275, (AIC for lm: -8.4365)
## LM test for residual autocorrelation
## test value: 0.0099139, p-value: 0.92069
## 
## 
## [[3]]
## 
## Call:lagsarlm(formula = as.formula(olsModels$pneu11), data = pneuShp$pneu11, 
##     listw = nb2listw(kn4W$pneu11), Durbin = ~CPM + CVV + IPSE +         VAC)
## 
## Residuals:
##       Min        1Q    Median        3Q       Max 
## -0.301052 -0.054457  0.006944  0.069554  0.262058 
## 
## Type: mixed 
## Coefficients: (asymptotic standard errors) 
##                Estimate  Std. Error z value  Pr(>|z|)
## (Intercept) -4.0567e+01  1.3106e+01 -3.0953 0.0019661
## CPM         -1.9567e-02  4.9148e-03 -3.9813 6.853e-05
## CVV         -6.9515e-01  1.5512e-01 -4.4815 7.413e-06
## IPSE         1.3613e-03  9.0062e-04  1.5115 0.1306666
## VAC          4.9714e-02  1.3180e-02  3.7719 0.0001620
## ACU          5.3614e-01  1.4688e-01  3.6503 0.0002619
## lag3SMR      1.1838e+00  3.3809e-01  3.5014 0.0004627
## lag.CPM     -2.5836e-02  9.1506e-03 -2.8234 0.0047511
## lag.CVV     -1.0087e+00  2.4765e-01 -4.0730 4.642e-05
## lag.IPSE    -9.2824e-03  2.5826e-03 -3.5942 0.0003254
## lag.VAC      5.8532e-02  3.1943e-02  1.8324 0.0668992
## 
## Rho: 0.70102, LR test value: 6.0355, p-value: 0.014021
## Asymptotic standard error: 0.15546
##     z-value: 4.5093, p-value: 6.5027e-06
## Wald statistic: 20.334, p-value: 6.5027e-06
## 
## Log likelihood: 10.491 for mixed model
## ML residual variance (sigma squared): 0.017425, (sigma: 0.13201)
## Nagelkerke pseudo-R-squared: 0.87678 
## Number of observations: 19 
## Number of parameters estimated: 13 
## AIC: 5.018, (AIC for lm: 9.0535)
## LM test for residual autocorrelation
## test value: 0.70649, p-value: 0.40061
## 
## 
## [[4]]
## 
## Call:errorsarlm(formula = as.formula(olsModels$pneu14), data = pneuShp$pneu14, 
##     listw = nb2listw(queenW$pneu14))
## 
## Residuals:
##        Min         1Q     Median         3Q        Max 
## -0.2217362 -0.1008818 -0.0071495  0.0946907  0.2418540 
## 
## Type: error 
## Coefficients: (asymptotic standard errors) 
##                Estimate  Std. Error  z value  Pr(>|z|)
## (Intercept)  1.5428e+01  1.1188e+00  13.7897 < 2.2e-16
## IDD          4.9093e-03  2.6274e-02   0.1868  0.851781
## TEM          1.4335e-01  4.5656e-02   3.1398  0.001691
## CPM         -9.9509e-03  1.7192e-03  -5.7881 7.117e-09
## NUT         -1.0078e-01  2.5222e-02  -3.9955 6.455e-05
## CVV         -1.0633e+00  9.5253e-02 -11.1629 < 2.2e-16
## IPSE         4.1316e-04  6.0324e-05   6.8491 7.433e-12
## VAC         -2.5497e-02  3.9153e-03  -6.5121 7.410e-11
## 
## Lambda: -1.0938, LR test value: 5.4716, p-value: 0.019327
## Asymptotic standard error: 0.25247
##     z-value: -4.3326, p-value: 1.4739e-05
## Wald statistic: 18.771, p-value: 1.4739e-05
## 
## Log likelihood: 8.791543 for error model
## ML residual variance (sigma squared): 0.017827, (sigma: 0.13352)
## Nagelkerke pseudo-R-squared: 0.94196 
## Number of observations: 19 
## Number of parameters estimated: 10 
## AIC: 2.4169, (AIC for lm: 5.8886)
```

### SAR assumptions

As we corrected for spatial autocorrelation, we must ensure our models do not violate the autocorrelation assumption. We define a function to plot the residuals of our models and corroborate randomness.


```r
# Residual plot function for SarLm objects
residual.plot <- function(model, year) {
  plot(residuals(model) ~ model$fitted.values, xlab = "fitted values", ylab = "residuals", main = year)
  abline(h=0, lty="dotted")
  lines(lowess(model$fitted.values, residuals(model)), col="red")
}
```

We also check other assumptions as in our linear models.


```r
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
  if (class(saryear)[1] == 'SlX') {homos <- bptest(saryear)} else { homos <- bptest.Sarlm(saryear)}
  if (homos$p.value >= 0.05) {cat('Homocedastic\n')} else {cat('Heteroscedastic\n')}
  # autocorrelation
  moran <- moran.test(residuals(saryear), nb2listw(get(spatialW[[x]])[[x]]))
  if (moran$p.value >= 0.05) {cat('Non-spatially autocorrelated\n')} else {cat('Spatially autocorrelated\n')}
})) 
```

```
## 
## ***** 2004 *****
```

```
## Normal
## Homocedastic
## Non-spatially autocorrelated
## 
## ***** 2007 *****
```

```
## Normal
## Homocedastic
## Non-spatially autocorrelated
## 
## ***** 2011 *****
```

```
## Normal
## Homocedastic
## Non-spatially autocorrelated
## 
## ***** 2014 *****
```

![](SAR-Pneumonia_files/figure-html/unnamed-chunk-38-1.png)<!-- -->

```
## Normal
## Homocedastic
## Non-spatially autocorrelated
```

# Inference

## Impacts

If you are familiar with SAR models, you know they are difficult to interpret because the coefficients are a combination of direct and indirect effects. Direct effects are the effects of the spatial unit on itself. Indirect effects are the effects spatial units have on other spatial units, also known as spillover effects.

We assess the contribution of our variables trought impacts. We compute the impacts of the models' spatial lagged covariates using a simulation scheme in which the distributions for the impact measures are calculated. We also compute z-values and $p-values$ for the impacts based on the simulations.


```r
# Impacts
impactModels <- lapply(1:3, function (x) {summary(impacts(sarReg[[x]], listw = nb2listw(get(spatialW[[x]])[[x]]), R = 1000), zstats = T, short=TRUE)}) # 2014: error models doesn't have impacts
names(impactModels) <- pneuNames[1:3]
impactModels
```

```
## $pneu04
## Impact measures (mixed, exact):
##               Direct   Indirect       Total
## TEM      0.033057637 0.17619252  0.20925016
## NUT      0.033995320 0.04955298  0.08354830
## CVD      0.041353052 0.05873489  0.10008794
## ESC     -0.061147358 0.04054303 -0.02060433
## VAC      0.004716873 0.01414365  0.01886052
## lag3SMR -1.129283310 0.54226790 -0.58701541
## ========================================================
## Simulation results ( variance matrix):
## ========================================================
## Simulated standard errors
##              Direct    Indirect       Total
## TEM     0.022248261 0.027534259 0.029606795
## NUT     0.016067636 0.023060364 0.014468576
## CVD     0.020195578 0.030602317 0.015785720
## ESC     0.009979834 0.022505547 0.021916316
## VAC     0.001603976 0.002691512 0.002824084
## lag3SMR 0.223932181 0.200449213 0.128285762
## 
## Simulated z-values:
##            Direct Indirect      Total
## TEM      1.358113 6.483604  7.0503036
## NUT      1.986803 2.235221  5.7689328
## CVD      1.929253 1.995514  6.3367241
## ESC     -6.104175 1.787794 -0.9437433
## VAC      2.847645 5.396623  6.7606446
## lag3SMR -5.198515 2.840103 -4.6366684
## 
## Simulated p-values:
##         Direct     Indirect   Total     
## TEM     0.1744279  8.9557e-11 1.7852e-12
## NUT     0.0469443  0.0254028  7.9775e-09
## CVD     0.0536994  0.0459868  2.3470e-10
## ESC     1.0333e-09 0.0738093  0.3453    
## VAC     0.0044044  6.7907e-08 1.3738e-11
## lag3SMR 2.0089e-07 0.0045099  3.5407e-06
## 
## $pneu07
## Impact measures (lag, exact):
##                Direct     Indirect       Total
## log(DEP)   0.42138834 -0.025242917  0.39614542
## log(ESC)  -2.99465844  0.179392521 -2.81526592
## log(IPSE) -0.02832316  0.001696675 -0.02662648
## ========================================================
## Simulation results ( variance matrix):
## ========================================================
## Simulated standard errors
##               Direct  Indirect      Total
## log(DEP)  0.21917004 0.1213802 0.24709263
## log(ESC)  0.75793881 0.7906949 0.96796486
## log(IPSE) 0.03528195 0.0124761 0.03906829
## 
## Simulated z-values:
##               Direct   Indirect      Total
## log(DEP)   1.9017945 -0.0944859  1.6404684
## log(ESC)  -4.0034551  0.1524941 -3.0102309
## log(IPSE) -0.7945346 -0.2055323 -0.7831664
## 
## Simulated p-values:
##           Direct     Indirect Total    
## log(DEP)  0.057198   0.92472  0.1009078
## log(ESC)  6.2424e-05 0.87880  0.0026105
## log(IPSE) 0.426884   0.83716  0.4335294
## 
## $pneu11
## Impact measures (mixed, exact):
##                Direct    Indirect       Total
## CPM     -0.0277853070 -0.12407381 -0.15185911
## CVV     -1.0059640917 -4.69271943 -5.69868352
## IPSE    -0.0003690531 -0.02612434 -0.02649339
## VAC      0.0691156352  0.29292868  0.36204431
## ACU      0.6142287905  1.17897811  1.79320690
## lag3SMR  1.3562236712  2.60319614  3.95941981
## ========================================================
## Simulation results ( variance matrix):
## ========================================================
## Simulated standard errors
##             Direct   Indirect      Total
## CPM     0.05679750  1.0119438  1.0685379
## CVV     2.33278067 41.6178129 43.9451767
## IPSE    0.01113114  0.1981132  0.2092065
## VAC     0.12004322  2.1059007  2.2250464
## ACU     0.73059172 12.6560927 13.3719096
## lag3SMR 1.75485684 30.3561065 32.0796863
## 
## Simulated z-values:
##             Direct   Indirect      Total
## CPM     -0.6085231 -0.2419341 -0.2614660
## CVV     -0.5483569 -0.2306920 -0.2475834
## IPSE    -0.1542665 -0.2533476 -0.2481218
## VAC      0.7054217  0.2689983  0.2926522
## ACU      0.9465266  0.2080068  0.2485866
## lag3SMR  0.8952363  0.2049149  0.2428774
## 
## Simulated p-values:
##         Direct  Indirect Total  
## CPM     0.54284 0.80883  0.79373
## CVV     0.58345 0.81755  0.80446
## IPSE    0.87740 0.80000  0.80404
## VAC     0.48055 0.78793  0.76979
## ACU     0.34388 0.83522  0.80368
## lag3SMR 0.37066 0.83764  0.80810
```


