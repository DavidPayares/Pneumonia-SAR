assign("getis.cluster",
   function(x, listw, zero.policy = FALSE, shp, significant=T, ...){
  require(RColorBrewer)
  require(spdep)
get.dat <- localG(x, listw)

get.dat1<-data.frame(as.vector(get.dat))
names(get.dat1)<-c("Z.Gi")


get.dat1$cluster<-"UN"
get.dat1$cluster[get.dat1[,"Z.Gi"]>=1.644854]<-"HH"     # both z scores are "high"
get.dat1$cluster[get.dat1[,"Z.Gi"]<=-1.644854]<-"LL"     # both z scores are "low"


cols<-brewer.pal(5, "RdBu")
get.dat1$col[get.dat1$cluster=="UN"]<-cols[3]
get.dat1$col[get.dat1$cluster=="HH"]<-cols[1]
get.dat1$col[get.dat1$cluster=="LL"]<-cols[5]
get.dat1

par(pty="s",mar=c(0,0,0,0))
P1 <- plot(shp, col=get.dat1$col, ...)
legend(locator(1), legend=c(paste("Not Significant  ",  "(",length(get.dat1$cluster[get.dat1$cluster=="UN"]),")",sep="",collapse=""), 
                            paste("High-High  ",  "(",length(get.dat1$cluster[get.dat1$cluster=="HH"]),")",sep="",collapse=""), 
                            paste("Low-Low  ",  "(",length(get.dat1$cluster[get.dat1$cluster=="LL"]),")",sep="",collapse="")),
       fill=c(cols[3],cols[1],cols[5]), title = "Gi* Cluster Map",, bty="n", cex=1.2, y.intersp=0.8)

if (significant) {
get.dat1$prob<-"UN"
get.dat1$prob[get.dat1[,"Z.Gi"]>=1.644854|get.dat1[,"Z.Gi"]<=-1.644854]<-"5%"
get.dat1$prob[get.dat1[,"Z.Gi"]>=2.326348|get.dat1[,"Z.Gi"]<=-2.326348]<-"1%"
get.dat1$prob[get.dat1[,"Z.Gi"]>=3.090232|get.dat1[,"Z.Gi"]<=-3.090232]<-"0.1%"
get.dat1$prob[get.dat1[,"Z.Gi"]>=3.719016|get.dat1[,"Z.Gi"]<=-3.719016]<-"0.01%"


colsp<-brewer.pal(5, "Greens")
get.dat1$col1[get.dat1$prob=="UN"]<-colsp[1]
get.dat1$col1[get.dat1$prob=="5%"]<-colsp[2]
get.dat1$col1[get.dat1$prob=="1%"]<-colsp[3]
get.dat1$col1[get.dat1$prob=="0.1%"]<-colsp[4]
get.dat1$col1[get.dat1$prob=="0.01%"]<-colsp[5]
get.dat1

x11()
par(pty="s",mar=c(0,0,0,0))
P2 <- plot(shp, col=get.dat1$col1, ...)
legend(locator(1), legend=c(paste("Not Significant  ", "(",length(get.dat1$prob[get.dat1$prob=="UN"]),")",sep="",collapse=""), 
                            paste("p=0.05  ",  "(",length(get.dat1$prob[get.dat1$prob=="5%"]),")",sep="",collapse=""), 
                            paste("p=0.01  ",  "(",length(get.dat1$prob[get.dat1$prob=="1%"]),")",sep="",collapse=""), 
                            paste("p=0.001  ",  "(",length(get.dat1$prob[get.dat1$prob=="0.1%"]),")",sep="",collapse=""),
                            paste("p=0.0001  ",  "(",length(get.dat1$prob[get.dat1$prob=="0.01%"]),")",sep="",collapse="")), 
                            bty="n", fill=colsp, title = "Gi* Significance Map", cex=1.2, y.intersp=0.8)
}
})