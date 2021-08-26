assign("moran.cluster",
   function(x, listw, zero.policy = FALSE, shp, significant=T, ...){
  require(RColorBrewer)
  require(spdep)
mor.dat <- localmoran(x, listw)
wx<-lag.listw(listw, x)
lag.z<-scale(wx, center=T, scale=T)
dat.z<-scale(x, center=T, scale=T)

mor.dat1<-data.frame(mor.dat,lag.z,dat.z)
names(mor.dat1)<-c("Ii","E.Ii","Var.Ii","Z.Ii","Pr.Zi","WY","Y")

mor.dat1$cluster<-"UN"
mor.dat1$cluster[mor.dat1[,"Z.Ii"]>=1.644854&mor.dat1[,"Y"]>0&mor.dat1[,"WY"]>0]<-"HH"     # both z scores are "high"
mor.dat1$cluster[mor.dat1[,"Z.Ii"]>=1.644854&mor.dat1[,"Y"]<0&mor.dat1[,"WY"]<0]<-"LL"     # both z scores are "low"
mor.dat1$cluster[mor.dat1[,"Z.Ii"]<=-1.644854&mor.dat1[,"Y"]>0&mor.dat1[,"WY"]<0]<-"HL"    # one z score "high", the other "low"
mor.dat1$cluster[mor.dat1[,"Z.Ii"]<=-1.644854&mor.dat1[,"Y"]<0&mor.dat1[,"WY"]>0]<-"LH"    # one z score "low", the other "high"

cols<-brewer.pal(5, "RdBu")
mor.dat1$col[mor.dat1$cluster=="UN"]<-cols[3]
mor.dat1$col[mor.dat1$cluster=="HH"]<-cols[1]
mor.dat1$col[mor.dat1$cluster=="LL"]<-cols[5]
mor.dat1$col[mor.dat1$cluster=="HL"]<-cols[2]
mor.dat1$col[mor.dat1$cluster=="LH"]<-cols[4]
mor.dat1

par(pty="s", mfrow=c(1,2),mar=c(0,0,0,0))
P1 <- plot(shp, col=mor.dat1$col, ...)
legend("bottomright", legend=c(paste("Not Significant  ","(",length(mor.dat1$cluster[mor.dat1$cluster=="UN"]),")",sep="",collapse=""), 
                            paste("High-High  ","(",length(mor.dat1$cluster[mor.dat1$cluster=="HH"]),")",sep="",collapse=""), 
                            paste("Low-Low  ","(",length(mor.dat1$cluster[mor.dat1$cluster=="LL"]),")",sep="",collapse=""), 
                            paste("Low-High  ","(",length(mor.dat1$cluster[mor.dat1$cluster=="LH"]),")",sep="",collapse=""),
                            paste("High-Low  ","(",length(mor.dat1$cluster[mor.dat1$cluster=="HL"]),")",sep="",collapse="")),             
fill=c(cols[3],cols[1],cols[5],cols[4],cols[2]), title = "LISA Cluster Map", bty="n", cex=0.6, y.intersp=0.8)

if (significant) {
mor.dat1$prob<-"UN"
mor.dat1$prob[mor.dat1[,"Pr.Zi"]<=0.05&mor.dat1[,"Pr.Zi"]>0.01|mor.dat1[,"Pr.Zi"]>0.95&mor.dat1[,"Pr.Zi"]<=0.99]<-"5%"
mor.dat1$prob[mor.dat1[,"Pr.Zi"]<0.01&mor.dat1[,"Pr.Zi"]>0.001|mor.dat1[,"Pr.Zi"]>0.99&mor.dat1[,"Pr.Zi"]<=0.999]<-"1%"
mor.dat1$prob[mor.dat1[,"Pr.Zi"]<0.001&mor.dat1[,"Pr.Zi"]>0.0001|mor.dat1[,"Pr.Zi"]>0.999&mor.dat1[,"Pr.Zi"]<=0.9999]<-"0.1%"
mor.dat1$prob[mor.dat1[,"Pr.Zi"]<0.0001&mor.dat1[,"Pr.Zi"]>0|mor.dat1[,"Pr.Zi"]>0.9999&mor.dat1[,"Pr.Zi"]<=1]<-"0.01%"


colsp<-brewer.pal(5, "Greens")
mor.dat1$col1[mor.dat1$prob=="UN"]<-colsp[1]
mor.dat1$col1[mor.dat1$prob=="5%"]<-colsp[2]
mor.dat1$col1[mor.dat1$prob=="1%"]<-colsp[3]
mor.dat1$col1[mor.dat1$prob=="0.1%"]<-colsp[4]
mor.dat1$col1[mor.dat1$prob=="0.01%"]<-colsp[5]
mor.dat1

P2 <- plot(shp, col=mor.dat1$col1, ...)
legend("bottomright", legend=c(paste("Not Significant  ","(",length(mor.dat1$prob[mor.dat1$prob=="UN"]),")",sep="",collapse=""),
                            paste("p=0.05  ","(",length(mor.dat1$prob[mor.dat1$prob=="5%"]),")",sep="",collapse=""), 
                            paste("p=0.01 ","(",length(mor.dat1$prob[mor.dat1$prob=="1%"]),")",sep="",collapse=""), 
                            paste("p=0.001  ","(",length(mor.dat1$prob[mor.dat1$prob=="0.1%"]),")",sep="",collapse=""),
                            paste("p=0.0001  ","(",length(mor.dat1$prob[mor.dat1$prob=="0.01%"]),")",sep="",collapse="")), 
       bty="n", fill=colsp, title = "LISA Significance Map", cex=0.6, y.intersp=0.8)
}
})