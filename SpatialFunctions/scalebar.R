assign("scalebar",
function(loc,length,unit="km",division.cex=.8,...) {
if(missing(loc)) stop("loc is missing")
if(missing(length)) stop("length is missing")
x <- c(0,length/c(4,2,4/3,1),length*1.1)+loc[1]
y <- c(0,length/(10*3:1))+loc[2]
cols <- rep(c("black","white"),2)
for (i in 1:4) rect(x[i],y[1],x[i+1],y[2],col=cols[i])
for (i in 1:5) segments(x[i],y[2],x[i],y[3])
labels <- x[c(1,3)]-loc[1]
labels <- append(labels,paste(x[5]-loc[1],unit))
text(x[c(1,3,5)],y[4],labels=labels,adj=.5,cex=division.cex)
}
)