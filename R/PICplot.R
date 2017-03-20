PICplot = function(PICs,index)
{
  PIC <- PICs$PICs[[index]]
  rt <- PIC[,1]
  inte <- PIC[,2]
  plot(rt,inte,wd=2,pch=15,col='blue',cex.lab=1.5,cex.axis=1.3,font=2,
       type="l",xlab="Retention time",ylab="Intensity",
       main="Pure Ion Chromatogram",cex.main=1.5)
}

stem <- function(x,y,pch=5,linecol=1,col='blue',cex.lab=1.2,cex.axis=1.3,font=2,...){
  if (missing(y)){
    y = x
    x = 1:length(x) }
  plot(x,y,pch=pch,type='n',col=col,cex.lab=cex.lab,cex.axis=cex.axis,font=font,,xlab="M/Z",ylab="Intensity")
  for (i in 1:length(x)){
    lines(c(x[i],x[i]), c(0,y[i]),pch=pch,col=col,cex.lab=cex.lab,cex.axis=cex.axis,font=font)
  }
  lines(c(x[1]-2,x[length(x)]+2), c(0,0),pch=pch,col=col,cex.lab=cex.lab,cex.axis=cex.axis,font=font)
}
