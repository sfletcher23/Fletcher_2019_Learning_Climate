library(base)
setwd("C:/Works/Vietnam/Model4/OutputData/data/sim")
Allf <- matrix(scan("Allout.txt"),1200,1254)



dd= (1:1200)*0
for( i in 0:56) 
{
dd=cbind(dd,Allf[,(i*22)+7]) # [,(i*22)+DD]) DD is the subbasin number ranging from 1:22 
}
dd=dd[,-1]

summ=matrix(0,100,57)
for(i in 1:100)
summ[i,]=apply(dd[((i-1)*12+1):(i*12),],2,sum)

plot(summ[,1],type="l",col="red",xlim=c(0,100),ylim=c(min(summ),max(summ)))
for(i in 1:56)
lines(summ[,i],type="l")
lines(summ[,1],type="l",col="RED")
win.graph()
boxplot(t(summ))