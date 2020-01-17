#########################################################################################
# Draw a time aggregated random walk (4 subperiods)
library(ggplot2)
library(gridExtra)

figures_dir = "./Figures/"
figures_dir = "C:/Users/edmun/OneDrive/Documents/Research/BPP_PSID_TimeAgg/Code/Figures/"


N_sub =2
N_period =3
N_subperiods =N_sub*N_period
time = (1:(N_subperiods))/N_sub
shocks = c(0,0,0,1,0,0)
underlying = 50000+50000*cumsum(shocks)
time_agg = underlying*0.0
for (T in 1:N_period){
  t = (T-1)*N_sub+1
  time_agg[t:(t+N_sub-1)]=mean(underlying[t:(t+N_sub-1)])
}

time_agg_func=stepfun(time[1:(N_subperiods-1)], time_agg)
underlying_func=stepfun(time[1:(N_subperiods-1)], underlying)


pdf(paste(figures_dir, "TimeAggAER_simple.pdf",sep=""))
par(mar=c(5,7,2,2), mfcol=c(1,1),cex.axis=1.5,cex.lab=2)
plot(underlying_func, time, xlim= c(0,N_period),ylim=c(0,110000),col="black",lty="solid",lwd=4, col.points=FALSE, verticals=TRUE,xlab="Time",main="",ylab="",yaxt = "n",las=1)
axis(side = 2, at = c(0,25000,50000,75000,100000), las=1,label=c("$0","$25,000","$50,000","$75,000","$100,000"))
lines(time_agg_func, time, xlim= c(0,N_period),ylim=c(0,110000),col="blue",lty="dashed",lwd=4, verticals=TRUE ,col.points=FALSE)
legend(0.3,40000,legend=c("Permanent Income Flow","Observed Annual Income"),lty=c("solid","dashed"),lwd=c(4,4),col=c("black","blue"),bty = "n", cex=1.7)
dev.off()
