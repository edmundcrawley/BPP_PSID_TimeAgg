#########################################################################################
# Draw a time aggregated random walk (4 subperiods)
library(ggplot2)
library(gridExtra)

figures_dir = "C:/Users/edmun/OneDrive/Documents/Research/BPP_PSID_TimeAgg/Code/Figures/"

N_sub =2
N_period =3
N_subperiods =N_sub*N_period
time = (1:(N_subperiods))/N_sub
#shocks = rnorm(N_subperiods,0,1)
shocks = c(0,0,1,0,0,0)
underlying = cumsum(shocks)
time_agg = underlying*0.0
for (T in 1:N_period){
  t = (T-1)*N_sub+1
  time_agg[t:(t+N_sub-1)]=mean(underlying[t:(t+N_sub-1)])
}



time_agg_func=stepfun(time[1:(N_subperiods-1)], time_agg)
underlying_func=stepfun(time[1:(N_subperiods-1)], underlying)

dev.new()
par(mar=c(5,5,5,5), mfcol=c(2,2),cex.axis=1.2,cex.lab=1.5)
plot(underlying_func, time, xlim= c(0,N_period),ylim=c(0,1),col="black",lty="solid", col.points=FALSE, verticals=TRUE,xlab="Time",ylab="Income Flow",main="Underlying with shock at time 1",yaxt = "n")
axis(side = 2, at = c(0.0,0.5,1.0))
plot(time_agg_func, time, xlim= c(0,N_period),ylim=c(0,1),col="black",lty="dashed", col.points=FALSE, verticals=TRUE,xlab="Time",ylab="Observed Income",main="Time aggregated with shock at time 1",yaxt = "n")
axis(side = 2, at = c(0.0,0.5,1.0))

shocks = c(0,0,0,1,0,0)
underlying = cumsum(shocks)
time_agg = underlying*0.0
for (T in 1:N_period){
  t = (T-1)*N_sub+1
  time_agg[t:(t+N_sub-1)]=mean(underlying[t:(t+N_sub-1)])
}

time_agg_func=stepfun(time[1:(N_subperiods-1)], time_agg)
underlying_func=stepfun(time[1:(N_subperiods-1)], underlying)


plot(underlying_func, time, xlim= c(0,N_period),ylim=c(0,1.05),col="black",lty="solid", col.points=FALSE, verticals=TRUE,xlab="Time",ylab="Income Flow",main="Underlying with shock at time 1.5",yaxt = "n")
axis(side = 2, at = c(0.0,0.5,1.0))
plot(time_agg_func, time, xlim= c(0,N_period),ylim=c(0,1),col="black",lty="dashed", col.points=FALSE, verticals=TRUE,xlab="Time",ylab="Observed Income",main="Time aggregated with shock at time 1.5",yaxt = "n")
axis(side = 2, at = c(0.0,0.5,1.0))
dev.copy(pdf, paste(figures_dir, "TimeAggExample.pdf",sep=""))
dev.off()
#########################################################################################
#plot the above for Danish slides
figures_denmark_dir = "C:/Users/edmun/OneDrive/Documents/Research/Denmark/IncomeUncertaintyGit/Code/Rcode/Figures/"
dev.new()
par(mar=c(5,5,5,5), mfcol=c(1,1),cex.axis=1.5,cex.lab=2)
plot(underlying_func, time, xlim= c(0,N_period),ylim=c(0,1.05),col="black",lty="solid", col.points=FALSE, verticals=TRUE,xlab="Time",main="Time Aggregation",ylab="Income",yaxt = "n")
axis(side = 2, at = c(0.0,0.5,1.0))
legend(0,0.8,legend=c("Income Flow",""),lty=c("solid","dashed"),bty = "n")
dev.copy(pdf, paste(figures_denmark_dir, "TimeAggExample1.pdf",sep=""))
dev.off()
dev.new()
par(mar=c(5,5,5,5), mfcol=c(1,1),cex.axis=1.5,cex.lab=2)
plot(underlying_func, time, xlim= c(0,N_period),ylim=c(0,1.05),col="black",lty="solid", col.points=FALSE, verticals=TRUE,xlab="Time",main="Time Aggregation",ylab="Income",yaxt = "n")
axis(side = 2, at = c(0.0,0.5,1.0))
lines(time_agg_func, time, xlim= c(0,N_period),ylim=c(0,1),col="black",lty="dashed", col.points=FALSE, verticals=TRUE)
legend(0,0.8,legend=c("Income Flow","Observed Income"),lty=c("solid","dashed"),bty = "n")
dev.copy(pdf, paste(figures_denmark_dir, "TimeAggExample2.pdf",sep=""))
dev.off()

#########################################################################################
# Now plot the autocorrelation of N-subperiod example and show how quickly it converges to 1/4
max_subperiods = 12
Num_subperiods = 1:max_subperiods
autocorr = (Num_subperiods**2-1)/(2*(2*Num_subperiods**2 +1 ))

dev.new()
par(mfcol=c(1,1))
plot(Num_subperiods,autocorr,ylim = c(0,0.3),xlab="Number of sub-periods",ylab="Autocorrelation",main="Induced Autocorrelation",yaxt = "n",xaxt = "n")
axis(side = 2, at = c(0.0,0.05,0.1,0.15,0.2,0.25,0.3))
axis(side = 1, at = 1:max_subperiods)
constant = Num_subperiods*0+0.25
lines(Num_subperiods,constant,lty="dotted")
dev.copy(pdf, paste(figures_dir, "InducedAutocorrelation.pdf",sep=""))
dev.off()

