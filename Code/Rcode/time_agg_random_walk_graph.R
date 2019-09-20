#########################################################################################
# Graphs to show how time aggregation can result in positive autocorrelation
# Figure 1 (last graph produced) is the one that appears in the paper
#########################################################################################

library(ggplot2)
library(gridExtra)

figures_dir = "../Figures/"


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

#dev.new()
pdf(paste(figures_dir, "TimeAggExample.pdf",sep=""))
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
#dev.copy(pdf, paste(figures_dir, "TimeAggExample.pdf",sep=""))
dev.off()
# #########################################################################################
# # #plot the above for Danish slides
# dev.new()
# par(mar=c(5,5,5,5), mfcol=c(1,1),cex.axis=1.5,cex.lab=2)
# plot(underlying_func, time, xlim= c(0,N_period),ylim=c(0,1.05),col="black",lty="solid", col.points=FALSE, verticals=TRUE,xlab="Time",main="Time Aggregation",ylab="Income",yaxt = "n")
# axis(side = 2, at = c(0.0,0.5,1.0))
# legend(0,0.8,legend=c("Income Flow",""),lty=c("solid","dashed"),bty = "n")
# dev.copy(pdf, paste(figures_dir, "TimeAggExample1.pdf",sep=""))
# dev.off()
# dev.new()
# par(mar=c(5,5,5,5), mfcol=c(1,1),cex.axis=1.5,cex.lab=2)
# plot(underlying_func, time, xlim= c(0,N_period),ylim=c(0,1.05),col="black",lty="solid", col.points=FALSE, verticals=TRUE,xlab="Time",main="Time Aggregation",ylab="Income",yaxt = "n")
# axis(side = 2, at = c(0.0,0.5,1.0))
# lines(time_agg_func, time, xlim= c(0,N_period),ylim=c(0,1),col="black",lty="dashed", col.points=FALSE, verticals=TRUE)
# legend(0,0.8,legend=c("Income Flow","Observed Income"),lty=c("solid","dashed"),bty = "n")
# dev.copy(pdf, paste(figures_dir, "TimeAggExample2.pdf",sep=""))
# dev.off()

#########################################################################################
# Now plot the autocorrelation of N-subperiod example and show how quickly it converges to 1/4
max_subperiods = 12
Num_subperiods = 1:max_subperiods
autocorr = (Num_subperiods**2-1)/(2*(2*Num_subperiods**2 +1 ))

#dev.new()
pdf(paste(figures_dir, "InducedAutocorrelation.pdf",sep=""))
par(mfcol=c(1,1))
plot(Num_subperiods,autocorr,ylim = c(0,0.3),xlab="Number of sub-periods",ylab="Autocorrelation",main="Induced Autocorrelation",yaxt = "n",xaxt = "n")
axis(side = 2, at = c(0.0,0.05,0.1,0.15,0.2,0.25,0.3))
axis(side = 1, at = 1:max_subperiods)
constant = Num_subperiods*0+0.25
lines(Num_subperiods,constant,lty="dotted")
#dev.copy(pdf, paste(figures_dir, "InducedAutocorrelation.pdf",sep=""))
dev.off()

#########################################################################################
# Figure 1 in the paper
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


pdf(paste(figures_dir, "TimeAgg_simple.pdf",sep=""))
par(mar=c(5,7,2,2), mfcol=c(1,1),cex.axis=1.5,cex.lab=2)
plot(underlying_func, time, xlim= c(0,N_period),ylim=c(0,110000),col="black",lty="solid",lwd=4, col.points=FALSE, verticals=TRUE,xlab="Time",main="",ylab="",yaxt = "n",las=1)
axis(side = 2, at = c(0,25000,50000,75000,100000), las=1,label=c("$0","$25,000","$50,000","$75,000","$100,000"))
lines(time_agg_func, time, xlim= c(0,N_period),ylim=c(0,110000),col="blue",lty="dashed",lwd=4, verticals=TRUE ,col.points=FALSE)
legend(0.3,40000,legend=c("Permanent Income Flow","Observed Annual Income"),lty=c("solid","dashed"),lwd=c(4,4),col=c("black","blue"),bty = "n", cex=1.7)
dev.off()
