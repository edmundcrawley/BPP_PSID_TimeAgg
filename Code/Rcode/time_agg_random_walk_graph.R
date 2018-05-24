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

par(mfcol=c(2,2))
plot(underlying_func, time, xlim= c(0,N_period),ylim=c(0,1),col="black",lty="solid", col.points=FALSE, verticals=TRUE,xlab="Time",ylab="Income",main="Underlying with shock at time 1",yaxt = "n")
axis(side = 2, at = c(0.0,0.5,1.0))
plot(time_agg_func, time, xlim= c(0,N_period),ylim=c(0,1),col="black",lty="dashed", col.points=FALSE, verticals=TRUE,xlab="Time",ylab="Income",main="Time aggregated with shock at time 1",yaxt = "n")
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


plot(underlying_func, time, xlim= c(0,N_period),ylim=c(0,1.05),col="black",lty="solid", col.points=FALSE, verticals=TRUE,xlab="Time",ylab="Income",main="Underlying with shock at time 1.5",yaxt = "n")
axis(side = 2, at = c(0.0,0.5,1.0))
plot(time_agg_func, time, xlim= c(0,N_period),ylim=c(0,1),col="black",lty="dashed", col.points=FALSE, verticals=TRUE,xlab="Time",ylab="Income",main="Time aggregated with shock at time 1.5",yaxt = "n")
axis(side = 2, at = c(0.0,0.5,1.0))
dev.copy(png, paste(figures_dir, "TimeAggExample.png",sep=""))
dev.off()



