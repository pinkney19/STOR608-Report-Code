library(rainbow)
X=ElNino_OISST_region_1and2
V=matrix(data=NA, nrow=37, ncol=12)
for (i in 1:nrow(V)){
  V[i,] = as.numeric(X$y[,i])
}

plot_basic <- function(dt, median_index = 9,
                       outlier_indices = c() ,
                       other_indices = c(),
                       title="",ylabel="y", xlabel="x"){
  plot(dt[median_index,],ylim=c( min(dt) , max(dt) ),
       main=title,ylab=ylabel,xlab=xlabel,type="l", lwd=4)
  for (ind in outlier_indices){
    lines(dt[ind,],lty="dashed",col="red",lwd=3)
  } 
  for (ind in other_indices){
    lines(dt[ind,],lty="dotted")
  } 
}
par(mfrow=c(2,2))

plot(X,ylab='Sea Surface Temperature', main='SST Data')


library(fdaoutlier)
source("my_functional_boxplot.R") 


#SEQ TRANSFORM ON SEA DATA
# Transformations 

dts=V
center_curves <- function(dt) dt - rowMeans(dt) 
normalize_curves <- function(dt) dt/sqrt(rowSums(dt^2))

### T0 
dts_T0 = dts 
out_T0 = my_functional_boxplot(dts_T0)

### T1 
dts_T1 = center_curves(dts_T0)
out_T1 = my_functional_boxplot(dts_T1)

### T2
dts_T2 = normalize_curves(dts_T1)
out_T2 = my_functional_boxplot(dts_T2)

plot_complete <- function(dts,out,title="",ylabel="",xlabel=""){
  plot_basic(dts, median_index = order(out$depth_values)[37],
             outlier_indices = out$outliers,
             other_indices = setdiff(c(1:37), out$outliers),
             title=title, ylabel = ylabel, xlabel = xlabel
  )
  lines( out$lower, col="blue",lwd=2 )
  lines(out$upper, col="blue",lwd=2 )
  lines( out$inf, col="magenta",lwd=5)
  lines(out$sup, col="magenta",lwd=5 )
  time= c(1:12)
  polygon( c(time,rev(time)), c(out$sup, rev(out$inf)), col=rgb(1,0,1,0.5) )
}

plot_complete(dts_T0,out_T0,title="T0", ylabel='Sea Surface Temperature', xlabel='Month')
plot_complete(dts_T1,out_T1,title="T1", ylabel='Sea Surface Temperature', xlabel='Month')
plot_complete(dts_T2,out_T2,title="T2", ylabel='Sea Surface Temperature', xlabel='Month')

#Functional Bagplot
par(mfrow=c(1,2))

fboxplot(X, plot.type = "functional", type = "bag", projmethod = "PCAproj",plotlegend = T,ylab='Sea surface temperature')

fboxplot(X, plot.type="bivariate", type="bag", projmethod = "PCAproj",
         ylim=c(-10,20), xlim=c(-10,20))



###TIMINGS###


Time1=NULL;
for(i in 1:100){
  Time1[i] = system.time(functional_boxplot(V) )[3] #extract elapsed time
}
Av_time_fb=mean(Time1)*1000 #get time in ms -> 18.4

#Timing for seq transform
Time2=NULL;

for(i in 1:100){
  Time2[i] = system.time(seq_transform(V, c('T0','T1','T2'), depth_method = 'mbd'))[3] #extract elapsed time
}
Av_time_s=mean(Time2)*1000 #get time in ms  

#Timing for bagplot
Time3=NULL;

for(i in 1:100){
  Time3[i] = system.time(fboxplot(X, plot.type = "functional", type = "bag", projmethod = "PCAproj")
  )[3] 
}
Av_time_m=mean(Time3)*1000 #get time in ms -> 66.0

print(Av_time_fb)
print(Av_time_s)
print(Av_time_m)