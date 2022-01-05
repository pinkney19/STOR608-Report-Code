#INITIAL
data <- read.csv("Group1.csv")
x <- data$x

period = 1440 
dts = t(matrix(x,nrow=period))

plot_basic <- function(dt, median_index = 12,
                       outlier_indices = c() ,
                       other_indices = c(),
                       title="",ylabel="y"){
  plot(dt[median_index,],ylim=c( min(dt) , max(dt) ),
       main=title,ylab=ylabel,type="l", lwd=4)
  for (ind in outlier_indices){
    lines(dt[ind,],lty="dashed",col="red",lwd=3)
  } 
  for (ind in other_indices){
    lines(dt[ind,],lty="dotted")
  } 
}
par(mfrow=c(2,2))
plot_basic(dts,other_indices=c(1:50), title='Original Data', ylabel = 'Y') #median curve in bold


library(fdaoutlier)
source("my_functional_boxplot.R") #FROM GEORGE


#SEQ TRANSFORM
####Sanity Check for seq_transform!####
# Transformations 
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

plot_complete <- function(dts,out,title="",ylabel=""){
  plot_basic(dts, median_index = order(out$depth_values)[50],
             outlier_indices = out$outliers,
             other_indices = setdiff(c(1:50), out$outliers),
             title=title, ylabel = ylabel  
  )
  lines( out$lower, col="blue",lwd=2 )
  lines(out$upper, col="blue",lwd=2 )
  lines( out$inf, col="magenta",lwd=5)
  lines(out$sup, col="magenta",lwd=5 )
  time= c(1:period)
  polygon( c(time,rev(time)), c(out$sup, rev(out$inf)), col=rgb(1,0,1,0.5) )
}

plot_complete(dts_T0,out_T0,title="T0", ylabel='Y')
plot_complete(dts_T1,out_T1,title="T1", ylabel='Y')
plot_complete(dts_T2,out_T2,title="T2", ylabel='Y')



#FAST MUOD
fast_muod3 <- function(){
  N=50
  a <- rep(0,N)
  b <- rep(0,N)
  r <- rep(0,N)
  Y_tilda = colMeans(dts) #point-wise mean 
  Y_tilda_bar = mean(Y_tilda)
  sy_tilda = sd(Y_tilda)
  for(i in 1:N){
    s_yi = sd(dts[i,])
    sy_tilda_yi = cov(Y_tilda,dts[i,])
    r[i] =  1 -  sy_tilda_yi/(sy_tilda*s_yi)
    b[i] =  sy_tilda_yi/sy_tilda^2
    a[i] = abs ( mean(dts[i,]) - b[i] * Y_tilda_bar )
    b[i] = abs(1 - b[i])
  }
  return(list( a=a,b=b,r=r ))
}

m_fast = fast_muod3()
a = m_fast$a
b = m_fast$b
r = m_fast$r

#cut-off values
#shape
Q3_s = quantile(r, 0.75) #Q3
IQR_s = quantile(r, 0.75)-quantile(r, 0.25)#IQR
cut_off_shape = Q3_s + 1.5*IQR_s 
r[r>cut_off_shape]

#magnitude
Q3_m = quantile(a, 0.75) #Q3
IQR_m = quantile(a, 0.75)-quantile(a, 0.25)#IQR
cut_off_mag = Q3_m + 1.5*IQR_m 

#amp
Q3_a = quantile(b, 0.75) #Q3
IQR_a = quantile(b, 0.75)-quantile(b, 0.25)#IQR
cut_off_amp = Q3_a + 1.5*IQR_a 

par(mfrow=c(2,2))
hist(a, main='', xlab='Magnitude Penalty', xlim=c(0,20),col = 'magenta')
abline(v=cut_off_mag, col='blue')
hist(b, main="", xlab='Amplitude Penalty', xlim=c(0,1), col='magenta')
abline(v=cut_off_amp, col='blue')
hist(r, breaks=10,main="", xlab='Shape Penalty', col='magenta')
abline(v=cut_off_shape,col='blue')

plot(dts[1,],type="l",ylim=c(21.5,25),ylab='Y')
for (p in c(1:50)){
  lines(dts[p,],type="l")
}
lines(dts[19,],type="l",col="magenta",lwd=2) 
lines(dts[43,],type="l",col="magenta",lwd=2)
lines(dts[37,],type="l",col="magenta",lwd=2)

#computation times
#Timing for Functional boxplot
Time1=NULL;
for(i in 1:100){
  Time1[i] = system.time(functional_boxplot(dts) )[3] #extract elapsed time
}
Av_time_fb=mean(Time1)*1000 #get time in ms -> 18.4

#Timing for seq transform
Time2=NULL;

for(i in 1:100){
  Time2[i] = system.time(seq_transform(dts) )[3] #extract elapsed time
}
Av_time_s=mean(Time2)*1000 #get time in ms -> 66.0

#Timing for fast-muod
Time3=NULL;

for(i in 1:100){
  Time3[i] = system.time(fast_muod3() )[3] #extract elapsed time
}
Av_time_m=mean(Time3)*1000 #get time in ms -> 66.0

print(Av_time_fb)
print(Av_time_s)
print(Av_time_m)
