# Oppgave c) (Posterior) ------------------------------------------------------

mu1 = 1.1
mu2 = 2.2
tau0 = 66
n = 100
yy = c(rnorm(tau0,mu1,1),rnorm(n-tau0,mu2,1))

dens = function(yy,tau) {
  prod(dnorm(yy[1:tau],mu1,1))*prod(dnorm(yy[(tau+1):n],mu2,1))
}

C = sum(sapply(1:(n-1),function(tau) dens(yy,tau)))

plot(45:90,sapply(45:90,function(tau) dens(yy,tau)/C),bty="l",
     xlab=expression(tau),ylab="Posterior",pch=20)
abline(v=66,lty=2,col="grey")

# Oppgave c) (Frekventistiske egenskaper) -------------------------------------
N = 100000

maxes = c()
means = c()

for (i in 1:N) {
  yy    = c(rnorm(tau0,mu1,1),rnorm(n-tau0,mu2,1))
  vals  = sapply(1:(n-1),function(tau) dens(yy,tau))
  C     = sum(vals)
  means = c(means,sum((1:(n-1))*vals/C))
  maxes = c(maxes,which.max(vals))
}

is = min(maxes[maxes>0]):max(maxes[maxes>0])
ts = seq(1,100,by=0.01)

plot(is,sapply(is,function(i) mean(maxes==i)),type="h",bty="l",
     col="blue",xlab=expression(widehat(tau)[mode]),ylab="Posterior",
     xlim=c(50,82),ylim=c(0,mean(maxes==66)))
lines(density(means,adjust=2),type="l",bty="l",col="red",main=NA,
      xlab=expression(widehat(tau)),ylab="Tetthet / Sannsynlighet",
      xlim=c(50,82),ylim=c(0,mean(maxes==66)))

b = mean(abs(means-66))
lines(ts,1/(2*b)*exp(-abs(ts-66)/b),lty=2,col="black")
legend("topright",c(expression(widehat(tau)[mean]),
                    expression(widehat(tau)[mode]),"Laplace"),bty="n",
       col=c("red","blue","black"),lty=c(1,1,2))

# Oppgave d) (Maximum likelihood) ---------------------------------------------

data = c(4, 5, 3, 1, 4, 4, 1, 5, 5, 3, 4, 2,
         5, 2, 2, 3, 4, 2, 1, 3, 2, 2, 1, 1,
         1, 1, 3, 0, 0, 1, 0, 1, 1, 0, 0, 3,
         1, 0, 3, 2, 2, 0, 1, 1, 1, 0, 1, 0,
         1, 0, 0)

n = length(data)

taus = 1:n

ml = function(tau) {
  thetaL = mean(data[1:tau])
  thetaR = mean(data[(tau+1):n])  
  loglik = sum(dpois(data[1:tau],thetaL)) + 
           sum(dpois(data[(tau+1):(n+1)],thetaR),na.rm=TRUE)
  loglik
}

thetas = function(tau) {
  thetaL = mean(data[1:tau])
  thetaR = mean(data[(tau+1):(n+1)],na.rm=TRUE)  
  c(thetaL,thetaR)
}

mls = sapply(taus,ml)
thetaLs = sapply(taus,thetas)[1,]
thetaRs = sapply(taus,thetas)[2,]

plot(taus,mls,bty="l",ylab="Log likelihood",xlab=expression(tau),pch=20)
abline(v=22,lty=2,col="grey")
plot(taus,thetaLs,col="blue",ylim=c(0,4.5),bty="l",
     ylab=expression(theta),xlab=expression(tau),pch=20)
points(taus,thetaRs,col="red",pch=20)
abline(h=min(thetaLs),lty=2,col="blue")
abline(h=max(thetaRs,na.rm=TRUE),lty=2,col="red")
legend("topright",c(expression(theta[L]),expression(theta[R])),bty="n",
       col=c("blue","red"),pch=c(20,20))

# Opppgave e ------------------------------------------------------------------

a = 1
b = 1

posterior = function(tau) {
  sum1 = sum(data[1:tau])
  sum2 = sum(data[(tau+1):n+1],na.rm=TRUE)
  b^(2*a)/(n*prod(factorial(data))*gamma(a)^2)*
    gamma(sum1+a)/(tau+b)^(sum1+a)*
    gamma(sum2+a)/(n-tau+b)^(sum2+a)  
}


prop = sum(sapply(taus,posterior))

posts = sapply(taus,posterior)/prop

plot(taus,sapply(taus,posterior)/prop,xlab=expression(tau),pch=20,
     ylab="Posterior",bty="l")

# -----------------------------------------------------------------------------

library("locfit")
N = 10000000
a = 1
b = 1

sampleGammas = function(tau,count) {
  sum1 = sum(data[1:tau])
  sum2 = sum(data[(tau+1):n+1],na.rm=TRUE)
  thetaR = rgamma(count,sum2+a,n-tau+b) 
  thetaL = rgamma(count,sum1+a,tau+b)
  thetaR/thetaL
}

counts = floor(posts*N)
res = unlist(sapply(2:50,function(i)  sampleGammas(i,counts[i])))
quants = round(quantile(res,c(0.05,0.95)),3)

plot(1,type="n",bty="l",xlab=expression(theta[R]/theta[L]),ylab="Posterior",
     xlim=c(min(res),max(res)),ylim=c(0,6),
     sub=paste0("CI: (",quants[1],", ",quants[2],")"))
lines(locfit(~res),main=NA,sub=NA,bty="l")
abline(v=median(res),lty=2,col="grey")
abline(v=quantile(res,0.05),col="purple",lty=3)
abline(v=quantile(res,0.95),col="purple",lty=3)