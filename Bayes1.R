# Oppgave f) ------------------------------------------------------------------
# ML, posterior forventing, posterior topppunkt
insects  = c(0.03, 0.15, 0.18, 0.44, 0.58, 0.63, 
             1.37, 1.91, 2.60, 3.58, 4.10, 7.46)
n        = length(insects)
MLest    = sum(sqrt(insects))^-1*length(insects)
postMean = (A+n)/(B+sum(sqrt(insects)))
postMode = (A+n-1)/(B+sum(sqrt(insects)))

# Oppgave g) ------------------------------------------------------------------
# Late Bayesianere

ys = insects
tt = seq(-0.1,2.2,by=0.001)
ys = insects
n  = length(ys)
thetaML = n/sum(sqrt(ys))
top = (A+n-1)/(B+sum(sqrt(ys)))
plot(tt,dgamma(tt,A+n,B+n/thetaML),type="l",ylim=c(0,1.8),bty="l",
     xlab=expression(theta),ylab="Tetthet")
lines(tt,dnorm(tt,top,sqrt(top^2/(A+n-2))),type="l",lty=2,col="red" )
lines(tt,dnorm(tt,top,top/(sqrt(n)))  ,type="l",lty=3,col="blue")
legend("topright",c("Posterior","Lat","Ekstra lat"),
       col=c("black","red","blue"),lty=c(1,2,3),bty="n")

# Oppgave h) ------------------------------------------------------------------
# Posterior predictive density 

ys = insects
a = A
b = B
n = length(ys)

alpha = a + n
beta = b + sum(sqrt(ys))

dppd = function(x,alpha,beta) {
  alpha/(2*sqrt(x))*(beta/(beta+sqrt(x)))^alpha*1/(beta+sqrt(x))
}

pppd = function(x,alpha,beta) {
  1 - (beta/(beta+sqrt(x)))^alpha
}

qppd = function(x,alpha,beta) {
  (beta*(1-x)^(-1/alpha)-beta)^2
}

plot(yy,dppd(yy,alpha,beta),type="l",col="blue",bty="l",
     xlab=expression(y[13]),ylab="Tetthet",ylim=c(0,2),xlim=c(0,5))
lines(yy,dexp(yy,1/mean(sqrt(insects))),lty=2)
lines(yy,dppd(yy,4+n,4+sum(sqrt(ys))),col="red",lty=3)
legend("topright",c("(5.31,2.66)","(4,4)","Exp"),
       col=c("blue","red","black"),lty=c(1,3,2),bty="n")

yy = seq(0,20,by=0.01)
q01 = qppd(0.1,alpha,beta)
q05 = qppd(0.5,alpha,beta)
q09 = qppd(0.9,alpha,beta)
plot(yy,pppd(yy,alpha,beta),type="l",col="blue",bty="l",
     xlab=expression(y[13]),ylab="Sannsynlighet",ylim=c(0,1))
lines(yy,pppd(yy,4+n,4+sum(sqrt(ys))),col="red",lty=3)
lines(yy,pexp(yy,1/mean(sqrt(ys))),lty=2)
points(q01,pppd(q01,alpha,beta),pch=4)
points(q05,pppd(q05,alpha,beta),pch=4)
points(q09,pppd(q09,alpha,beta),pch=4)

plot(xs,qppd(xs,alpha,beta),col="blue",bty="l",
     xlab="Sannsynlighet",ylab="Kvantil",type="l",ylim=c(0,20))
lines(xs,qppd(xs,4+n,4+sum(sqrt(ys))),col="red",lty=3)
lines(xs,qexp(xs,1/mean(sqrt(ys))),lty=2)

round(cbind(qppd(c(0.1,0.5,0.9),alpha,beta),
qppd(c(0.1,0.5,0.9),n+4,4+sum(sqrt(ys))),
qexp(c(0.1,0.5,0.9),1/mean(sqrt(ys)))),4)