insects = c(0.03, 0.15, 0.18, 0.44, 0.58, 0.63, 1.37, 1.91, 2.60, 3.58, 4.10, 7.46)
n = length(insects)

A = nlm(f,c(4,4))$estimate[1]
B = nlm(f,c(4,4))$estimate[2]

MLest = sum(sqrt(insects))^-1*length(insects)
postMean = (A+n)/(B+sum(sqrt(insects)))
postMode = (A+n-1)/(B+sum(sqrt(insects)))

plot(xs,dgamma(xs,a,b),type="l",ylim=c(0,1.6),bty="l",ylab="Tetthet",xlab=expression(theta))
lines(xs,dgamma(xs,a+12,(b+14.01658)))
abline(v=MLest)
abline(v=postMode)

ys = 1:1000000
sum(1/((ys)*(ys+1)))

integrate(function(x) x^-1*(1-x)^(-1/2),lower=0,upper=0.1)
xs = seq(0,1,by=0.001)
plot(xs,xs^-1*(1-xs)^(-1/2),type="l",xlim=c(0.01,0.99))


y = 20
plot(xs,dbeta(xs,1,y+0.5),type="l")
abline(v=1/y,col="blue")
abline(v=1/(1+y+0.5),col="red")

ys = seq(1:10)
plot(ys,1/ys)
points(ys,2/(2+ys),col="blue")
points(ys,1/(1+ys+0.5),col="red")


# Late Bayesianere ------------------------------------------------------------

ys = insects
tt = seq(-0.1,2.2,by=0.001)
#ys = (rexp(1000,3))^2
ys = insects
n = length(ys)
thetaML = n/sum(sqrt(ys))
top = (A+n-1)/(B+sum(sqrt(ys)))
plot(tt,dgamma(tt,A+n,B+n/thetaML),type="l",ylim=c(0,1.8),bty="l",
     xlab=expression(theta),ylab="Tetthet")
lines(tt,dnorm(tt,top,sqrt(top^2/(A+n-2))),type="l",lty=2,col="red" )
lines(tt,dnorm(tt,top,top/(sqrt(n)))  ,type="l",lty=3,col="blue")
legend("topright",c("Posterior","Lat","Ekstra lat"),
       col=c("black","red","blue"),lty=c(1,2,3),bty="n")


# Risikofunksjoner ------------------------------------------------------------

thetas = seq(0.001,0.999,by=0.01)

risky = function(theta,f) {
  ys = 1:10000
  sum((theta-f(ys))^2*(1-theta)^(ys-1)*theta)
}

A = 10
B = 10

f_ML = function(x) 1/x
f_UP = function(x) 1/(1+x/2)
f_JP = function(x) 1/(1/2+x)

risk_ML = sapply(thetas,function(theta) risky(theta,f_ML))
risk_UP = sapply(thetas,function(theta) risky(theta,f_UP))
risk_JP = sapply(thetas,function(theta) risky(theta,f_JP))
plot(thetas,risk_ML,type="l",ylim=c(0,0.17),bty="l",ylab="Risiko",xlab=expression(theta))
lines(thetas,risk_UP,col="blue")
lines(thetas,risk_JP,col="red")
legend("topright",c("ML","Uniform","Jeffreys"),bty="n",
       col=c("black","blue","red"),lty=c(1,1,1))

# Risikofunksjoner, større n --------------------------------------------------

thetas = seq(0.001,0.999,by=0.01)

N = 10000
n = 50

riski = function(theta,f,ys) {
 (theta-f(ys))^2
}

risk = function(theta,f,bb) {
  mean(apply(bb,1,function(ys) riski(theta,f,ys)))
}

A = 10
B = 10

f_ML = function(ys) n/sum(ys)
f_UP = function(ys) (n+1)/(2+sum(ys))

risk_ML = sapply(thetas,function(theta) risk(theta,f_ML,matrix(rnbinom(n*N, 1, theta) + 1,ncol=n,nrow=N)))
risk_UP = sapply(thetas,function(theta) risk(theta,f_UP,matrix(rnbinom(n*N, 1, theta) + 1,ncol=n,nrow=N)))
plot(0,type="n",bty="l",ylab="Risiko",xlab=expression(theta),xlim=c(0,1),ylim=c(0,max(c(risk_ML,risk_UP))))
lines(locfit(risk_ML~thetas),col="black")
lines(locfit(risk_UP~thetas),col="blue")
legend("topleft",c("ML","Uniform"),bty="n",
       col=c("black","blue","red"),lty=c(1,1,1))



# Gamma-parametere ------------------------------------------------------------
ys = c(1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 5, 6, 6, 14, 15, 15, 20, 75)

llik = function(p,ys){
  a = p[1]
  b = p[2]
  sum(lbeta(a+1,b+ys-1) - lbeta(a,b))
}

llik2 = function(p,ys){
  a = p[1]
  b = p[2]
  lbeta(a+n,b+sum(ys)-n) - lbeta(a,b)
}

ab = nlm(function(p) -llik(p,ys),p=c(1,1))$estimate
ab = nlm(function(p) -llik2(p,ys),p=c(1,11000))$estimate

plot(as,sapply(as,function(a) -llik2(c(a,11),ys)),type="l",col="red")
plot(as,sapply(as,function(a) -llik(c(a,11000),ys)))

a = ab[1]
b = ab[2]


ys = c(1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 5, 6, 6, 14, 15, 15, 20, 75)

llik2 = function(p,ys){
  a = p[1]
  b = p[2]
  lbeta(a+n,b+sum(ys)-n) - lbeta(a,b)
}

ab = nlm(function(p) -llik(p,ys),p=c(1310,11236))$estimate
ab = nlm(function(p) -llik2(p,ys),p=c(1310,11236))$estimate

a = ab[1]
b = ab[2]

plot(xs,dbeta(xs,a,b),type="l")


# Sample from posterior -------------------------------------------------------

ys = c(1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 5, 6, 6, 14, 15, 15, 20, 75)

a1 = 78
b1 = 4
thetas = rbeta(100,a1,b1)
#ys = rnbinom(100, 1, thetas) + 1

n = length(ys)

sampleit = function() {
  a = runif(1,0,100)
  b = runif(1,0,100)
  c(a,b,prod(beta(a+1,b+ys-1)/beta(a,b)))
}

samps = replicate(N,sampleit())
const = mean(samps[3,])

(expected_a = mean(samps[1,]*samps[3,]/const))
(expected_b = mean(samps[2,]*samps[3,]/const))

sd_a =  sqrt(mean(samps[1,]^2*samps[3,]/const) - expected_a^2)
sd_b =  sqrt(mean(samps[2,]^2*samps[3,]/const) - expected_b^2)

corr = (mean(samps[1,]*samps[2,]*samps[3,]/const) - expected_a*expected_b)/(sd_a*sd_b)

# Table of vals ---------------------------------------------------------------

parti = function(i,fun) {
  fun(sapply(1:N,function(j) dbeta()*samps[3,j])
}

a = ab[1]
b = ab[2]

col1 = ys
col2 = 1/ys
col3 = (a+1)/(a+b+ys)



# Posterior predictive --------------------------------------------------------

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

plot(yy,dppd(yy,alpha,beta),type="l",col="blue",bty="l",xlab=expression(y[13]),
     ylab="Tetthet",ylim=c(0,2),xlim=c(0,5))
lines(yy,dexp(yy,1/mean(sqrt(insects))),lty=2)
lines(yy,dppd(yy,4+n,4+sum(sqrt(ys))),col="red",lty=3)
legend("topright",c("(5.31,2.66)","(4,4)","Exp"),col=c("blue","red","black"),lty=c(1,3,2),bty="n")


yy = seq(0,20,by=0.01)
q01 = qppd(0.1,alpha,beta)
q05 = qppd(0.5,alpha,beta)
q09 = qppd(0.9,alpha,beta)
plot(yy,pppd(yy,alpha,beta),type="l",col="blue",bty="l",xlab=expression(y[13]),
     ylab="Sannsynlighet",ylim=c(0,1))
lines(yy,pppd(yy,4+n,4+sum(sqrt(ys))),col="red",lty=3)
lines(yy,pexp(yy,1/mean(sqrt(ys))),lty=2)
points(q01,pppd(q01,alpha,beta),pch=4)
points(q05,pppd(q05,alpha,beta),pch=4)
points(q09,pppd(q09,alpha,beta),pch=4)

plot(xs,qppd(xs,alpha,beta),col="blue",bty="l",xlab="Sannsynlighet",ylab="Kvantil",type="l",
     ylim=c(0,20))
lines(xs,qppd(xs,4+n,4+sum(sqrt(ys))),col="red",lty=3)
lines(xs,qexp(xs,1/mean(sqrt(ys))),lty=2)

round(cbind(qppd(c(0.1,0.5,0.9),alpha,beta),
qppd(c(0.1,0.5,0.9),n+4,4+sum(sqrt(ys))),
qexp(c(0.1,0.5,0.9),1/mean(sqrt(ys)))),4)
