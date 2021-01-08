# Oppgave d) ------------------------------------------------------------------
# Beregning av risikofunksjoner
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
plot(thetas,risk_ML,type="l",ylim=c(0,0.17),bty="l",
     ylab="Risiko",xlab=expression(theta))
lines(thetas,risk_UP,col="blue")
lines(thetas,risk_JP,col="red")
legend("topright",c("ML","Uniform","Jeffreys"),bty="n",
       col=c("black","blue","red"),lty=c(1,1,1))

# Oppgave f) ------------------------------------------------------------------
# Beregning av a og b fra et empirisk bayes perspektiv.
ys = c(1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 
       3, 3, 5, 6, 6, 14, 15, 15, 20, 75)

llik = function(p,ys){
  a = p[1]
  b = p[2]
  sum(lbeta(a+1,b+ys-1) - lbeta(a,b))
}

ab = nlm(function(p) -llik(p,ys),p=c(1,1))$estimate
a  = ab[1]
b  = ab[2]

# Oppgave g) ------------------------------------------------------------------
# Simulering fra posterioren p(a,b|y).

N = 1000000
n = length(ys)

sampleit = function() {
  a = runif(1,0,100)
  b = runif(1,0,100)
  c(a,b,prod(beta(a+1,b+ys-1)/beta(a,b)))
}

samples = replicate(N,sampleit())
samples[3,] = samples[3,]/mean(samples[3,])

(expected_a = mean(samples[1,]*samples[3,]))
(expected_b = mean(samples[2,]*samples[3,]))

sd_a =  sqrt(mean(samples[1,]^2*samples[3,]) - expected_a^2)
sd_b =  sqrt(mean(samples[2,]^2*samples[3,]) - expected_b^2)

corr = (mean(samples[1,]*samples[2,]*samples[3,]) - 
          expected_a*expected_b)/(sd_a*sd_b)

# Oppgave h) ------------------------------------------------------------------
# Tabell over verdier.

parti = function(i,x) {
  mean(sapply(1:N,function(j) dbeta(x,samples[1,j]+1,
                    samples[2,j]+ys[i]-1)*samples[3,j]))
}

thetas = seq(0,1,by=0.001)

cols = function(i) {
  means  = sapply(thetas,function(x) parti(i,x))
  quants = cumsum(means)/length(means)
  q005   = thetas[which.min(abs(quants-0.05))]
  q050   = thetas[which.min(abs(quants-0.50))]
  q095   = thetas[which.min(abs(quants-0.95))]
  c(q005,q050,q095)
}

results = t(sapply(1:20,cols))
col1    = ys
col2    = 1/ys
col3    = (a+1)/(a+b+ys)
res     = cbind(col1,col2,col3,results)
round(res,3)