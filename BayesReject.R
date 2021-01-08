ys = c(1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 5, 6, 6, 14, 15, 15, 20, 75)
n  = length(ys)

as = seq(0,100)
b  = 66

plot(as,beta(as+n,b + sum(ys) - n)/beta(as,b))


# Importance sampling ---------------------------------------------------------

sampleit = function() {
  a = runif(1,0,100)
  b = runif(1,0,100)
  c(a,b,prod(beta(a+1,b+ys-1)/beta(a,b)))
}

N = 10
samples = replicate(N,sampleit()) 
samples[3,] = samples[3,]/mean(samples[3,])

# Oppgave 2h ------------------------------------------------------------------

parti = function(i,x) {
  mean(sapply(1:N,function(j) dbeta(x,samples[1,j]+1,samples[2,j]+ys[i]-1)*samples[3,j]))
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
col1 = ys
col2 = 1/ys
col3 = (a+1)/(a+b+ys)
res = cbind(col1,col2,col3,results)
round(res,3)