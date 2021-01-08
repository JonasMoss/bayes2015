dabd = function(a,alpha) {
  log(alpha)*exp((a-100)*log(alpha))
}

qabd = function(q,alpha) {
  log(q)/log(alpha) + 100
}

alpha = 1.1
qs = seq(0,1,by=0.001)
plot(qs,qabd(qs,alpha))
as = seq(0,100,by=0.1)
plot(as,dabd(as,alpha))


f = function(p) (qgamma(0.1,p[1],p[2])-1)^2 +
                 (qgamma(0.9,p[1],p[2])-sqrt(10))^2

nlm(f,c(4,4))$estimate
