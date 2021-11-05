### Simulation d'un mouvement bownien standard
# Discretisation du temps
temps = seq(0,1,length=501)
pas.temps = 1/500
#simulation des accroissements
B.acc = rnorm(500, sd = sqrt(pas.temps))
B.sim = c(0, cumsum(B.acc))
plot(temps, B.sim, type = "l", xlab="Temps", main = "Mvt brownien" )


n.sim = 300
n.point = 201
temps = seq(0,1,length = n.point)
pas.temps = 1/(n.point-1)
B.acc = matrix(rnorm((n.point-1)*n.sim, sd = sqrt(pas.temps)), nrow=n.sim)
B.sim = matrix(NA, ncol=n.point,nrow=n.sim)
for(i in 1:n.sim){
  B.sim[i,] = c(0,cumsum(B.acc[i,]))
}
plot(temps, B.sim[1,], type = "l", xlab="Temps")
title("un echantillon de trajectoires")
for(i in 1:15){lines(temps, B.sim[i+1,])}



mu = 0.08;
sigma = 0.2;
T=1;
n=2^10
X0=0.1;
dt = T/n
t=seq(0,T, by = dt)
x = c(X0, mu*dt+sigma*sqrt(dt)*rnorm(n, mean = 0, sd = 1))
Xt = cumsum(x)
plot(t, Xt, type = "l", xlab = "time")



B = gbm(x0 = 1, mu = 1, sigma = 0.5, t0 = 0, t = 1, n = 1000)
plot(B)

start_time = Sys.time()

gbm = function(x0, mu, sigma, t0, t, n){
  delta = (t-t0)/n
  W = numeric(n+1)
  tseq = seq(t0,t,length = n+1)
  for(i in 2:(n+1))
    W[i] = W[i-1]+rnorm(1)*sqrt(delta)
  S = x0 * exp((mu-sigma^2)/2*(tseq-t0)+sigma*W)
  X = ts(S, start = t0, deltat = delta)
  return(X)
}

mu = 0.16
sigma = 0.2
P0 = 40
T = 1/12
t0=0
nt = 50
n=2^8

# generate n trajectoires
dt = T/n
t = seq(0,T, by = dt)
X = matrix(rep(0,length(t)*nt), nrow=nt)
for(i in 1:nt){X[i,] = gbm(x0=P0, mu = mu, sigma = sigma, t0=t0, t= T, n = n)}
##PLot
ymax=max(X)
ymin=min(X)
plot(t, X[1,], t = "l", ylim = c(ymin, ymax), col=1, ylab ="price P(t)", xlab = "temps")
for(i in 2:nt){lines(t, X[i,], t="l", ylim = c(ymin, ymax), col = i)}
mean(X[,n+1])

end_time = Sys.time()
end_time-start_time
