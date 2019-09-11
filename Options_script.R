####TIME-Series Analisys####

###Download TS
library(quantmod)
library(tseries)
getSymbols('JPM', src = "yahoo", from = '2000-01-01',to='2019-08-23')
Stock<- JPM$JPM.Close #Stock
sma_50 <-TTR::SMA(Stock,n = 50)
sma_200 <-TTR::SMA(Stock,n = 200)
#Visualize Data
plot(Stock,main= "JP Morgan Stock Price",type='l')
lines(sma_50, col='red')
lines(sma_200, col='green')


chartSeries(Stock, name = "JP Morgan & Chase", TA= "addBBands()",theme = 'white')
zoomChart("2006-10::2008-10")


#Log Return
S_log_ret <- na.remove(diff(log(Stock))) #na.remove from tseries


mean_S_log <- round(mean(S_log_ret), 2)
sd_S_log <- round(sd(S_log_ret), 2)
plot(S_log_ret,type='l', ylab='Log Returns')


#####Change point detection####
require(sde)
library(tseries)
S <- get.hist.quote("JPM", start = "2000-01-01", end = "2019-08-23") 
chartSeries(S, name = "JP Morgan & Chase", TA= "addBBands()",theme = 'white')
S <- S$Close
cpoint(S)
addVLine = function(dtlist) plot(addTA(xts(rep(TRUE, NROW(dtlist)), dtlist), on = 1, col = "#FF0000"))
addVLine(cpoint(S)$tau0)
S <- as.numeric(S)
n <- length(S)
X <- log(S[-1]/S[-n])



plot(X, type = "l", main = "Log Returns")
abline(v = cpoint(X)$tau0, lty = 3, col = "red")









####POLISHED Dataset####

getSymbols('JPM', src = "yahoo", from = '2007-10-31',to='2019-08-23')
Stock<- JPM$JPM.Close #Stock
chartSeries(Stock, name = "JP Morgan & Chase", TA= "addBBands()",theme = 'white')

S <- get.hist.quote("JPM", start = "2007-10-31", end = "2019-08-23")
S <- na.remove(S$Close)
n <- length(S)
X <- na.omit(diff(log(S)))



S_log_ret <- na.remove(diff(log(Stock))) #na.remove from tseries

mean_S_log <- round(mean(S_log_ret), 8)
sd_S_log <- round(sd(S_log_ret), 8)
plot(S_log_ret,type='l', ylab='Log Returns')


# Plot the histogram along with a legend 
hist(S_log_ret, breaks = 100, prob=T, cex.main = 0.9)
abline(v = mean_S_log, lwd = 2, col='red')
legend("topright", cex = 0.8, border = NULL, bty = "n",
       paste("mean=", mean_S_log, "; sd=", sd_S_log)) ##ADD GREENLEGEND



x <- seq(-5 * sd_S_log, 5 * sd_S_log, length = nrow(S_log_ret))
lines(x, dnorm(x, mean_S_log, sd_S_log), col = "green", lwd = 2)
legend("topleft", legend=c("Gaussian distribution"),col=c('green'), lty=1:2, cex=0.8)

# Right Tail zoom with respect Normal dist
plot(density(S_log_ret), main='Return EDF - upper tail', xlim = c(0.1, 0.2),
     ylim=c(0,2));
curve(dnorm(x, mean=mean(S_log_ret),sd=sd(S_log_ret)), from = -0.3, to = 0.2, add=TRUE, col="green")
legend("topleft", legend=c("Gaussian right-tail"),col=c('green'), lty=1:2, cex=0.8)


####Test QQ plot Normal####
library(car)
qqPlot(coredata(S_log_ret), distribution = "norm", ylab = "JP Returns", envelope = FALSE)

####t distribution####
library(fGarch)
fit_T_dist <-fitdistr( S_log_ret, 't')
params.T = fit_T_dist$estimate
mu_T = as.numeric(params.T[1])
sigma_T = as.numeric(params.T[2] * sqrt(params.T[3] / (params.T[3] - 2)))
nu.T = params.T[3]
x <- seq(-5 * sd_S_log, 5 * sd_S_log, length = nrow(S_log_ret))
hist(S_log_ret, 80, freq = FALSE)
lines(x, dstd(x, mean = mu_T, sd = sigma_T, nu = nu.T),
      lwd = 2, lty = 2, col = 'blue')
y <- seq(-5 * sd_S_log, 5 * sd_S_log, length = nrow(S_log_ret))
lines(x=y, dnorm(x, mean_S_log, sd_S_log), col = "green", lwd = 2)
legend("topright", legend=c( "Gaussian distribution ","T-distribution"),col=c( "green","blue"), lty=1:2, cex=0.8)

qqPlot(coredata(S_log_ret), distribution = "t",df = 3, ylab = "JP Ret", envelope = FALSE)

####Parameters Estimation####

#Alfa & sigma hat
Delta <- 1/252 
alpha.hat<-mean(X,na.rm=TRUE)/Delta 
sigma.hat <- sqrt(var(X,na.rm=TRUE)/Delta) 
mu.hat <- alpha.hat + 0.5*sigma.hat^2
sigma.hat
mu.hat
alpha.hat

##MLE function##

library(MASS)
library("stats4")

LogL <- function(mu=1, sigma=1) { 
  A =suppressWarnings( dnorm(S_log_ret, mu, sigma))
  -sum(log(A))} 
 
fit1 <- mle(LogL, start = list(mu = 1, sigma = 1),method = "BFGS") 
summary(fit1)
confint(fit1)


mu_mle <-as.numeric(fit1@coef[1])
sigma_mle<- as.numeric(fi1t@coef[2])

###QUASI MLE###
library(yuima)

Delta_t <- 1/252

gBm <- setModel(drift="mu*x", diffusion="sigma*x",solve.var = "S")

mod <- setYuima(model=gBm, data=setData(S, delta=Delta_t))

mle1 <- qmle(mod, start = list(mu = 0.02, sigma = 0.05),
             lower = list(mu=-0.01, sigma=0.0), upper = list(mu=1, sigma=1),
             method = "L-BFGS-B")
summary(mle1)

mle2 <- qmle(mod, start = list(mu = 0.001, sigma = 0.02),
             lower = list(mu=0.00002, sigma=0.015), upper = list(mu=.001, sigma=0.028),
             method = "L-BFGS-B")#### Attempt 4
summary(mle2)


mle3 <- qmle(mod, start = list(mu = 0.01, sigma = 0.06),
             lower = list(mu=0, sigma=0), upper = list(mu=0.01, sigma=0.09),
             method = "L-BFGS-B")
summary(mle3)



##########################
##European Option Pricing: 
##########################
####BS CALL####

EU_call_bs = function(S, K, r, sigma,t=0, T){
  d1 = (log(S/K)+(r+((sigma)^2)/2)*(T-t))/(sigma*sqrt(T-t))
  d2 = d1-(sigma*sqrt(T-t))
  return((S*pnorm(d1))-(K*exp(-r*(T-t))*pnorm(d2)))
}

####BS-Call-PRICE-PLOT####

lb     = 0         # lower limit of x axis 
ub     = 170       # upper limit of x axis
tau <- 1
r = 0.0199

K <- 100
call1   = EU_call_bs(c(lb:ub), K, r, T=tau, sigma = sigma.hat)
payoff <- ifelse(c(lb:ub<100),0,c(lb:ub)-K)
plot(x = c(lb:ub), y = call1, main = "Black-Scholes price", xlab = "S", ylab = expression(paste("C(S,", tau, 
                                                                                                ")")), xlim = c(lb, ub), ylim = c(0, max(call1)), type = "l", col = "blue", 
     lwd = 2)
lines(payoff,col='red')
legend("topleft", legend=c("Call Price", "Payoff"),col=c("blue", "red"), lty=1:1, cex=0.8)



####BS_PUT####
EU_put_bs <- function(x,t=0,T,r,K,sigma){
  d2 <- (log(x/K) + (r - 0.5* sigma^2)*(T-t))/ (sigma* sqrt(T-t))
  d1 <- d2 + sigma* sqrt(T-t)
  P <- exp(-r * (T-t))*K*pnorm(-d2) - x*pnorm(-d1)
  return(P)
}
####Montecarlo Method CALL & PUT####
f <- function(x) max(0, x - K)
g<- function(x) max(0, K - x)
set.seed(123)

MCprice <- function(x,t=0,T=1,r,sigma,M=1000,type){
  h <- function(m) {
    u <- rnorm(m/2)
    tmp <- c(x * exp((r - 0.5*sigma^2)*(T-t)+ sigma*sqrt(T-t)*u),#dynamics of Stock price Brownian Motio
             x * exp((r - 0.5*sigma^2)*(T-t)+ sigma*sqrt(T-t)*(-u)))
    if (type=='call') {
      mean(sapply(tmp, function(xx) f(xx)))
    } else if (type == 'put') {
      mean(sapply(tmp, function(xx) g(xx)))
    }
    
  }
  p <- h(M)
  p * exp(-r * (T-t))
}
####Parameter BS##########
#Call at December 30, 2019
estimationDate <- as.Date(Sys.Date(), format="%Y-%m-%d") 
optionExpiryDate <- as.Date("2019-12-30", format="%Y-%m-%d")

T <- as.numeric((optionExpiryDate - estimationDate)/252) #to compute actual T-t

threeM_rates <- (getSymbols("DGS3MO", src = "FRED", auto.assign = FALSE))# get actual r from FRED
r <-as.numeric(threeM_rates[as.Date(end(threeM_rates), format="%Y-%m-%d")])/100

S <- as.numeric(Stock[as.Date(end(Stock), format="%Y-%m-%d")])
K <- 100
T<-.4801587
sigma <- sigma.hat
r <- 0.0199
set.seed(123)
mc_call <- MCprice(x = S, t = 0, T = T, r = r,sigma=as.numeric(sigma.hat), M = 1000, type ="call" )
mc_call
set.seed(123)
mc_put <- MCprice(x = S, t = 0, T = T, r = r,sigma=as.numeric(sigma.hat), M = 100000, type ="put" )
mc_put
P <- EU_put_bs(S,K=K,t=0,T=T,r=r,sigma=as.numeric(sigma))
C <- EU_call_bs(S = S, t = 0, T = T, r = r, K = K, sigma = as.numeric(sigma.hat) )
P
C
####Put-Call parity#### 
C - S + K* exp(-r*(T)) #this proves that put-call parity is kept


####FFT######
FFTcall.price <- function(phi, S0, K, r, T, alpha = 1, N = 2^12, eta = 0.25) {
  m <- r - log(phi(-(0+1i)))
  phi.tilde <- function(u) (phi(u) * exp((0+1i) * u * m))^T
  psi <- function(v) exp(-r * T) * phi.tilde((v - (alpha + 1) *(0+1i)))/(alpha^2 + alpha - v^2 + (0+1i) * (2 * alpha + 1) * v)
  lambda <- (2 * pi)/(N * eta)
  b <- 1/2 * N * lambda
  ku <- -b + lambda * (0:(N - 1))
  v <- eta * (0:(N - 1))
  tmp <- exp((0+1i) * b * v) * psi(v) * eta *(3 + (-1)^(1:N) - ((1:N) - 1 == 0))/3
  ft <- fft(tmp)
  res <- exp(-alpha * ku) * ft/pi
  inter <- spline(ku, Re(res), xout = log(K/S0)) 
  return(inter$y * S0)
}
phiBS <- function(u) exp((0+1i) * u * (mu - 0.5 * sigma^2) -
                           0.5 * sigma^2 * u^2)
mu=mu.hat
sigma =sigma.hat

FFTcall.price(phiBS, S0 = S, K = K, r = r, T = T)
#########Quality of approximation####
par(mar=c(4,4,1,1))
K.seq <- seq(100, 130, length=100)
exactP <- NULL
fftP <- NULL
for(K in K.seq){ 
  exactP  <- c(exactP , fOptions::GBSOption(TypeFlag = "c", S = S, X = K, Time = T, r = r, b = r, sigma = sigma.hat)@price )
  fftP  <- c(fftP , FFTcall.price(phiBS, S0=S, K=K, r=r, T=T) )
}
plot(K.seq, exactP-fftP, type="l",xlab="strike price K")
########MC -BS convergence PLOT######

# settings
s.0   = S    # current stock price
k     = 100   # strike price
tau   = T   # time to maturity in years
r     = r  # annual interest rate
sigma = sigma.hat  # volatility
set.seed(123)

# function for simulating stock price by geometric brownian motion
MCPath = function(s.0, tau, n.t, r, sigma, n.mc) {
  dt = tau/n.t
  t = seq(0, tau, by = dt)
  s.t = s.0 * exp((r - sigma^2/2) * t) * 
    exp(replicate(n.mc, c(0, cumsum(sigma * sqrt(dt) * rnorm(n.t)))))
  return(s.t)
}

# function for calculating MC estimate and confidence intervall of european option price
MCEurop = function(s.0, k, tau, n.t, r, sigma, n.mc) {
  
  s.t = MCPath(s.0, tau, n.t, r, sigma, n.mc)
  
  # payoff
  c = matrix(pmax(0, s.t[nrow(s.t), ] - k), 1, ncol(s.t))  # call
  p = matrix(pmax(0, k - s.t[nrow(s.t), ]), 1, ncol(s.t))  # put
  v = exp(-r * tau) * rbind(c, p)                          # discounted
  
  # MC estimate
  v.mc = rowSums(v)/ncol(s.t)
  
  # 95% confidence intervall for MC-estimate
  conf95 = 1.96 * 1/sqrt(ncol(s.t)) * sqrt(1/(ncol(s.t) - 1) * rowSums((v - v.mc)^2))
  
  return(cbind(v.mc, conf95))
}

# function for calculating the Black-Scholes price of European option
BSEurop = function(s0, k, tau, r, sigma) {
  d1 = (log(s0/k) + (r + (sigma^2/2)) * tau)/(sigma * sqrt(tau))
  d2 = d1 - (sigma * sqrt(tau))
  c  = s0 * pnorm(d1) - k * exp(-r * tau) * pnorm(d2)     # call
  p  = -s0 * pnorm(-d1) + k * exp(-r * tau) * pnorm(-d2)  # put
  return(c(c, p))
}

# main

j      = seq(2, 6, 0.5)
n.mc   = 10^j
v.call = matrix(0, length(n.mc), 2)
v.put  = matrix(0, length(n.mc), 2)

# calculate MC prices for several n.mc
for (i in 1:length(n.mc)) {
  
  mc = MCEurop(s.0, k, tau, n.t = 1, r, sigma, n.mc = n.mc[i])
  
  # MC estimate
  v.call[i, 1] = mc[1, 1]
  v.put[i, 1]  = mc[2, 1]
  
  # 95% confidence
  v.call[i, 2] = mc[1, 2]
  v.put[i, 2]  = mc[2, 2]
}

# calculate Black-Scholes prices
bs = BSEurop(s.0, k, tau, r, sigma)

# plotting

# up = value + 95%CI, dw = value - 95%CI
v.call.up = v.call[, 1] + v.call[, 2]
v.call.dw = v.call[, 1] - v.call[, 2]
v.put.up  = v.put[, 1] + v.put[, 2]
v.put.dw  = v.put[, 1] - v.put[, 2]

# y-axis limit for plot
ylimit = c(min(v.call.dw, v.put.dw), max(v.call.up, v.put.up))

par(mar = c(5, 5, 2, 1))

# plot MC values
plot(n.mc, v.call[, 1], 
     log     = "x", type = "p", pch = 16, ylim = ylimit, col = "black",
     main    = "Convergence of MC Simulation",
     ylab    = "Price of European option", 
     xlab    = "Number of MC samples", 
     cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.2, cex.sub = 1.2)
lines(n.mc, v.put[, 1], type = "p", pch = 16, col = "darkred")

# plot Black-Scholes values
lines(n.mc, matrix(bs[1], length(n.mc), 1)[, 1], 
      type = "l", lty = 2, lwd = 2, col = "black")
lines(n.mc, matrix(bs[2], length(n.mc), 1)[, 1], 
      type = "l", lty = 2, lwd = 2, col = "darkred")

# plot errorbars
arrows(n.mc, v.call.dw, n.mc, v.call.up, 
       code = 3, angle = 90, length = 0.1, lwd = 2, col = "black")
arrows(n.mc, v.put.dw, n.mc, v.put.up, 
       code = 3, angle = 90, length = 0.1, lwd = 2, col = "darkred")

# plot legend
legend(2*10^3, max(bs[1], bs[2]) - 0.25 * min(bs[1], bs[2]), 
       legend = c("MC-Estimation with 95% CI", "Black-Scholes-Model", "Call", "Put"), 
       lwd = 2, lty = c(0, 2, 0, 0), pch = c(16, NA, 15, 15), cex = 1.2, 
       col = c("black", "black", "black", "darkred"))
##################
####Greeks####

T=tau
K <- 100
##Delta


#tau = seq(tau_min, tau_max, by = (tau_max - tau_min)/(steps - 1))
#S   = seq(S_max, S_min, by = -(S_max - S_min)/(steps - 1))
Delta_BS <- function(tau, S, K, r, sig){
  d1 <- (log(S/K)+ (r+0.5*sig^2)*(tau))/(sqrt(tau)*sig)
  pnorm(d1)
}

fOptions::GBSCharacteristics(TypeFlag = 'c', S=S,X=K,Time=T,r=r,b=r,sigma=sigma.hat)

#Montecarlo Delta
library(foreach)
MCdelta <- function(x = 1, t = 0, T = 1, r = 1, sigma = 1, M = 1000, f){ 
  h <- function(m){
    u <- rnorm(M/2)
    tmp <- c(x * exp((r - 0.5 * sigma^2) * (T - t) + sigma *sqrt(T - t) * u),
             x * exp((r - 0.5 * sigma^2) * (T - t) + sigma *sqrt(T - t) * (-u)))
    g <- function(z) f(z) * (log(z/x) - (r - 0.5 * sigma^2) * (T - t))/(x * sigma^2 * (T - t)) 
    mean(sapply(tmp, function(z) g(z)))
  }
  nodes <- getDoParWorkers()
  p <- foreach(m = rep(M/nodes, nodes), .combine = "c") %dopar%
    h(m)
  p <- mean(p)
  p * exp(-r * (T - t))
}
set.seed(123)
MCdelta(x=S,T=T,r=r,sigma=sigma.hat,f=f,M=10000)

#The numerical approximation -  use the centered derivative
h <- 0.01

delta.num <- function(x){ (EU_call_bs(S = S+h, t = 0, T = T, r = r, K = K, sigma = as.numeric(sigma.hat))-
                             EU_call_bs(S = S-h, t = 0, T = T, r = r, K = K, sigma = as.numeric(sigma.hat) ))/( 2*h)}



delta.num(S)

MCdelta_mix <- function(x = 1, t = 0, T = 1, r = 1, sigma = 1, M = 1000, f, dx = 0.001){
  h <- function(m) {
    u <- rnorm(M/2)
    tmp1 <- c((x + dx) * exp((r - 0.5 * sigma^2) * (T - t) +sigma * sqrt(T - t) * u), 
              (x + dx) * exp((r - 0.5 * sigma^2) * (T - t) + sigma * sqrt(T - t) * (-u)))
    tmp2 <- c((x - dx) * exp((r - 0.5 * sigma^2) * (T - t) +sigma * sqrt(T - t) * u), 
              (x - dx) * exp((r - 0.5 * sigma^2) * (T - t) + sigma * sqrt(T - t)* (-u))) 
    mean(sapply(tmp1, function(x) f(x)) - sapply(tmp2, function(x) f(x)))/(2 * dx)
  }
  nodes <- getDoParWorkers()
  p <- foreach(m = rep(M/nodes, nodes), .combine = "c") %dopar%
    h(m)
  p <- mean(p)
  p * exp(-r * (T - t))
}
set.seed(123)
MCdelta_mix(x=S,T=T,r=r,sigma=sigma.hat,f=f)

#PLOT
# parameter settings for plots
S_min   = 50          # lower bound of Asset Price
S_max   = as.numeric(max(Stock))         # upper bound of Asset Price 
tau_min = 0.01        # lower bound of Time to Maturity
tau_max = 1           # upper bound of Time to Maturity
K       = 100         # exercise price
r       = 0.02         # riskfree interest rate                  
sig     = 0.3       # volatility               
d       = 0        # dividend rate                
steps   = 60          # steps 
q=0
tau = seq(tau_min, tau_max, by = (tau_max - tau_min)/(steps - 1))
S   = seq(S_max, S_min, by = -(S_max - S_min)/(steps - 1))
mesh = outer(tau, sort(S), Delta_BS, K = K, r = r, sig = sig)
title = bquote(expression(paste("Strike price is ", .(K), ", interest rate is ", 
                                .(r), ", annual volatility is ", .(sig))))



lattice::wireframe(mesh, drape = T, main = expression(paste("Delta as function of the time to maturity ", 
                                                            tau, " and the asset price S")), sub = title, scales = list(arrows = FALSE,col = "black", distance = 1, tick.number = 8, cex = 0.7,
                                                                                                                        x = list(labels = round(seq(tau_min,  tau_max, length = 7), 1)),
                                                                                                                        y = list(labels = round(seq(S_min, S_max, length = 7),  1))), 
                   xlab = list(expression(paste("Time to Maturity  ", tau)), rot = 30, cex = 1.2), 
                   ylab = list("Asset Price S", rot = -40, cex = 1.2), zlab = list("Delta",   cex = 1.1))


##Vega

Vega_BS <- function(S,sig,Tau,K,r) {
  d1 <- (log(S/K)+ (r+0.5*sig^2)*(Tau))/(sqrt(Tau)*sig)
  
  return(sqrt(Tau)*S*dnorm(d1))
}

MCvega<- function(x=1, t=0, T=1, r=1, sigma=1, M=1000, f, ds=0.01){ 
  u <- rnorm(M)
  h <- x*exp((r-0.5*(sigma+ds)^2)*(T-t)+(sigma+ds)*sqrt(T-t)*u)
  p <- sapply(h, function(x) f(x))
  p <- mean(p)
  p1 <- p*exp(-r*(T-t))
  h <- x*exp((r-0.5*(sigma-ds)^2)*(T-t)+(sigma-ds)*sqrt(T-t)*u) 
  p <- sapply(h, function(x) f(x))
  p <- mean(p)
  p2 <- p*exp(-r*(T-t))
  (p1-p2)/(2*ds)
}
#plot
meshgrid = function(a, b) {
  list(x = outer(b * 0, a, FUN = "+"), y = outer(b, a * 0, FUN = "+"))
}
first = meshgrid(seq(tau_min, tau_max, (tau_max - tau_min)/(steps - 1)), seq(tau_min, 
                                                                             tau_max, (tau_max - tau_min)/(steps - 1)))

tau    = first$x
dump   = first$y
second = meshgrid(seq(S_max, S_min, -(S_max - S_min)/(steps - 1)), seq(S_max, S_min, 
                                                                       -(S_max - S_min)/(steps - 1)))

dump2 = second$x
S     = second$y 
d1   = (log(S/K) + (r  + sig^2/2) * tau)/(sig * sqrt(tau))
Vega = S  * dnorm(d1) * sqrt(tau)
title = bquote(expression(paste("Strike price is ", .(K), ", interest rate is ", 
                                .(r), ", dividend rate is ", .(d), ", annual volatility is ", .(sig))))
lattice::wireframe(Vega ~ tau * S, drape = T, ticktype = "detailed", main = expression(paste("Vega as function of the time to maturity ", 
                                                                                             tau, " and the asset price S")), sub = title, scales = list(arrows = FALSE, 
                                                                                                                                                         col = "black", distance = 1, tick.number = 8, cex = 0.7, x = list(labels = round(seq(tau_min, 
                                                                                                                                                                                                                                              tau_max, length = 11), 1)), y = list(labels = round(seq(S_min, S_max, length = 11), 
                                                                                                                                                                                                                                                                                                  1))), xlab = list(expression(paste("Time to Maturity  ", tau)), rot = 30, 
                                                                                                                                                                                                                                                                                                                    cex = 1.2), ylab = list("Asset Price S", rot = -40, cex = 1.2), zlab = list("Vega", 
                                                                                                                                                                                                                                                                                                                                                                                                cex = 1.1))



##Gamma
Gamma_BS <- function(tau, S, K, r, sig) {
  d1 <- (log(S/K)+ (r+0.5*sig^2)*(tau))/(sqrt(tau)*sig)
  return(dnorm(d1)/(sig*S*sqrt(tau)))
}

tau = seq(tau_min, tau_max, by = (tau_max - tau_min)/(steps - 1))
S   = seq(S_max, S_min, by = -(S_max - S_min)/(steps - 1))

MCgamma <- function(x=1, t=0, T=1, r=1, sigma=1, M=1000, f, dS=0.01){ 
  u <- rnorm(M)
  z <- (x+dS)*exp((r-0.5*sigma^2)*(T-t)+sigma*sqrt(T-t)*u)
  g <- function(t,x,z)f(z) * (log(z/x)-(r-0.5*sigma^2)*(T-t))/(x*sigma^2*(T-t))
  p <- sapply(z, function(z) g(t,(x+dS),z) )
  p <- mean(p)
  p1 <- p*exp(-r*(T-t))
  z <- (x-dS)*exp((r-0.5*sigma^2)*(T-t)+sigma*sqrt(T-t)*u) 
  g <- function(t,x,z)
    f(z) * (log(z/x)-(r-0.5*sigma^2)*(T-t))/(x*sigma^2*(T-t)) 
  p <- sapply(z, function(z) g(t,(x-dS),z) )
  p <- mean(p)
  p2 <- p*exp(-r*(T-t))
  (p1-p2)/(2*dS)
}

gamma = function(tau, S, K, r, sig) {
  d1 = (log(S/K) + (r + sig^2/2) * tau)/(sig * sqrt(tau))
  return(dnorm(d1)/(S * (sig * sqrt(tau))))
}

mesh = outer(tau, sort(S), gamma, K = K, r = r, sig = sig)

lattice::wireframe(mesh, drape = T, main = expression(paste("Gamma as function of the time to maturity ", 
                                                            tau, " and the asset price S")), aspect = c(1, 0.75), sub = "Strike price is 100 and annual volatility is 0.30", 
                   scales = list(arrows = FALSE, col = "black", distance = 1, tick.number = 8, 
                                 cex = 0.7, x = list(labels = round(seq(tau_min, tau_max, length = 7), 1)), 
                                 y = list(labels = round(seq(S_min, S_max, length = 7), 1))), xlab = list(expression(paste("Time to Maturity  ", 
                                                                                                                           tau)), rot = 30, cex = 1.2), ylab = list("Asset Price S", rot = -30, cex = 1.2), 
                   zlab = list("Gamma", rot = 90, cex = 1.2), screen = list(z = 45, x = -70))

#Theta
Theta_BS <- function(S,K,r,sig,tau,Type){
  d1 <- (log(S/K)+ (r+0.5*sig^2)*(tau))/(sqrt(tau)*sig)
  d2 <- d1-(sig*sqrt(tau))
  if (Type=="call") {
    return((-r*K*exp(-r*(tau))*pnorm(d2)) - ((sig*S)/(2*sqrt(tau)))*dnorm(d1))
  } else if (Type =="put") {
    
    return((r*K*exp(-r*(tau))*pnorm(-d2)) - ((sigma*S)/(2*sqrt(tau)))*dnorm(d1))
  }
  
}

MCtheta <- function(x=1, t=0, T=1, r=1, sigma=1, M=1000, f,dt=0.01){ 
  u <- rnorm(M)
  h <- x*exp((r-0.5*sigma^2)*(T-(t+dt))+sigma*sqrt(T-(t+dt))*u)
  p <- sapply(h, function(x) f(x))
  p <- mean(p)
  p1 <- p*exp(-r*(T-(t+dt)))
  h <- x*exp((r-0.5*sigma^2)*(T-(t-dt))+sigma*sqrt(T-(t-dt))*u) 
  p <- sapply(h, function(x) f(x))
  p <- mean(p)
  p2 <- p*exp(-r*(T-(t-dt)))
  (p1-p2)/(2*dt)
}

meshgrid = function(a, b) {
  list(x = outer(b * 0, a, FUN = "+"), y = outer(b, a * 0, FUN = "+"))
}
first = meshgrid(seq(tau_min, tau_max, (tau_max - tau_min)/(steps - 1)), seq(tau_min, 
                                                                             tau_max, (tau_max - tau_min)/(steps - 1)))

tau    = first$x
dump   = first$y
second = meshgrid(seq(S_max, S_min, -(S_max - S_min)/(steps - 1)), seq(S_max, S_min, 
                                                                       -(S_max - S_min)/(steps - 1)))

dump2 = second$x
S     = second$y


d1    = (log(S/K) + (r + sig^2/2) * tau)/(sig * sqrt(tau))
d2 = d1-(sig*sqrt(tau))

theta = (-r*K*exp(-r*(tau))*pnorm(d2)) - ((sig*S)/(2*sqrt(tau)))*dnorm(d1)
# Plot
title = bquote(expression(paste("Strike price is ", .(K), ", interest rate is ", 
                                .(r),  ", annual volatility is ", .(sig))))
lattice::wireframe(theta ~ tau * S, drape = T, ticktype = "detailed", main = expression(paste("Theta as function of the time to maturity ", 
                                                                                              tau, " and the asset price S")), sub = title, scales = list(arrows = FALSE, 
                                                                                                                                                          col = "black", distance = 1, tick.number = 8, cex = 0.7, x = list(labels = round(seq(tau_min, 
                                                                                                                                                                                                                                               tau_max, length = 11), 1)), y = list(labels = round(seq(S_min, S_max, length = 11), 
                                                                                                                                                                                                                                                                                                   1))), xlab = list(expression(paste("Time to Maturity  ", tau)), rot = 30, 
                                                                                                                                                                                                                                                                                                                     cex = 1.2), ylab = list("Asset Price S", rot = -40, cex = 1.2), zlab = list("Theta", 
                                                                                                                                                                                                                                                                                                                                                                                                 cex = 1.1))





#Rho
Rho_BS <- function(x,Type,sigma,T,t=0){
  d1 <- (log(x/K)+ (r+0.5*sigma^2)*(T-t))/(sqrt(T-t)*sigma)
  d2 <- d1-(sigma*sqrt(T-t))
  if (Type=="call") {
    return((K*(T-t)*exp(-r*(T-t))*pnorm(d2))) 
  } else if (Type =="put") {
    return((-K*(T-t)*exp(-r*(T-t))*pnorm(-d2))) 
  }
  
}

rho =(K*(tau)*exp(-r*(tau))*pnorm(d2))
lattice::wireframe(rho ~ tau * S, drape = T, ticktype = "detailed", main = expression(paste("Rho as function of the time to maturity ", 
                                                                                            tau, " and the asset price S")), sub = title, scales = list(arrows = FALSE, 
                                                                                                                                                        col = "black", distance = 1, tick.number = 8, cex = 0.7, x = list(labels = round(seq(tau_min, 
                                                                                                                                                                                                                                             tau_max, length = 11), 1)), y = list(labels = round(seq(S_min, S_max, length = 11), 
                                                                                                                                                                                                                                                                                                 1))), xlab = list(expression(paste("Time to Maturity  ", tau)), rot = 30, 
                                                                                                                                                                                                                                                                                                                   cex = 1.2), ylab = list("Asset Price S", rot = -40, cex = 1.2), zlab = list("Rho", 
                                                                                                                                                                                                                                                                                                                                                                                               cex = 1.1))





#Volga
Volga_BS <- function(x,q=0,sigma,t=0,T){
  d1 <- (log(x/K)+ (r+0.5*sigma^2)*(T-t))/(sqrt(T-t)*sigma)
  d2 <- d1-(sigma*sqrt(T-t))
  
  return(exp(-q*t)*sqrt(T-t)*dnorm(d1)*((d1*d2)/(sigma)))
}



q       = 0           # dividend rate
steps   = 100         # steps

meshgrid = function(a, b) {
  list(x = outer(b * 0, a, FUN = "+"), y = outer(b, a * 0, FUN = "+"))
}

first = meshgrid(seq(tau_min, tau_max, -(tau_min - tau_max)/(steps - 1)), seq(tau_min, tau_max, -(tau_min - tau_max)/(steps - 
                                                                                                                        1)))

tau  = first$x
dump = first$y

second = meshgrid(seq(S_min, S_max, -(S_min - S_max)/(steps - 1)), seq(S_min, S_max, -(S_min - S_max)/(steps - 1)))

dump2 = second$x
S     = second$y

d1    = (log(S/K) + (r + sig^2/2) * tau)/(sig * sqrt(tau))
volga = (S * sqrt(tau) * exp(-q * tau) * dnorm(d1) * d1 * (d1 - sig * sqrt(tau)))

# Plot
title = bquote(expression(paste("Strike price is ", .(K), ", interest rate is ", 
                                .(r), ", annual volatility is ", .(sig))))
lattice::wireframe(volga ~ tau * S, drape = T, ticktype = "detailed", main = expression(paste("Volga as function of the time to maturity ", 
                                                                                              tau, " and the asset price S")), sub = title, scales = list(arrows = FALSE, 
                                                                                                                                                          col = "black", distance = 1, tick.number = 8, cex = 0.7, x = list(labels = round(seq(tau_min, 
                                                                                                                                                                                                                                               tau_max, length = 11), 1)), y = list(labels = round(seq(S_min, S_max, length = 11), 
                                                                                                                                                                                                                                                                                                   1))), xlab = list(expression(paste("Time to Maturity  ", tau)), rot = 30, 
                                                                                                                                                                                                                                                                                                                     cex = 1.2), ylab = list("Asset Price S", rot = -40, cex = 1.2), zlab = list("Volga", 
                                                                                                                                                                                                                                                                                                                                                                                                 cex = 1.1))

##Vanna
Vanna<- function(x,q=0,t=0,T,sigma){
  d1 <- (log(x/K)+ (r+0.5*sigma^2)*(T-t))/(sqrt(T-t)*sigma)
  d2 <- d1-(sigma*sqrt(T-t))
  
  return(exp(-q*t)*sqrt(T-t)*dnorm(d1)*((d2)/(sigma)))
}

exp(-q*t)*sqrt(T-t)*dnorm(d1)*((d1*d2)/(sigma))

q       = 0           # dividend rate
steps   = 100         # steps
b=r-q
meshgrid = function(a, b) {
  list(x = outer(b * 0, a, FUN = "+"), y = outer(b, a * 0, FUN = "+"))
}

first = meshgrid(seq(tau_min, tau_max, -(tau_min - tau_max)/(steps - 1)), seq(tau_min, tau_max, -(tau_min - tau_max)/(steps - 
                                                                                                                        1)))

tau  = first$x
dump = first$y

second = meshgrid(seq(S_min, S_max, -(S_min - S_max)/(steps - 1)), seq(S_min, S_max, -(S_min - S_max)/(steps - 1)))

dump2 = second$x
S     = second$y

d1    = (log(S/K) + (r - q - sig^2/2) * tau)/(sig * sqrt(tau))
d2    = d1 - sig * sqrt(tau)
Vanna = -(exp((b - r) * tau) * d2)/sig * dnorm(d1)

# Plot
title = bquote(expression(paste("Strike price is ", .(K), ", interest rate is ", 
                                .(r), ", dividend rate is ", .(q), ", annual volatility is ", .(sig))))
lattice::wireframe(Vanna ~ tau * S, drape = T, ticktype = "detailed", main = expression(paste("Vanna as function of the time to maturity ", 
                                                                                              tau, " and the asset price S")), sub = title, scales = list(arrows = FALSE, 
                                                                                                                                                          col = "black", distance = 1, tick.number = 8, cex = 0.7, x = list(labels = round(seq(tau_min, 
                                                                                                                                                                                                                                               tau_max, length = 11), 1)), y = list(labels = round(seq(S_min, S_max, length = 11), 
                                                                                                                                                                                                                                                                                                   1))), xlab = list(expression(paste("Time to Maturity  ", tau)), rot = 30, 
                                                                                                                                                                                                                                                                                                                     cex = 1.2), ylab = list("Asset Price S", rot = -40, cex = 1.2), zlab = list("Vanna", 
                                                                                                                                                                                                                                                                                                                                                                                                 cex = 1.1))






####GARCH####

library(moments)
acf(S_log_ret,type = 'correlation',main='Log Returns')
acf(S_log_ret^2,type = 'correlation',main='Squared Log Returns')
acf(abs(S_log_ret),type = 'correlation',main='|Log Returns|')



Box.test(abs(S_log_ret), lag = 10, type = c("Ljung-Box")) 
Box.test(S_log_ret^2, lag = 10, type = c("Ljung-Box"))

library(rugarch)
ret.jpm <- dailyReturn(Cl(JPM), type='log')

garch11.spec = ugarchspec(variance.model = list(model="sGARCH",garchOrder=c(1,1)), mean.model = list(armaOrder=c(0,0)))
jpm.garch11.fit = ugarchfit(spec=garch11.spec, data=ret.jpm)
show(jpm.garch11.fit )
#plot(jpm.garch11.fit)
jpm.garch11.fit


egarch11.spec = ugarchspec(variance.model = list(model="eGARCH",garchOrder=c(1,1)), mean.model = list(armaOrder=c(0,0)))
jpm.egarch11.fit = ugarchfit(spec=egarch11.spec, data=ret.jpm)
ni.egarch11 <- newsimpact(jpm.egarch11.fit)
jpm.egarch11.fit
plot(ni.egarch11$zx, ni.egarch11$zy, type="l", lwd=2, col="blue",
     main="EGARCH(1,1) - News Impact",
     ylab=ni.egarch11$yexpr, xlab=ni.egarch11$xexpr)


tgarch11.spec = ugarchspec(variance.model = list(model="fGARCH",
                                                 submodel="TGARCH", garchOrder=c(1,1)),
                           mean.model = list(armaOrder=c(0,0)))
jpm.tgarch11.fit = ugarchfit(spec=tgarch11.spec, data=ret.jpm)
ni.tgarch11 <- newsimpact(jpm.tgarch11.fit)
jpm.tgarch11.fit
plot(ni.tgarch11$zx, ni.tgarch11$zy, type="l", lwd=2, col="blue",
     main="TGARCH(1,1) - News Impact",
     ylab=ni.tgarch11$yexpr, xlab=ni.tgarch11$xexpr)

#plot(jpm.tgarch11.fit)

jpm.garch11.fit = ugarchfit(spec=garch11.spec, data=ret.jpm, out.sample=20)
jpm.garch11.fcst = ugarchforecast(jpm.garch11.fit, n.ahead=10,n.roll=10)
#plot(jpm.garch11.fcst)

####Levy#####
##Retunrs & ECDF
 stock.rtn=na.remove(diff(log(JPM$JPM.Close)))
 returns <- as.vector(stock.rtn)
 m=mean(returns,na.rm=TRUE)
 s=sd(returns,na.rm=TRUE)
 times=index(stock.rtn)
 n = sum(is.na(returns))+sum(!is.na(returns))
 x=seq(1,n)
 y=rnorm(n, m, s)
 plot(times,returns,pch=19,xaxs="i",cex=0.03,col="blue", ylab="X", xlab="n", main = '')
 segments(x0 = times, x1 = times, y0 = 0, y1 = returns,col="blue")
 points(times,y,pch=19,cex=0.3,col="red", ylab="X", xlab="n", main = '')
 abline(h = 3*s, col="black", lwd =1)
 abline(h = -3*s, col="black", lwd =1)

 stock.ecdf=ecdf(as.vector(stock.rtn))
 x <- seq(-0.25, 0.25, length=100)
 px <- pnorm((x-m)/s)
 plot(stock.ecdf, xlab = 'Sample Quantiles', col="blue",ylab = '', main = 'Empirical Cumulative
 Distribution')
 lines(x, px, type="l", lty=2, col="red",xlab="x value",ylab="Density", main="Gaussian cdf")
 legend("topleft", legend=c("Empirical cdf", "Gaussian cdf"),col=c("blue", "red"), lty=1:2, cex=0.8)

##FITTING LEVY
par(mfrow=c(2, 2)) 
par(mar=c(3, 3, 3, 1)) 
grid <- NULL 
library(fBasics)
nFit(X)
nigFit(X, trace=FALSE) 
hypFit(X, trace=FALSE)
ghFit(X, trace=FALSE)

#FFT method for VG
FFTcall.price <- function(phi, S0, K, r, T, alpha = 1, N = 2^12, eta = 0.25) {
   m <- r - log(phi(-(0+1i)))
   phi.tilde <- function(u) (phi(u) * exp((0+1i) * u * m))^T
   psi <- function(v) exp(-r * T) * phi.tilde((v - (alpha + 
                                                      1) * (0+1i)))/(alpha^2 + alpha - v^2 + (0+1i) * (2 * 
                                                                                                         alpha + 1) * v)
   lambda <- (2 * pi)/(N * eta)
   b <- 1/2 * N * lambda
   ku <- -b + lambda * (0:(N - 1))
   v <- eta * (0:(N - 1))
   tmp <- exp((0+1i) * b * v) * psi(v) * eta * (3 + (-1)^(1:N) - 
                                                  ((1:N) - 1 == 0))/3
   ft <- fft(tmp)
   res <- exp(-alpha * ku) * ft/pi
   inter <- spline(ku, Re(res), xout = log(K/S0))
   return(inter$y * S0)
 }



#Variance-Gamma method
#Parameters fitting
library(VarianceGamma)

vg_param <- vgFit(S_log_ret)$param #esitmate VG parameters on the sample


c <- as.numeric(vg_param[1])
sigma <- as.numeric(vg_param[2])
theta <- as.numeric(vg_param[3])
nu <- as.numeric(vg_param[4])

phiVG <- function(u) {
  omega <- (1/nu) * (log(1 - theta * nu - sigma^2 * nu/2))
  tmp<-1-(0+1i)*theta*nu*u+0.5*sigma^2*u^2*nu
  tmp <- tmp^(-1/nu)
  exp((0+1i)*u*log(S)+u*(r+omega)*(0+1i))*tmp }
#setting parameters
S <- as.numeric(Stock[as.Date(end(Stock), format="%Y-%m-%d")])
K <- 100
T<-.4801587

r <- 0.0199

price_nel_pa <-FFTcall.price(phiVG, S0 = S, K = K, r = r, T = T)
price_nel_pa
#MC VG 
MCpriceVG <- function(x,t=0,T=1,r,sigma,M=1000,type){
  h <- function(m) {
    u <- rnorm(m/2)
    t <- rgamma(n, shape = T/nu, scale = nu)
    N <- rnorm(n, 0, 1)
    X <- theta * t + N * sigma * sqrt(t)
    omega <- (1/nu) * (log(1 - theta * nu - sigma^2 * nu/2)) 
    tmp <- c(x * exp(r * T + omega * T + X),#dynamics of Stock price Brownian Motio
             x * exp(r * T + omega * T + X))
    if (type=='call') {
      mean(sapply(tmp, function(xx) f(xx)))
    } else if (type == 'put') {
      mean(sapply(tmp, function(xx) g(xx)))
    }
    
  }
  p <- h(M)
  p * exp(-r * (T-t))
} 
f <- function(x) max(0, x - K)
g<- function(x) max(0, K - x)
set.seed(123)
MC_VG_price <- MCpriceVG(x=S,T=T,sigma=sigma,type='call',r=r,M=n) 
MC_VG_price


set.seed(123)
#Quasi-Newton BFGS method
vg.QN <-vgFit(S_log_ret, method = "BFGS", hessian = TRUE)


vg_parm_QN <- vgFit(S_log_ret, method = "BFGS", hessian = TRUE)$param
c <- as.numeric(vg_parm_QN[1])
sigma <- as.numeric(vg_parm_QN[2])
theta <- as.numeric(vg_parm_QN[3])
nu <- as.numeric(vg_parm_QN[4])

price_QN_par <-FFTcall.price(phiVG, S0 = S, K = K, r = r, T = T)
price_QN_par
################# American options: Explicit finite-difference method####


AmericanPutExp <- function(Smin=0, Smax,  T=1, N=10, M=10, K, r=0.05, sigma=0.01){
  Dt = T/N 
  DS = (Smax-Smin)/M
  t <- seq(0, T, by =Dt) 
  S <- seq(Smin, Smax, by=DS)
  A <- function(j) (-0.5*r*j*Dt + 0.5*sigma^2*j^2*Dt)/(1+r*Dt) 
  B <- function(j) (1-sigma^2*j^2*Dt)/(1+r*Dt) 
  C <- function(j) (0.5*r*j*Dt + 0.5*sigma^2*j^2*Dt)/(1+r*Dt)
  P <- matrix(, M+1, N+1)
  colnames(P) <- round(t,2)
  rownames(P) <- round(rev(S),2)
  P[M+1, ] <- K   # C(,j=0) = K
  P[1,] <- 0   # C(,j=M) = 0
  P[,N+1] <- sapply(rev(S), function(x) max(K-x,0))
  optTime <- matrix(FALSE, M+1, N+1)
  optTime[M+1,] <- TRUE
  optTime[which(P[,N+1]>0),N+1] <- TRUE
  
  for(i in (N-1):0){
    for(j in 1:(M-1)){
      J <- M+1-j
      I <- i+1
      P[J,I] <- A(j)*P[J+1,I+1] + B(j)*P[J,I+1] + C(j)*P[J-1,I+1]
      if(P[J,I] < P[J,N+1])
        optTime[J,I] <- TRUE
    }
  }
  colnames(optTime) <- colnames(P)
  rownames(optTime) <- rownames(P)
  ans <- list(P=P, t=t, S=S, optTime=optTime,N=N,M=M)
  class(ans) <- "AmericanPut"
  return(invisible(ans))
}


### plot method for American put

plot.AmericanPut <- function( obj ){
  plot(range(obj$t),range(obj$S),type="n",axes=F,xlab="t", ylab="S")
  axis(1,obj$t,obj$t)
  axis(2,obj$S,obj$S)
  abline(v = obj$t, h = obj$S, col = "darkgray", lty = "dotted")
  for(i in 0:obj$N){
    for(j in 0:obj$M){
      J <- obj$M+1-j
      I <- i+1
      cl <- "grey"; 
      if(obj$optTime[J,I])
        cl <- "black"
      text(obj$t[i+1],obj$S[j+1], round(obj$P[J,I],2),cex=0.75, col=cl)
    }
  }
  DS <- mean(obj$S[1:2])
  y <- as.numeric(apply(obj$optTime,2, function(x) which(x)[1]))
  lines(obj$t, obj$S[obj$M+2-y]+DS, lty=2, col='red')
}



put <- AmericanPutExp(Smax = 65, sigma = sigma.hat, K = 60, T=1, Smin = 0, r=r,N=10,M=10)
round(put$P,2)

par(mar=c(3,3,1,1))
plot(put)



S0 <- 36
myval <- round(put$P[which(rownames(put$P)==S0),1],2)


### instability

put.bad <- AmericanPutExp(Smax = 65, sigma = sigma.hat, K = 60, M=15)
round(put.bad$P,2)
plot(put.bad)

####################### American options: Implicit finite-difference method###################
AmericanPutImp <- function( Smin=0, Smax,  T=1, N=10, M=10, K, r=0.05, sigma=0.01){
  Dt = T/N 
  DS = (Smax-Smin)/M
  t <- seq(0, T, by =Dt) 
  S <- seq(Smin, Smax, by=DS)
  
  A <- function(j) 0.5*r*j*Dt - 0.5*sigma^2*j^2*Dt 
  B <- function(j) 1+sigma^2*j^2*Dt+r*Dt
  C <- function(j) -0.5*r*j*Dt - 0.5*sigma^2*j^2*Dt
  
  a <- sapply(0:M, A)
  b <- sapply(0:M, B)
  c <- sapply(0:M, C)
  
  P <- matrix(, M+1, N+1)
  colnames(P) <- round(t,2)
  rownames(P) <- round(rev(S),2)
  
  P[M+1, ] <- K   # C(,j=0) = K
  P[1,] <- 0   # C(,j=M) = 0
  P[,N+1] <- sapply(rev(S), function(x) max(K-x,0))
  
  AA <- matrix(0, M-1, M-1)
  for(j in 1:(M-1)){
    if(j>1) AA[j,j-1] <- A(j)
    if(j<M) AA[j,j] <- B(j)
    if(j<M-1) AA[j,j+1] <- C(j)
  }
  
  optTime <- matrix(FALSE, M+1, N+1)
  for(i in (N-1):0){
    I <- i+1
    bb <- P[M:2,I+1]
    bb[1] <- bb[1]-A(1)*P[M+1-0,I+1]
    bb[M-1] <- bb[M-1]-C(M-1)*P[M+1-M,I+1] 
    P[M:2,I] <- solve(AA,bb)
    idx <- which(P[,I] < P[,N+1])
    P[idx,I] <- P[idx,N+1] 
    optTime[idx, I] <- TRUE
  }
  optTime[M+1,] <- TRUE 
  optTime[which(P[,N+1]>0),N+1] <- TRUE
  colnames(optTime) <- colnames(P)
  rownames(optTime) <- rownames(P)
  ans <- list(P=P, t=t, S=S, optTime=optTime,N=N,M=M)
  class(ans) <- "AmericanPut"
  return(invisible(ans))
}



put <- AmericanPutImp(Smax = 65, sigma = sigma.hat, K = 60, T=1, Smin = 0, r=r,N=10,M=10)
round(put$P,2)
putputbad <- AmericanPutImp(Smax = 65, sigma = sigma.hat, K = 60, M=15)

par(mar=c(3,3,1,1))
plot(put)
plot(putputbad)



###### Broadie and Glasserman Monte Carlo method#####

simTree <- function(b,d, S0, sigma, T, r){
  tot <- sum(b^(1:(d-1)))
  S <- numeric(tot+1) 
  S[1] <- S0
  dt <- T/d
  for(i in 0:(tot - b^(d-1))){
    for(j in 1:b){
      S[i*b+j +1] <- S[i+1]*exp((r - 0.5*sigma^2)*dt + sigma*sqrt(dt)*rnorm(1))
    }
  }
  S
}



upperBG <- function(S, b, d, f){
  tot <- sum(b^(1:(d-1)))
  start <- tot - b^(d-1) +1
  end <- tot +1
  P <- S
  P[start:end] <- f(S[start:end])
  tot1 <- sum(b^(1:(d-2)))
  for(i in tot1:0){
    m <- mean(P[i*b+1:b+1])
    v <- f(S[i+1])
    P[i+1] <- max(v,m)
  }
  P
}

lowerBG <- function(S, b, d, f){
  tot <- sum(b^(1:(d-1)))
  start <- tot - b^(d-1) +1
  end <- tot +1
  p <- S 
  p[start:end] <- f(S[start:end])
  tot1 <- sum(b^(1:(d-2)))
  
  m <- numeric(b)
  for(i in tot1:0){
    v <- f(S[i+1])
    for(j in 1:b){
      m[j] <- mean(p[i*b+(1:b)[-j]+1])
      m[j] <- ifelse( v>m[j], v, p[i*b+(1:b)[j]+1])
    }
    p[i+1] <- mean(m)
  }
  p
}

#setting parameters
set.seed(123)
b <- 3
d <- 3
K <- 100
f <- function(x) sapply(x, function(x) max(x-K,0))
g <- function(x) sapply(x, function(x) max(K-x,0))
T <- 1
r <- 0.0199
sigma <- sigma.hat
S0 <- 100

low <- 0
upp <- 0
M <- 1e5
for(i in 1:M){
  S <- simTree(b,d, S0, sigma, T, r)
  low <- low + lowerBG(S, b,d,f)[1]
  upp <- upp + upperBG(S, b,d,f)[1]
}
low/M
upp/M


#################Longstaff and Schwartz least squares method#####

LSM <- function(n, d, S0, K, sigma, r, T){
  s0 <- S0/K
  dt <- T/d
  z <- rnorm(n)
  s.t <- s0*exp((r-1/2*sigma^2)*T+sigma*z*(T^0.5))
  s.t[(n+1):(2*n)] <- s0*exp((r-1/2*sigma^2)*T-sigma*z*(T^0.5))
  CC <- pmax(1-s.t, 0)
  payoffeu <- exp(-r*T)*(CC[1:n]+CC[(n+1):(2*n)])/2*K
  euprice <- mean(payoffeu)
  
  for(k in (d-1):1){
    z <- rnorm(n)
    mean <- (log(s0) + k*log(s.t[1:n]))/(k+1)
    vol <- (k*dt/(k+1))^0.5*z
    s.t.1 <- exp(mean+sigma*vol)
    mean <- (log(s0) + k*log( s.t[(n+1):(2*n)] )) / ( k + 1 )
    s.t.1[(n+1):(2*n)] <- exp(mean-sigma*vol)
    CE <- pmax(1-s.t.1,0)
    idx<-(1:(2*n))[CE>0]
    discountedCC<- CC[idx]*exp(-r*dt)
    basis1 <- exp(-s.t.1[idx]/2)
    basis2 <- basis1*(1-s.t.1[idx])
    basis3 <- basis1*(1-2*s.t.1[idx]+(s.t.1[idx]^2)/2)
    
    p <- lm(discountedCC ~ basis1+basis2+basis3)$coefficients
    estimatedCC <- p[1]+p[2]*basis1+p[3]*basis2+p[4]*basis3
    EF <- rep(0, 2*n)
    EF[idx] <- (CE[idx]>estimatedCC)
    CC <- (EF == 0)*CC*exp(-r*dt)+(EF == 1)*CE
    s.t <- s.t.1
  }
  
  payoff <- exp(-r*dt)*(CC[1:n]+CC[(n+1):(2*n)])/2
  usprice <- mean(payoff*K)
  error <- 1.96*sd(payoff*K)/sqrt(n)
  earlyex <- usprice-euprice
  data.frame(usprice, error, euprice)
}
set.seed(123)
LSM(100, 3, S0, K, sigma, r, T)
set.seed(123)
LSM(1000, 3, S0, K, sigma, r, T)
set.seed(123)
LSM(1e5, 3, S0, K, sigma, r, T)




#####gBM & BASKET####

GBM <- function(N, sigma, mu, S0, Wt = NULL) {  
  if (is.null(Wt)) {
    Wt <- cumsum(rnorm(N, 0, 1))
  }
  t <- (1:N)/252
  p1 <- (mu - 0.5*(sigma*sigma)) * t
  p2 <- sigma * Wt
  St = S0 * exp(p1 + p2)
  return(St)
}

CorrelatedGBM <- function(N, S0, mu, sigma, cor.mat) {
  mu <- as.matrix(mu)
  sigma <- as.matrix(sigma)
  GBMs <- matrix(nrow = N, ncol = nrow(mu))
  Wt <- matrix(rnorm(N * nrow(mu), 0, 1), ncol = nrow(mu))
  Wt <- apply(Wt, 2, cumsum)
  chol.mat <- chol(cor.mat) # upper triangular cholesky decomposition
  Wt <- Wt %*% chol.mat   # key trick for creating correlated paths
  for (i in 1:nrow(mu)) {
    GBMs[,i] <- GBM(N, sigma[i], mu[i] , S0[i], Wt[, i])
  }
  return(GBMs)
}

GetPrices <- function(tickers, startDate, endDate) {
  prices <- get.hist.quote(instrument = tickers[1], start = startDate,end=endDate, 
                           quote = 'AdjClose')
  # download the rest of the prices
  for (tik in 2:length(tickers)) {
    tmp <- get.hist.quote(instrument = tickers[tik], 
                          start = start, end=end, quote = 'AdjClose')
    prices <- merge(prices, tmp)
  }
  return(prices)
}

set.seed(123)
N <- 4 * 252 # 4 years, each with 252 trading days
t <- (1:N)/252
start <- '2015-08-23'
end = '2019-08-23'
tickers <- c('JPM', 'MSFT')
prices <- GetPrices(tickers, start,end)

# get the cc returns and vectors of average returns and volaitiliy
returns.mat <- as.matrix(na.omit(diff(log(prices))))
mean.vec <- as.numeric(colMeans(returns.mat))
sigma.vec <- as.numeric(sqrt(apply(returns.mat, 2, var)))
prices.vec <- as.numeric(prices[nrow(prices)])
cor.mat <-as.matrix(cor(returns.mat))

paths <- CorrelatedGBM(N, prices.vec , mean.vec, sigma.vec, cor.mat)

colors <- c('darkblue', 'darkgreen', 'darkgoldenrod')
plot(t, paths[,1], type = 'l', ylim = c(50, max(paths)), xlab = "Year", 
     ylab = "Price", main = "Simulated Asset Prices", col = colors[1])
for (i in 2:ncol(paths)) {
  lines(t, paths[, i], col = colors[i])
}
legend(x = 0., y = 80, c('JPM', 'MSFT'), lty = c(1,1,1), col = colors, cex = 0.9)


BSMulti <- function( N, T=1, S0=1, r=0.01, z,f, M=1000){
  require(tcltk)
  X <- numeric(M)
  Z <- numeric(M)
  pb <- tkProgressBar(min = 1, max = M)
  for(i in 1:M){
    X[i] <- z( CorrelatedGBM(N, prices.vec , mean.vec, sigma.vec, cor.mat))
    setTkProgressBar(pb, i) 
  }
  for(i in 1:M){ 
    Z[i] <- f(X[i])
  }
  close(pb)
  exp(-r*T)*mean(Z)
}
z <- function(x) sapply(x, function(x) max(x))
f <- function(x) sapply(x, function(x) max(x-K,0))

set.seed(123)
K<-100
r=0.0199
T=1
BSMulti(S0=last(prices),z=z,f=f,N=N,M=1000,r=r)


###############COPULA & BASKET##############
require(copula)
require(quantmod)
getSymbols("JPM", from="2015-01-01", to="2019-08-23")
attr(JPM, "src")
S <- JPM [,"JPM.Close"]
X <- na.omit(diff(log(S)))
A <- as.numeric(X)

getSymbols("MSFT", from="2015-01-01", to="2019-08-23")
attr(MSFT, "src")
Z <- MSFT [,"MSFT.Close"]
W <- na.omit(diff(log(Z)))
B <- as.numeric(W)

dataset = cbind(A,B)
dataset = as.data.frame(dataset)
prices_ <- cbind(S,Z)

#estimate Kendall tau from data, Kendall Tau measure the association of two time series. 
tau <- cor.test(x=A, y=B, method = 'kendall')$estimate

#Generator as a function of t
t <- seq(0, 2, length.out = 257) # evaluation points

psi. <- cbind(Pi = exp(-t), # Pi generator
              C  = copClayton@psi(t, theta = iTau(claytonCopula(), tau)),
              F  = copFrank@psi  (t, theta = iTau(frankCopula(),   tau)),
              GH = copGumbel@psi (t, theta = iTau(gumbelCopula(),  tau))
)
plot(t, type = "l", lwd = 2,
     xlim = range(t), ylim = range(psi.), col = 1, ylab = "",
     xlab = quote(psi(t)~"as a function of t"))
for(j in 2:ncol(psi.)) lines(t, psi.[,j], col = j, lwd = 2)
legend("topright", bty = "n", lty = 1, lwd = 2, col = 2:ncol(psi.),
       legend = c( "Clayton", "Frank",
                   "Gumbel"))

#From Tau to Theta manually
gumbel.tau2theta  = function(tau){
  1/(1 - tau)
}

clayton.tau2theta = function(tau){
  2 * tau / (1-tau)
}

normal.tau2theta  = function(tau){
  sin(tau * pi/2)
}


means = colMeans(dataset)
sds   = sapply(dataset, sd)
dataset2 <- ordered(dataset$A)


params.margins = list(list(mean = means[[1]], sd = sds[[1]]), 
                      list(mean = means[[2]], sd = sds[[2]]))
# marginal distributions
marg = rep("norm", length(params.margins)) 

gumbel.cop = gumbelCopula(gumbel.tau2theta(tau))                             
out.gumbel        = mvdc(gumbel.cop, marg, params.margins)

normal.cop = normalCopula(normal.tau2theta(tau))
out.normal        = mvdc(normal.cop, marg, params.margins)

clayton.cop = claytonCopula(clayton.tau2theta(tau)) 
out.clayton        = mvdc(clayton.cop, marg, params.margins)


frank.cop =frankCopula(iTau(frankCopula(),tau))
out.frank        = mvdc(frank.cop, marg, params.margins)

# Create perspective plots of densities
persp(frank.cop, dCopula, phi = 20, theta = 20, ticktype = "detailed", ylab = "", 
      xlab = "Frank Copula", zlab = "", shade = 0.1)
persp(gumbel.cop, dCopula, phi = 20, theta = 20, ticktype = "detailed", ylab = "", 
      xlab = "Gumbel Copula", zlab = "", shade = 0.1)
persp(clayton.cop, dCopula, phi = 20, theta = 20, ticktype = "detailed", ylab = "", 
      xlab = "Clayton Copula", zlab = "", shade = 0.1)

### rClayton MO algorithm ################################################################


rClaytonMO <- function(u, theta)
{
  if(!is.matrix(u)) u <- rbind(u)
  
  dim. <- dim(u)
  n <- dim.[1]
  d <- dim.[2] - 1
  V <- qgamma(u[,d+1], shape=1/theta) # n-vector, frailty component
  E <- qexp(u[,seq_len(d)]) # (n,d)-matrix
  
  (1 + E/matrix(rep(V, d), ncol=d))^(-1/theta) # (n,d)-matrix
}

#####################rGBM################################
rGeoBM <- function(u, S0, mu, sigma, T)
{
  stopifnot(0 < u, u < 1, length(mu) == (d <- length(S0)), mu >= 0,
            length(sigma) == d, sigma >= 0, T > 0)
  log.diff <- qnorm(u) * matrix(rep(sigma, each=n), ncol=d) # (n,d)-matrix; or t(t(qnorm(u))*sigma)
  log.drft <- (mu - sigma^2 / 2) * T # d-vector
  log. <- matrix(rep(log.drft, n), ncol=d, byrow=TRUE) + log.diff # (n,d)-matrix
  matrix(rep(S0, n), ncol=d, byrow=TRUE) * exp(log.) # S_t, t in 1,..,T; (n,d)-matrix
}

##############PAYOFF#######################################

payoff <- function(K, N, S0, S, type = c("call", "put"),
                   method = c("basket", "worst.of", "best.of"))
{
  stopifnot(K >= 0, N >= 0, S0 >= 0, S >= 0, length(S0) == ncol(S))
  type <- match.arg(type)
  method <- match.arg(method)
  perf <- switch(method,
                 "basket" = {
                   rowMeans(t(t(S)/S0))
                 },
                 "worst.of" = {
                   apply(t(t(S)/S0), 1, min)
                 },
                 "best.of" = {
                   apply(t(t(S)/S0), 1, max)
                 },
                 stop("Wrong 'method'"))
  N * pmax(0, if(type=="call") perf - K else K - perf)
}


### Define parameters ######################################################

n <- 1e5 # Monte Carlo sample size
d <- 2 # dimension

## Stochastic process parameters
sigma <- as.numeric(sds) # volatilities
r <- 0.0199 # continuously compounded short rate
S0 <- as.numeric(last(prices_)) # initial stocks' levels
K <- 1 # option strike
N <- 1000 # option notional
T <- 1 # time horizon

## Copulas
## Clayton
family.C <- "Clayton"
th.C <- iTau(getAcop(family.C), tau) # corresponding parameter
clayton.cop <- onacopulaL(family.C, nacList=list(th.C, 1:d)) # Clayton copula
## t_3
family.t <- "t"
nu <- 3 # degrees of freedom
th.t <- iTau(ellipCopula(family.t, df=nu), tau) # corresponding parameter
t.cop <- ellipCopula(family.t, param=th.t, dim=d, df=nu) # define copula object


### Sampling ###############################################################
library(qrng)
## Uniform samples for CDM
set.seed(271)
U.CDM  <- matrix(runif(n*d), ncol=d) # pseudo
set.seed(271)
U.CDM. <- ghalton(n, d=d) # quasi

## Uniform samples for MO
set.seed(271)
U.MO  <- matrix(runif(n*(d+1)), ncol=d+1) # pseudo
set.seed(271)
U.MO. <- ghalton(n, d=d+1) # quasi

## t samples via CDM
U.t.CDM  <- cCopula(U.CDM,  cop=t.cop, inverse=TRUE) # pseudo
U.t.CDM. <- cCopula(U.CDM., cop=t.cop, inverse=TRUE) # quasi

## Clayton samples via CDM
U.C.CDM  <- cCopula(U.CDM,  cop=clayton.cop, inverse=TRUE) # pseudo
U.C.CDM. <- cCopula(U.CDM., cop=clayton.cop, inverse=TRUE) # quasi

## Clayton samples via MO
U.C.MO  <- rClaytonMO(U.MO,  theta=th.C) # pseudo
U.C.MO. <- rClaytonMO(U.MO., theta=th.C) # quasi


## Geometric Brownian Motion samples
S.t.CDM <- rGeoBM(U.t.CDM, S0=S0, mu=rep(r, d), sigma=sigma, T=T)
S.C.CDM <- rGeoBM(U.C.CDM, S0=S0, mu=rep(r, d), sigma=sigma, T=T)
S.C.MO  <- rGeoBM(U.C.MO,  S0=S0, mu=rep(r, d), sigma=sigma, T=T)

## Quasi-geometric Brownian Motion samples
S.t.CDM. <- rGeoBM(U.t.CDM., S0=S0, mu=rep(r, d), sigma=sigma, T=T)
S.C.CDM. <- rGeoBM(U.C.CDM., S0=S0, mu=rep(r, d), sigma=sigma, T=T)
S.C.MO.  <- rGeoBM(U.C.MO.,  S0=S0, mu=rep(r, d), sigma=sigma, T=T)


### 2.3 Functional Calculation #################################################

erT <- exp(-r*T)

## Using pseudo-random samples

## Basket call
basket.t.CDM <- erT * mean(payoff(K, N=N, S0=S0, S=S.t.CDM))
basket.C.CDM <- erT * mean(payoff(K, N=N, S0=S0, S=S.C.CDM))
basket.C.MO  <- erT * mean(payoff(K, N=N, S0=S0, S=S.C.MO))

##Best-of
best.of.t.CDM <- erT * mean(payoff(K, N=N, S0=S0, S=S.t.CDM,method="best.of"))
best.of.C.CDM <- erT * mean(payoff(K, N=N, S0=S0, S=S.C.CDM,method="best.of"))
best.of.C.MO  <- erT * mean(payoff(K, N=N, S0=S0, S=S.C.MO,method="best.of"))

## Worst of call
worst.of.t.CDM <- erT * mean(payoff(K, N=N, S0=S0, S=S.t.CDM, method="worst.of"))
worst.of.C.CDM <- erT * mean(payoff(K, N=N, S0=S0, S=S.C.CDM, method="worst.of"))
worst.of.C.MO  <- erT * mean(payoff(K, N=N, S0=S0, S=S.C.MO,  method="worst.of"))



## Using quasi-random samples

## Basket call
basket.t.CDM. <- erT * mean(payoff(K, N=N, S0=S0, S=S.t.CDM.))
basket.C.CDM. <- erT * mean(payoff(K, N=N, S0=S0, S=S.C.CDM.))
basket.C.MO.  <- erT * mean(payoff(K, N=N, S0=S0, S=S.C.MO.))

##Best-of
best.of.t.CDM. <- erT * mean(payoff(K, N=N, S0=S0, S=S.t.CDM.,method="best.of"))
best.of.C.CDM. <- erT * mean(payoff(K, N=N, S0=S0, S=S.C.CDM.,method="best.of"))
best.of.C.MO.  <- erT * mean(payoff(K, N=N, S0=S0, S=S.C.MO.,method="best.of"))


## Worst of call
worst.of.t.CDM. <- erT * mean(payoff(K, N=N, S0=S0, S=S.t.CDM., method="worst.of"))
worst.of.C.CDM. <- erT * mean(payoff(K, N=N, S0=S0, S=S.C.CDM., method="worst.of"))
worst.of.C.MO.  <- erT * mean(payoff(K, N=N, S0=S0, S=S.C.MO.,  method="worst.of"))

### 2.4 Results ################################################################


res <- array(, dim=c(3,2,2), dimnames=list(type=c("basket", "worst.of",'best.of'),
                                           copula=c("Clayton", paste0("t", nu)),
                                           method=c("CDM", "MO")))
res["basket",,]     <- matrix(c(basket.C.CDM., NA, basket.C.MO., basket.t.CDM.), ncol=2)
res["worst.of",,]   <- matrix(c(worst.of.C.CDM., NA, worst.of.C.MO., worst.of.t.CDM.), ncol=2)
res["best.of",,]     <- matrix(c(best.of.C.CDM., NA, best.of.C.MO., best.of.t.CDM.), ncol=2)
res


#####ANNs####
library(keras)
library(rsample)

sample<-sde::GBM(x=100,N=1000000,r = r,sigma = 0.8,T = 1)


#Generate Dataset
mydata <- NULL
mydata$Stock <- sample
mydata$Strike <- sample*runif(length(sample),min = 0.4,max=1)
mydata$Time <- runif(n=length(sample))
mydata$sigma <- runif(n=length(sample), min = 0.1, max = 0.8)
mydata$r <-runif(n=length(sample),min = 0.01, max=0.05)
mydata$BS <-  EU_call_bs(S = mydata$Stock, t = 0, T = mydata$Time, r = mydata$r, K = mydata$Strike, sigma = as.numeric(mydata$sigma));

mydata <- as.data.frame(mydata)

#Split dataset
split <- rsample::initial_split(mydata, prop = .7, strata='BS' )

train <- rsample::training(split)
test  <- rsample::testing(split)

# Create & standardize feature sets
# training features
train_x <- train %>% dplyr::select(-BS)
mean    <- colMeans(train_x)
std     <- apply(train_x, 2, sd)
train_x <- scale(train_x, center = mean, scale = std)
# testing features
test_x <- test %>% dplyr::select(-BS)
test_x <- scale(test_x, center = mean, scale = std)

# Create & transform response sets
train_y <- log(train$BS)
test_y  <- log(test$BS)


#MODEL 1
model_1 <- keras_model_sequential() %>%
  layer_dense(units =  200,activation = "relu", input_shape = ncol(train_x)) %>%
  layer_dense(units = 130,activation = "relu") %>%
  layer_dense(units = 50,activation = "relu") %>%
  layer_dense(units = 1)  %>%
  
  # backpropagation
  compile(
    optimizer = "rmsprop",
    loss = "mse",
    metrics = c("mae")
  )

learn_1 <- model_1 %>% fit(
  x = train_x,
  y = train_y,
  epochs = 45,
  batch_size = 256,
  validation_split = .2,
  verbose = TRUE,
)


#MODEL 2 

model_2 <- keras_model_sequential() %>%
  layer_dense(units =  200,activation = "relu", input_shape = ncol(train_x),kernel_regularizer = regularizer_l2(0.001)) %>%
  layer_batch_normalization() %>%
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 130,activation = "relu",kernel_regularizer = regularizer_l2(0.001)) %>%
  layer_batch_normalization() %>%
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 50,activation = "relu",kernel_regularizer = regularizer_l2(0.001)) %>%
  layer_batch_normalization() %>%
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 1)  %>%
  
  compile(
    optimizer = "rmsprop", #optimizer_adam(lr = 0.01) 
    loss = "mse",
    metrics = c("mae")
  )
learn_2 <- model_2 %>% fit(
  x = train_x,
  y = train_y,
  epochs = 45,
  batch_size = 256,
  validation_split = .2,
  verbose = TRUE,
  callbacks = list(
    #callback_early_stopping(patience = 10),
    callback_reduce_lr_on_plateau(patience = 5))
)

#Storing Prediction
risultati <- NULL
risultati$predcted_values_model_1<- model_1 %>% predict(test_x)
risultati$true_value <- test_y
risultati$predcted_values_model_2<- model_2 %>% predict(test_x)
risultati$pred_value_model_1_converted <- exp(risultati$predcted_values_model_1)
risultati$true_converted <- exp(risultati$true_value)
risultati$pred_value_model_2_converted <- exp(risultati$predcted_values_model_2)
risultati$S_K <- test[,1] / test[,2]
risultati$Err_model_1 <- abs(risultati$true_converted - risultati$pred_value_model_1_converted )
risultati$Err_model_2 <- abs(risultati$true_converted - risultati$pred_value_model_2_converted )
risultati <- as.data.frame(risultati)


plot(x=risultati$true_converted, y=risultati$pred_value_model_1_converted , cex= 0.001, xlab='Actual', ylab='Predicted Model 1')
plot(x=risultati$true_converted, y=risultati$pred_value_model_2_converted , cex= 0.001, xlab='Actual', ylab='Predicted Model 2')
plot(x=risultati$S_K, risultati$Err_model_1, xlab='S / K', ylab='|BS Price - Predicted model 1|',cex=0.01)
plot(x=risultati$S_K, risultati$Err_model_2, xlab='S / K', ylab='|BS Price - Predicted model 2|',cex=0.01)

learn_1
learn_2

(result_1<- model_1 %>% evaluate(test_x, test_y))
(result_2<- model_2 %>% evaluate(test_x, test_y))