library(easypackages)
libraries('yuima')

#Brownian motion
BM <- setModel(drift = "0", diffusion = "1", state.variable = "x", time.variable = "t", solve.variable = "x", xinit = 100)
x <- simulate(BM) #sampling 100 times by default 
plot(x)

grid <- setSampling(0, 1, 1000) #sampling a 1000 times between 0, 1
x1 <- simulate(BM, sampling = grid)
plot(x1)

#Simple SDE
SDE <- setModel(drift = "mu", diffusion = "sigma", state.variable = "x", time.variable = "t", solve.variable = "x", xinit = 100)
x <- simulate(SDE, true.parameter = list(mu = 0.1, sigma = 0.2))
plot(x)

#Geometric brownian motion
GBM <- setModel(drift = "mu*x", diffusion = "sigma*x", state.variable = "x", time.variable = "t", solve.variable = "x", xinit = 100)
x <- simulate(GBM, true.parameter = list(mu = 0.1, sigma = 0.2))
plot(x)

#Probability distribution of GBM
simn <- 1000
sim <- matrix(NA, 101, 1000)
for (i in 1:1000) {
  cat("iteration", i, "\n")
  sim[,i] <- simulate(GBM, true.parameter = list(mu = 0.1, sigma = 0.2))@data@original.data
}

simmean <- colMeans(sim)
hist(simmean, probability = T); lines(density(simmean))
plot.ts(sim[,1:100], plot.type = "single")

#Ornstein–Uhlenbeck process
OU <- setModel(drift = "theta*(mu - x)", diffusion = "sigma", solve.variable = "x", xinit = 0.5)
x <- simulate(OU, true.parameter = list(theta = 2, mu = 0.1, sigma = 0.2))
plot(x)

#CIR Model
CIR <- setModel(drift = "theta*(mu - x)", diffusion = "sigma*(x^0.5)", solve.variable = "x", xinit = 0.15)
x <- simulate(CIR, true.parameter = list(theta = 2, mu = 0.1, sigma = 0.2))
plot(x)

#CKLS Models
CKLS <- setModel(drift = "alpha1 + alpha2*x", diffusion = "alpha3*(x^alpha4)", solve.variable = "x", xinit = 0.5)
grid <- setSampling(0, 1, 1000)
x <- simulate(CKLS, true.parameter = list(alpha1 = 0.1, alpha2 = 0.2, alpha3 = 0.3, alpha4 = 0.4), sampling = grid)
plot(x)

#Hyperbolic Process
HP = setModel(drift = "theta*x/((1 + x^2)^{0.5})", diffusion = "1", solve.var = "x", xinit = 0.5)
x = simulate(HP, true.parameter = list(theta = 1))
plot(x)

#Fitting SDE to data
sim_data <- setModel(drift = "theta*(mu - x)", diffusion = "sigma", solve.var = "x", xinit = 0.5)
x <- simulate(sim_data, true.parameter = list(mu = 0.1, sigma = 0.2, theta = 2))
ini = list(mu = 0.05, sigma = 0.5, theta = 1)
low = list(mu = 0, sigma = 0, theta = 0)
up = list(mu = 0.2, sigma = 2, theta = 3)
mle = qmle(x, start = ini, lower = low, upper = up, method = "L-BFGS-B")
summary(mle)

#Jumps in SDEs
grid <- setSampling(0, 1, 1000)
jump <- list(intensity = "7", df = list("dnorm(z, 0, 0.2)"))
OU <- setModel(drift = "theta*(mu - x)", diffusion = "sigma", solve.variable = "x", xinit = 0.2, jump.coeff = "1", measure = jump, measure.type = "CP")
x <- simulate(OU, true.parameter = list(theta = 2, mu = 0.1, sigma = 0.2), sampling = grid)
plot(x)

# Fractional Brownian Motion
grid <- setSampling(0, 1, 1000)
OU0 <- setModel(drift = "theta*(mu - x)", diffusion = "sigma", solve.variable = "x", xinit = 0.2, hurst = 0.3)
OU1 <- setModel(drift = "theta*(mu - x)", diffusion = "sigma", solve.variable = "x", xinit = 0.2, hurst = 0.7)
x0 <- simulate(OU0, true.parameter = list(theta = 2, mu = 0.1, sigma = 0.2), sampling = grid)
x1 <- simulate(OU1, true.parameter = list(theta = 2, mu = 0.1, sigma = 0.2), sampling = grid)
plot(cbind(x0, x1), plot.type = "m", col = c(2, 3))

#Correlated Brownian Motion
sol_var <- c("x1", "x2", "x3")
dr <- c("b1*x1", "b2*x2", "b3*x3")
cv <- c(2, 1, 3, 1, 4, 2, 3, 2, 5)
cv <- matrix(cv, 3, 3)
cv == t(cv)
cv <- chol(cv)
CBM <- setModel(drift = dr, diffusion = cv, solve.variable = sol_var, xinit = c(1, 2, 3))
x <- simulate(CBM, true.parameter = list(b1 = 0.5, b2 = 0.6, b3 = 0.7))
plot(x)

#Multidimensional Brownian Motion
sol_var <- c("r", "u")
dr <- c("theta - alpha*r + u", "-b*u")
diff <- matrix(c("sigma1", "0", "0", "sigma2"), 2, 2)
MBM <- setModel(drift = dr, diffusion = diff, solve.variable = sol_var, xinit = c(0.1, 0.2))
x <- simulate(MBM, true.parameter = list(theta = 1, alpha = 1, b = 1, sigma1 = 2, sigma2 = 2))
plot(x)

# Heston Model
sol_var <- c("s", "sigma")
dr <- c("mu*s", "k*(theta-sigma)")
diff <- c("c1*s*sigma^0.5", "c2*s*sigma^0.5", "0", "c3*eta*sigma^0.5")
diff <- matrix(diff, 2, 2, byrow = T)
cov <- matrix(c(2,0.7,0.7,5), 2, 2)
cov <- chol(cov) 
HM = setModel(drift = dr, diffusion = diff, solve.variable = sol_var, xinit = c(50,5))
x = simulate(HM, true.param = list(theta = 1, eta = 1, mu = 1, k = 2, c1 = cov[1,1], c2 = cov[1,2], c3 = cov[2,2]))
plot(x)








