library(tidyverse)
library(readxl)
library(xts)
library(PerformanceAnalytics)
library(stats4)
library(MASS)
library(gridExtra)
library(Rfast)
library(spatstat)
library(copula)
library(VineCopula)
library(mvtnorm)
library(sn)

################################
###          SETUP           ###
################################

# Set working directory
setwd("~/Desktop/QRM")

# Import data.
data <- read_excel("data.xlsx")


# Extract coredata and index.
first <- as.matrix(data[,-1])
second <- as.Date(data$Dates)


# Create xts object with raw data.
data <- xts(first, second)
rm(first)
rm(second)


# Create xts object with net returns.
rets <- na.omit(Return.calculate(data, method = "discrete"))
rets <- rets[,c(1,6)] # only Apple and Tesla returns

# Descriptive statistics
mu_A <- mean(rets$AAPL)
mu_T <- mean(rets$TSLA)
mu  <- rbind(mu_A, mu_T)
coVarMat <- cov(rets)
mu_P <- mean(rets_portfolio)
var_P <- var(rets_portfolio)

# Empirical distribution of the returns
histo_A <- hist(rets$AAPL)
hist_T <- hist(rets$TSLA)
hist_P <- hist(rets_portfolio)



################################
###       (i) MODELS         ###
################################

# M1: Portfolio returns with empirically distributed returns.
# no calculations necessary

# M2: Portfolio returns with bivariate Gaussian distributed returns.
test <- dnorm(rets)
# mle of the bivariate gaussian distribution
mle_M2 <- mvnorm.mle(rets)
mle_M2$mu
mle_M2$sigma

# M3: Portfolio returns with Gaussian distributed returns.

u1 <- dnorm((rets$AAPL - mle_M2$mu[1])/mle_M2$sigma[1,1]) %>% as.vector()
u2 <- dnorm((rets$TSLA - mle_M2$mu[2])/mle_M2$sigma[2,2]) %>% as.vector()
mle_M3 <- BiCopEst(u1, u2, family = 4, method = "mle")
summary(mle_M3)

# M4: Portfolio returns with t-distributed returns.
mle_M4 <- BiCopEst(u1, u2, family = 1, method = "mle")
summary(mle_M4)


################################
###      (ii) VaR M1         ###
################################

# Create 10k simulations using model M1.
retsA_10k <- sample(rets$AAPL, 10000, replace = TRUE)
retsT_10k <- sample(rets$TSLA, 10000, replace = TRUE)
rets_M1 <- 0.3*retsA_10k + 0.7*retsT_10k

# Estimate 1-day Value at Risk at confidence levels of 90%, 95% and 99%.
VAR_M190 <- VaR(rets_M1, 0.9, method = "historical")
VAR_M195 <- VaR(rets_M1, 0.95, method = "historical")
VAR_M199 <- VaR(rets_M1, 0.99, method = "historical")
VAR_M1 <- round(c(VAR_M190,VAR_M195,VAR_M199),4)

# Compute expected shortfalls at confidence levels of 90%, 95% and 99%.
CVAR_M190 <- CVaR(rets_M1, 0.9, method = "historical")
CVAR_M195 <- CVaR(rets_M1, 0.95, method = "historical")
CVAR_M199 <- CVaR(rets_M1, 0.99, method = "historical")
CVAR_M1 <- round(c(CVAR_M190,CVAR_M195,CVAR_M199),4)

# Display result
result_M1 <- data.frame(VaR_M1 = VAR_M1, ES_M1 = CVAR_M1)
row.names(result_M1) <- c("90%", "95%", "99%")
grid.table(result_M1)


################################
###    (iii) VaR M2          ###
################################

# Create 10k simulations using models M2
rets_M2 <- as.data.frame(mvrnorm(n=10000, mu = mle_M2$mu, Sigma = mle_M2$sigma, empirical = TRUE))
pairs.panels(rets_M2)
rets_M2 <- 0.3*rets_M2$AAPL + 0.7*rets_M2$TSLA

# Estimate 1-week Value at Risk using models M2
VAR_M290 <- VaR(rets_M2, 0.9, method = "gaussian")
VAR_M295 <- VaR(rets_M2, 0.95, method = "gaussian")
VAR_M299 <- VaR(rets_M2, 0.99, method = "gaussian")
VAR_M2 <- round(c(VAR_M290,VAR_M295,VAR_M299),4)

# Compute expected shortfalls at confidence levels of 90%, 95% and 99% using models M2
CVAR_M290 <- CVaR(rets_M2, 0.9, method = "gaussian")
CVAR_M295 <- CVaR(rets_M2, 0.95, method = "gaussian")
CVAR_M299 <- CVaR(rets_M2, 0.99, method = "gaussian")
CVAR_M2 <- round(c(CVAR_M290,CVAR_M295,CVAR_M299),4)

# Display result M2
result_M2 <- data.frame(VaR_M2 = VAR_M2, ES_M2 = CVAR_M2)
row.names(result_M2) <- c("90%", "95%", "99%")
grid.table(result_M2)


################################
###    (iii) VaR M4          ###
################################

# Create 10k simulations using models M4.

set.seed(1234)
set.seed(100)
m <- 2
n <- 10000
sigma <- matrix(c(1, rho,
                  rho, 1), 
                nrow=2)
z <- mvrnorm(n,mu=rep(0, m),Sigma=sigma,empirical=T)
pairs.panels(z)
u <- pnorm(z)
pairs.panels(u)

x1 <- qnorm(u[,1],mean = mle_M2$mu[1], sd = mle_M2$sigma[1,1])
x2 <- qnorm(u[,2],mean = mle_M2$mu[2], sd = mle_M2$sigma[2,2])


normal.cop <- normalCopula(param=rho, dim = 2)
m <- pobs(as.matrix(cbind(rets$AAPL,rets)))
fit <- fitCopula(normal.cop,m,method='ml')
coef(fit)
set.seed(500)
m <- pobs(as.matrix(cbind(cree,yahoo)))
fit <- fitCopula(t.cop,m,method='ml')
coef(fit)



ns <- 10000
rho <- mle_M4$par
cov <- rho*sqrt(mle_M2$sigma[1,1]*mle_M2$sigma[2,2])

sim.GC <- function(n, rho, qmarg1, qmarg2){
  R <- rbind(c(1,rho),c(rho,1))
  dat <- rmvnorm(n, mean = c(0,0), sigma = R)
  dat[,1] <- qmarg1(pnorm(dat[,1]))
  dat[,2] <- qmarg2(pnorm(dat[,2]))
  return(dat)
}

q1 <- function(p) qnorm(p, mean = mle_M2$mu[1], sd = mle_M2$sigma[1,1])
q2 <- function(p) qnorm(p, mean = mle_M2$mu[2], sd = mle_M2$sigma[2,2])

rets_M4 <- as.data.frame(sim.GC(n = ns,rho = rho,q1,q2))
rets_M4 <- 0.3*rets_M4$V1 + 0.7*rets_M4$V2

df <- data.frame(x = rets_M4[,1], y = rets_M4[,2])
p <- ggplot(df, aes(x, y)) + geom_point() + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  xlab("X") + ylab("Y") + ggtitle(expression(paste(rho, "=0.5")))
ggExtra::ggMarginal(p, type = "histogram")

# Estimate 1-week Value at Risk using models M4
VAR_M490 <- VaR(rets_M4, 0.9, method = "gaussian")
VAR_M495 <- VaR(rets_M4, 0.95, method = "gaussian")
VAR_M499 <- VaR(rets_M4, 0.99, method = "gaussian")
VAR_M4 <- c(VAR_M490,VAR_M495,VAR_M499)

# Compute expected shortfalls at confidence levels of 90%, 95% and 99% using models M4
CVAR_M490 <- CVaR(rets_M4, 0.9, method = "gaussian")
CVAR_M495 <- CVaR(rets_M4, 0.95, method = "gaussian")
CVAR_M499 <- CVaR(rets_M4, 0.99, method = "gaussian")
CVAR_M4 <- c(CVAR_M490,CVAR_M495,CVAR_M499)

# Display result M4
result_M4 <- data.frame(VaR_M4 = VAR_M4, ES_M4 = CVAR_M4)
row.names(result_M4) <- c("90%", "95%", "99%")
grid.table(result_M4)


################################
###    (iv) VaR & ES N       ###
################################


# Create subsets with last 100, 200, 300 and 400 observations.
rets_N100 <- rets[1:99,]
rets_N200 <- rets[1-199,]
rets_N300 <- rets[1-299,]
rets_N400 <- rets[1-399,]

# Estimate Value at Risk and ES (again use 1-week VaR?)
RM.N <- function(n, c, model){
  if (model == 2){
    mle <- mvnorm.mle(n)
    rets <- as.data.frame(mvrnorm(n=10000, mu = mle$mu, Sigma = mle$sigma, empirical = TRUE))
    rets <- 0.3*rets$AAPL + 0.7*rets$TSLA
  }else if (model == 4) {
    mle_M2 <- mvnorm.mle(n)
    u1 <- dnorm((rets$AAPL - mle_M2$mu[1])/mle_M2$sigma[1,1]) %>% as.vector()
    u2 <- dnorm((rets$TSLA - mle_M2$mu[2])/mle_M2$sigma[2,2]) %>% as.vector()
    mle <- BiCopEst(u1, u2, family = 1, method = "mle")
    
  }else {
    retsA <- sample(n$AAPL, 10000, replace = TRUE)
    retsT <- sample(n$TSLA, 10000, replace = TRUE)
    rets <- 0.3*retsA + 0.7*retsT
    VAR <- VaR(rets, c, method = "historical")
    ES <- CVaR(rets, c, method = "historical")
  }
  
  VAR <- VaR(rets, c, method = "gaussian")
  ES <- CVaR(rets, c, method = "gaussian")
  
  result <- data.frame(VaR = VAR, ES = ES)
  colnames(result) <- c(paste0("VaR M",model), paste0("ES M",model))
  rownames(result) <- paste0(c, "%")
  return(result)
}

RM.N(n=rets_N400, c=0.9, model = 2)


################################
###  (v) VaR ROLLING WINDOW  ###
################################

# Estimate 1-week Value at Risk over 100-day rolling window using models M2 and
# M4.



# Compute no. of violations using next out-of-sample portfolio return.