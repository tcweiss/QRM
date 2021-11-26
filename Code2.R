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


################################
###       (i) MODELS         ###
################################

# M1: Portfolio returns with empirically distributed returns.
# no calculations necessary
rets_M1 <- rets

# M2: Portfolio returns with bivariate Gaussian distributed returns.

# mle of the bivariate gaussian distribution

rets_M2_A <- (rets$AAPL-mean(rets$AAPL)/sd(rets$AAPL))
rets_M2_T <- (rets$TSLA-mean(rets$TSLA)/sd(rets$TSLA))
rets_M2 <- cbind(rets_M2_A, rets_M2_T)
mle_M2 <- mvnorm.mle(rets)

mle_M2$mu
mle_M2$sigma

# M3: Portfolio returns with Gaussian distributed returns.
test <- pnorm(rets_M2)
hist(test$TSA)

u1 <- pnorm((rets$AAPL - mle_M2$mu[1])/sqrt(mle_M2$sigma[1,1])) %>% as.vector() 
u2 <- pnorm((rets$TSLA - mle_M2$mu[2])/sqrt(mle_M2$sigma[2,2])) %>% as.vector()
plot(u1)
hist(u2)

mle_M3 <- BiCopEst(u1, u2, family = 4, method = "mle")
summary(mle_M3)

x1 <- qnorm(u1, mean = mle_M2$mu[1], sd = sqrt(mle_M2$sigma[1,1]))
x2 <- qnorm(u2, mean = mle_M2$mu[2], sd = sqrt(mle_M2$sigma[2,2]))
rets_M3 <- cbind(x1,x2)

# M4: Portfolio returns with t-distributed returns.
eps1 <- (rets$AAPL - mle_M2$mu[1])/sqrt(mle_M2$sigma[1,1]) %>% as.vector()
eps2 <- (rets$TSLA - mle_M2$mu[2])/sqrt(mle_M2$sigma[2,2]) %>% as.vector()

v1 <- fitdistr(eps1, "t")$estimate[3]
v2<- fitdistr(eps2, "t")$estimate[3]

u11 <- pt(eps1, df = v1)
u22 <- pt(eps2, df = v2)

mle_M4 <- BiCopEst(u11, u22, family = 1, method = "mle")
summary(mle_M4)


################################
###      (ii) VaR M1         ###
################################

# Create 10k simulations using model M1.

#prets_M1 <- Return.portfolio(rets_M1, weights = c(0.3,0.7), rebalance_on = "week", geometric = TRUE)
prets_M1 <- 0.3*rets$AAPL + 0.7*rets$TSLA
pretsM1 <- sample(prets_M1, 10000, replace = TRUE)

# Estimate 1-day Value at Risk at confidence levels of 90%, 95% and 99%.
VAR_M190 <- VaR(prets_M1, 0.9, method = "historical", invert = FALSE, operational = FALSE)
VAR_M195 <- VaR(prets_M1, 0.95, method = "historical", invert = FALSE, operational = FALSE)
VAR_M199 <- VaR(prets_M1, 0.99, method = "historical", invert = FALSE, operational = FALSE)
VAR_M1 <- c(VAR_M190,VAR_M195,VAR_M199)

# Compute expected shortfalls at confidence levels of 90%, 95% and 99%.
CVAR_M190 <- ES(prets_M1, 0.9, method = "historical", invert = FALSE, operational = FALSE)
CVAR_M195 <- ES(prets_M1, 0.95, method = "historical", invert = FALSE, operational = FALSE)
CVAR_M199 <- ES(prets_M1, 0.99, method = "historical", invert = FALSE, operational = FALSE)
CVAR_M1 <- c(CVAR_M190,CVAR_M195,CVAR_M199)

# Display result
result_M1 <- tibble("VaR" = VAR_M1, "ES" = CVAR_M1)
row.names(result_M1) <- c("90%", "95%", "99%")
grid.table(result_M1)


################################
###    (iii) VaR M2          ###
################################

# Create 10k simulations using models M2
rets_M2 <- as.data.frame(mvrnorm(n=10000, mu = mle_M2$mu, Sigma = mle_M2$sigma, empirical = TRUE))

# plot results
pairs.panels(rets_M2)

# simulated portfolio returs
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

n <- 10000
rho <- mle_M4$par

sim.GC <- function(n, rho, qmarg1, qmarg2){
  R <- rbind(c(1,rho),c(rho,1))
  dat <- rmvnorm(n, mean = c(0,0), sigma = R)
  dat[,1] <- qmarg1(pnorm(dat[,1]))
  dat[,2] <- qmarg2(pnorm(dat[,2]))
  return(dat)
}

set.seed(1234)

# gaussian marginals
q1 <- function(p) qnorm(p, mean = mle_M2$mu[1], sd = sqrt(mle_M2$sigma[1,1]) )
q2 <- function(p) qnorm(p, mean = mle_M2$mu[2], sd = sqrt(mle_M2$sigma[2,2]))

rets_M4 <- sim.GC(n,rho = rho,q1,q2)
rets_M4 <- data.frame(x = sim[,1], y = sim[,2])

# plot results
p <- ggplot(rets_M4, aes(x, y)) + geom_point() + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  xlab("X") + ylab("Y") + ggtitle(expression(paste(rho, "=0.5")))
ggExtra::ggMarginal(p, type = "histogram")

pairs.panels(rets_M4)

prets_M4 <- rets_M4$x*0.3 + rets_M4$y*0.7

# Estimate 1-week Value at Risk using models M4
VAR_M490 <- VaR(prets_M4, 0.9, method = "gaussian", invert = FALSE, operational = FALSE)
VAR_M495 <- VaR(prets_M4, 0.95, method = "gaussian", invert = FALSE, operational = FALSE)
VAR_M499 <- VaR(prets_M4, 0.99, method = "gaussian", invert = FALSE, operational = FALSE)
VAR_M4 <- c(VAR_M490,VAR_M495,VAR_M499)

# Compute expected shortfalls at confidence levels of 90%, 95% and 99% using models M4
CVAR_M490 <- CVaR(prets_M4, 0.9, method = "gaussian", invert = FALSE, operational = FALSE)
CVAR_M495 <- CVaR(prets_M4, 0.95, method = "gaussian", invert = FALSE, operational = FALSE)
CVAR_M499 <- CVaR(prets_M4, 0.99, method = "gaussian", invert = FALSE, operational = FALSE)
CVAR_M4 <- c(CVAR_M490,CVAR_M495,CVAR_M499)

# Display result M4
result_M4 <- tibble("VaR_M4" = VAR_M4, "ES_M4" = CVAR_M4)
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
    n <- 10000
    rho <- mle_M4$par
    
    sim.GC <- function(n, rho, qmarg1, qmarg2){
      R <- rbind(c(1,rho),c(rho,1))
      dat <- rmvnorm(n, mean = c(0,0), sigma = R)
      dat[,1] <- qmarg1(pnorm(dat[,1]))
      dat[,2] <- qmarg2(pnorm(dat[,2]))
      return(dat)
    }
    
    set.seed(1234)
    
    # gaussian marginals
    q1 <- function(p) qnorm(p, mean = mle_M2$mu[1], sd = sqrt(mle_M2$sigma[1,1]) )
    q2 <- function(p) qnorm(p, mean = mle_M2$mu[2], sd = sqrt(mle_M2$sigma[2,2]))
    
    rets <- sim.GC(n,rho = rho,q1,q2)
    rets  <- data.frame(x = sim[,1], y = sim[,2])
    rets <- rets$x*0.3 + rets$y*0.7
  }else {
    retsA <- sample(n$AAPL, 10000, replace = TRUE)
    retsT <- sample(n$TSLA, 10000, replace = TRUE)
    rets <- 0.3*retsA + 0.7*retsT
    VAR <- VaR(rets, c, method = "historical", invert = FALSE, operational = FALSE)
    ES <- CVaR(rets, c, method = "historical", invert = FALSE, operational = FALSE)
  }
  
  VAR <- VaR(rets, c, method = "gaussian", invert = FALSE, operational = FALSE)
  ES <- CVaR(rets, c, method = "gaussian", invert = FALSE, operational = FALSE)
  
  result <- data.frame(VaR = VAR, ES = ES)
  colnames(result) <- c(paste0("VaR M",model), paste0("ES M",model))
  rownames(result) <- paste0(c, "%")
  return(result)
}

RM.N(n=rets_N400, c=0.9, model = 4)


################################
###  (v) VaR ROLLING WINDOW  ###
################################

# Estimate 1-week Value at Risk over 100-day rolling window using models M2 and
# M4.



# Compute no. of violations using next out-of-sample portfolio return.