library(tidyverse)
library(magrittr)
library(readxl)
library(xts)
library(PerformanceAnalytics)
library(stats4)
library(VineCopula)
library(Rfast)
library(MASS)
library(psych)
library(gridExtra)
library(mvtnorm)
library(sn)
################################
###          SETUP           ###
################################

# Import data.

data <- read_excel("data.xlsx")


# Extract coredata and index.

first <- as.matrix(data[,-1])
second <- as.Date(data$Dates)


# Create xts object from raw data.

data <- xts(first, second)
rm(first)
rm(second)


# Create xts object with net returns and extract relevant stocks.

rets <- Return.calculate(data, method = "discrete")
rets <- rets[-1, c("AAPL", "TSLA")]



################################
###       (i) MODELS         ###
################################

# M1: Portfolio returns with empirically distributed returns.

M1 <- rets



# M2: Portfolio returns with bivariate Gaussian distributed returns.

M2 <- mvnorm.mle(rets) # parameter estimation with mle
M2$mu # mean vector
M2$sigma # covariance matrix



# M3: Portfolio returns with Gaussian distributed returns.

u1 <- pnorm((rets$AAPL-M2$mu[1])/sqrt(M2$sigma[1,1])) # return transformation with gaussian distributed marginals
u2 <- pnorm((rets$TSLA-M2$mu[2])/sqrt(M2$sigma[2,2])) # return transformation with gaussian distributed marginals

M3 <- BiCopEst(u1, u2, family = 4, method = "mle") # parameter estimation with mle
M3$par # Gumble-copula paramter theta



# M4: Portfolio returns with t-distributed returns.

eps1 <- (rets$AAPL-M2$mu[1])/sqrt(M2$sigma[1,1]) # epsilon estimation
eps2 <- (rets$TSLA-M2$mu[2])/sqrt(M2$sigma[2,2]) # epsilon estimation

v1 <- fitdistr(eps1, densfun = "t")$estimate[3] # df estimation with t-distributed epsilon
v2 <- fitdistr(eps2, densfun = "t")$estimate[3] # df estimation with t-distributed epsilon

u1 <- pt(eps1, v1) # return transformation with t-distributed epsilon
u2 <- pt(eps2, v2) # return transformation with t-distributed epsilon

M4 <- BiCopEst(u1, u2, family = 1, method = "mle") # parameter estimation with mle
M4$par # Gaussian-copula parameter rho


rm("eps1")
rm("eps2")
rm("u1")
rm("u2")
rm("v1")
rm("v2")
rm("m")
rm("n")



################################
###      (ii) VaR M1         ###
################################


# Create 10k simulations using model M1.

pf_1 <- Return.portfolio(M1, 
                         weights = c(0.3, 0.7), 
                         rebalance_on = "week", 
                         geometric = TRUE) # portfoöio returns


rets_1 <- sample(pf_1, 10000, replace = TRUE) # simulated portfolio returns



# Estimate 1-week VaR and ES at confidence levels of 90%, 95% and 99%.

results_1 <- tibble("alpha" = rep(NA_real_, 3),
                    "VaR" = rep(NA_real_, 3),
                    "ES" = rep(NA_real_, 3)) # initialize table for results

conf <- c(0.9, 0.95, 0.99) # define confidence intervals

for (i in 1:3) {
  
  results_1$alpha[i] <- conf[i]
  
  results_1$VaR[i] <- VaR(rets_1, 
                          p = conf[i], 
                          method = "historical", 
                          invert = FALSE) # value-at-risk estimation
  
  results_1$ES[i] <- ES(rets_1, 
                        p = conf[i],
                        method = "historical", 
                        invert = FALSE, 
                        operational = FALSE) # expected-shortfall
  
}

rm("i")
rm("conf")


################################
###    (iii) VaR M2 & M4     ###
################################


# Create 10k simulations using models M2 and M4. I add some fantasy dates in the
# second step since PerformanceAnalytics only accepts time series (will be
# deleted afterwards).

# Simulated returns M2
rets_2 <- mvrnorm(10000, mu = M2$mu, Sigma = M2$sigma, empirical = TRUE) %>% 
  as_tibble(.) # return simulation with bivariate distributed retruns

date <- seq(as.Date("2000/1/1"), by = "weeks", length.out = 10000)

rets_2 <- xts(rets_2, date)

rets_2 <- Return.portfolio(rets_2, 
                           weights = c(0.3, 0.7), 
                           rebalance_on = "week", 
                           geometric = TRUE)


# Simulated returns M4
n <- 10000 # number of simulations
rho <- M4$par # correlation parameter Gaussian-Copula

# Build the bivariate distribution
dist_M4 <- mvdc(normalCopula(param = rho, dim = 2), margins = c("norm","norm"),
                paramMargins = list(list(mean = M2$mu[1], sd = sqrt(M2$sigma[1,1])), list(mean = M2$mu[2], sd = sqrt(M2$sigma[2,2]))))

# Sample 1000 observations from the distribution
rets_M4 <- rMvdc(10000, dist_M4) %>% 
  as_tibble(.)
colnames(rets_M4) <- c("AAPL", "TSLA")

# plot results
p <- ggplot(rets_M4, aes(AAPL, TSLA)) + 
  geom_point() + 
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14)) +
  xlab("Apple") + ylab("Tesla") + 
  ggtitle(expression(paste(rho, "=0.5")))

ggExtra::ggMarginal(p, type = "histogram")

pairs.panels(rets_M4)

rets_M4 <- rets_M4$AAPL*0.3 + rets_M4$TSLA*0.7



# Estimate 1-week VaR and ES at confidence levels of 90%, 95% and 99%.


results_2 <- tibble("alpha" = rep(NA_real_, 3),
                    "VaR" = rep(NA_real_, 3),
                    "ES" = rep(NA_real_, 3))

conf <- c(0.9, 0.95, 0.99)

for (i in 1:3) {
  
  results_2$alpha[i] <- conf[i]
  
  results_2$VaR[i] <- VaR(rets_2, 
                          p = conf[i], 
                          method = "historical", 
                          invert = FALSE)
  
  results_2$ES[i] <- ES(rets_2, 
                        p = conf[i],
                        method = "historical", 
                        invert = FALSE, 
                        operational = FALSE)
  
}


results_3 <- tibble("alpha" = rep(NA_real_, 3),
                    "VaR" = rep(NA_real_, 3),
                    "ES" = rep(NA_real_, 3))

rets_M4 <- xts(rets_M4, date)


for (i in 1:3) {
  
  results_3$alpha[i] <- conf[i]
  
  results_3$VaR[i] <- VaR(rets_M4, 
                          p = conf[i], 
                          method = "historical", 
                          invert = FALSE)
  
  results_3$ES[i] <- ES(rets_M4, 
                        p = conf[i],
                        method = "historical", 
                        invert = FALSE, 
                        operational = FALSE)
  
}


rrm("i")
rm("conf")



################################
###        (iv) VaR N        ###
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
    
    date <- seq(as.Date("2000/1/1"), by = "weeks", length.out = 10000)
    
    rets <- xts(rets, date)
    
    rets <- Return.portfolio(rets, 
                               weights = c(0.3, 0.7), 
                               rebalance_on = "week", 
                               geometric = TRUE)
    
  }else if (model == 4) {
   mle <- mvnorm.mle(n)
   
   eps1 <- (n$AAPL-mle$mu[1])/sqrt(mle$sigma[1,1]) # epsilon estimation
   eps2 <- (n$TSLA-mle$mu[2])/sqrt(mle$sigma[2,2]) # epsilon estimation
   
   v1 <- fitdistr(eps1, densfun = "t")$estimate[3] # df estimation with t-distributed epsilon
   v2 <- fitdistr(eps2, densfun = "t")$estimate[3] # df estimation with t-distributed epsilon
   
   u1 <- pt(eps1, v1) # return transformation with t-distributed epsilon
   u2 <- pt(eps2, v2) # return transformation with t-distributed epsilon
   
   M4 <- BiCopEst(u1, u2, family = 1, method = "mle") # parameter estimation with mle
   dist_M4 <- mvdc(normalCopula(param = M4$par, dim = 2), margins = c("norm","norm"),
                    paramMargins = list(list(mean = mle$mu[1], sd = sqrt(mle$sigma[1,1])), list(mean = mle$mu[2], sd = sqrt(mle$sigma[2,2]))))
    
    # Sample 1000 observations from the distribution
    rets <- rMvdc(10000, dist_M4) %>% 
      as_tibble(.)
    colnames(rets) <- c("AAPL", "TSLA")
    rets <- xts(rets, date)
    rets <- Return.portfolio(rets, 
                             weights = c(0.3, 0.7), 
                             rebalance_on = "week", 
                             geometric = TRUE)
    
  }else {
   rets <- sample(rets, 10000, replace = TRUE)
    rets <- xts(rets, date)
    rets <- Return.portfolio(rets, 
                             weights = c(0.3, 0.7), 
                             rebalance_on = "week", 
                             geometric = TRUE)
  }
  
  VAR <- VaR(rets, c, method = "gaussian", invert = FALSE, operational = FALSE)
  ES <- CVaR(rets, c, method = "gaussian", invert = FALSE, operational = FALSE)
  
  result <- data.frame(VaR = VAR, ES = ES)
  colnames(result) <- c(paste0("VaR M",model), paste0("ES M",model))
  rownames(result) <- paste0(c, "%")
  return(result)
}

RM.N(n=rets_N100, c=0.95, model = 2)


################################
###  (v) VaR ROLLING WINDOW  ###
################################

# Estimate 1-week Value at Risk over 100-day rolling window using models M2 and
# M4.

date <- seq(as.Date("2013-01-01"), by = "week", length.out = 10000)

rets_2 <- xts(coredata(rets_2), date)


rets_4 <- xts(coredata(rets_M4), date)


rets_4 <- Return.portfolio(rets_4, 
                           weights = c(0.3, 0.7), 
                           rebalance_on = "week", 
                           geometric = TRUE)


roll_2 <- rollapply(rets_2, width = 100, FUN = VaR,  p = 0.95, method = "historical", invert = FALSE) %>%
  na.trim(.)

roll_4 <- rollapply(rets_4, width = 100, FUN = VaR,  p = 0.95, method = "historical", invert = FALSE) %>%
  na.trim(.)


# Compute no. of violations using next out-of-sample portfolio return.


rets_2 <- rets_2["2014-11-25/"] %>% 
  lag.xts(., -1)

roll_2 <- cbind(roll_2, rets_2)

colnames(roll_2) <- c("VaR", "pf")

rets_4 <- rets_4["2014-11-25/"] %>% 
  lag.xts(., -1)

roll_4 <- cbind(roll_4, rets_4)

names(roll_4) <- c("VaR", "pf")


viol_2 <- roll_2[-9901,] %>%
  coredata() %>% 
  as_tibble() %>% 
  mutate(Violation = (VaR + pf)<0) %>% 
  select(Violation) %>% 
  summarise(., "Perc. Violations" = sum(Violation)/n())

viol_4 <- roll_4[-9901,] %>%
  coredata() %>% 
  as_tibble() %>% 
  mutate(Violation = (VaR + pf)<0) %>% 
  select(Violation) %>% 
  summarise(., "Perc. Violations" = sum(Violation)/n())