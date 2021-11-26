
library(tidyverse)
library(magrittr)
library(readxl)
library(xts)
library(PerformanceAnalytics)
library(stats4)
library(VineCopula)
library(Rfast)
library(MASS)


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

M2 <- mvnorm.mle(rets)


# M3: Portfolio returns with Gaussian distributed returns.

u1 <- pnorm((rets$AAPL-M2$mu[1])/sqrt(M2$sigma[1,1]))
u2 <- pnorm((rets$TSLA-M2$mu[2])/sqrt(M2$sigma[2,2]))

M3 <- BiCopEst(u1, u2, family = 4, method = "mle")


# M4: Portfolio returns with t-distributed returns.

eps1 <- (rets$AAPL-M2$mu[1])/sqrt(M2$sigma[1,1])
eps2 <- (rets$TSLA-M2$mu[2])/sqrt(M2$sigma[2,2])

v1 <- fitdistr(eps1, densfun = "t")$estimate[3]
v2 <- fitdistr(eps2, densfun = "t")$estimate[3]

u1 <- pt(eps1, v1)
u2 <- pt(eps2, v2)

M4 <- BiCopEst(u1, u2, family = 1, method = "mle")


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
                         geometric = TRUE)


rets_1 <- sample(pf_1, 10000, replace = TRUE)



# Estimate 1-week VaR and ES at confidence levels of 90%, 95% and 99%.

results_1 <- tibble("alpha" = rep(NA_real_, 3),
                    "VaR" = rep(NA_real_, 3),
                    "ES" = rep(NA_real_, 3))

conf <- c(0.9, 0.95, 0.99)

for (i in 1:3) {
  
  results_1$alpha[i] <- conf[i]
  
  results_1$VaR[i] <- VaR(rets_1, 
                          p = conf[i], 
                          method = "historical", 
                          invert = FALSE)
  
  results_1$ES[i] <- ES(rets_1, 
                        p = conf[i],
                        method = "historical", 
                        invert = FALSE, 
                        operational = FALSE)
  
}

rm("i")
rm("conf")


################################
###    (iii) VaR M2 & M4     ###
################################


# Create 10k simulations using models M2 and M4.



# Estimate 1-week Value at Risk.



# Compute expected shortfalls at confidence levels of 90%, 95% and 99%.




################################
###        (iv) VaR N        ###
################################


# Create subsets with last 100, 200, 300 and 400 observations.



# Estimate Value at Risk (again use 1-week VaR?).




################################
###  (v) VaR ROLLING WINDOW  ###
################################

# Estimate 1-week Value at Risk over 100-day rolling window using models M2 and
# M4.



# Compute no. of violations using next out-of-sample portfolio return.






