
library(tidyverse)
library(magrittr)
library(readxl)
library(xts)
library(PerformanceAnalytics)
library(stats4)
library(copula)


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



# M2: Portfolio returns with bivariate Gaussian distributed returns.



# M3: Portfolio returns with Gaussian distributed returns.



# M4: Portfolio returns with t-distributed returns.




################################
###      (ii) VaR M1         ###
################################


# Create 10k simulations using model M1.



# Estimate 1-day Value at Risk.



# Compute expected shortfalls at confidence levels of 90%, 95% and 99%.




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






