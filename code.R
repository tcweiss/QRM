
library(tidyverse)
library(magrittr)
library(readxl)
library(xts)
library(PerformanceAnalytics)
library(stats4)
library(VineCopula)
library(copula)
library(Rfast)
library(MASS)
library(psych)
library(gridExtra)


###########################################################################
###########################################################################
###                                                                     ###
###                                SETUP                                ###
###                                                                     ###
###########################################################################
###########################################################################

# Import data from Excel into R.

data <- read_excel("data.xlsx")


# Extract stock prices (as a matrix) and dates (as a vector).

prices <- as.matrix(data[,-1])
dates <- as.Date(data$Dates)


# Create xts (time series) object using price matrix and date vector.

data <- xts(prices, dates)


# Create xts time series object net returns from stock price data. For a price
# P(t) in a given period t, Return.calculate() computes the net return using
# [P(t)-P(t-1)]/P/(t-1). The output is a xts object with net returns for all
# stocks in the data set. The second line then removes the first line from this
# object (= NA due to return calculation), and extracts two specific stocks
# which will be needed from here on (AAPL and TSLA).

rets <- Return.calculate(data, method = "discrete")
rets <- rets[-1, c("AAPL", "TSLA")]


# Function to compute VaR. The if() condition first checks the input data type
# and converts it to a vector. The ecdf() function then computes the empirical
# cumulative distribution of this vector. Finally, the last code chunk computes
# VaR by using the definition of the generalized inverse. It subsets the vector
# to only include values for which the cumulative distribution gives at least
# the desired confidence level, and then returns the smallest one.

VaR <- function(x, p) {
  
  if(is.xts(x) == TRUE) {
    
    x <- x %>% 
      coredata() %>% 
      as.vector() 
    
  } else {
    
    x <- as.vector(x)
    
  }
  
  Fn <- ecdf(x)
  
  x[Fn(x) >= 1-p] %>% 
    min(.)*(-1) %>% 
    return()
  
}


# Function to compute ES. The first three chunks are from the VaR formula above.
# However, instead of returning VaR, the function saves it in an object (VAR),
# which is the used to compute the tail conditional expectation (TCE). Finally,
# the last line applies the formula from chapter 3 to compute the expected
# shortfall (ES).

ES <- function(x, p) {
  
  if(is.xts(x) == TRUE) {
    
    x <- x %>% 
      coredata() %>% 
      as.vector() 
    
  } else {
    
    x <- as.vector(x)
    
  }
  
  Fn <- ecdf(x)
  
  VAR <- x[Fn(x) >= 1-p] %>% 
    min(.)*(-1)
  
  TCE <- -mean(x[x<=-VAR])
  
  return(TCE + ((1/(1-p))*(length(x[x<=-VAR])/length(x)) - 1)*(TCE-VAR))
  
}


# Clean up environment.

rm("data")
rm("prices")
rm("dates")



###########################################################################
###########################################################################
###                                                                     ###
###                             COMPUTATION                             ###
###                                                                     ###
###########################################################################
###########################################################################

################################
###       (i) Models         ###
################################

# M1: Stock returns with empirical distribution. No further
# calculation necessary.

M1 <- rets


# M2: Stock returns with bivariate Gaussian distribution using maximum
# likelihood estimation. If a bivariate set of returns X=(x1, x2) is
# (supposedly) normally distributed, then Y=exp(X) has a multivariate log-normal
# distribution with expectation E[Y_i] = exp(µ_i + 0.5*Σ_ii) and covariance
# matrix Var(Y_ij) = exp(µ_i+µ_j+0.5(Σ_ii+Σ_jj))*(exp(Σ_ij)-1). The mvnorm.mle()
# function applies this assumption to the stock return data and estimates the
# corresponding mean and covariance matrix using MLE-estimation.

M2 <- mvnorm.mle(rets)


# M3: Stock returns with bivariate Gaussian distribution and Gumbel-copula. The
# parameters of the bivariate Gaussian are given by M2 (see previous code
# chunk); we therefore still need to get the Gumbel-copula. The first step is to
# find u1 and u2, which we do by computing the normal distribution (using pnorm)
# of the standardizing net returns (using mu and sd from M2). These results are
# then passed to BiCopEst(), where we set family=4 (for Gumbel-copula) and
# method=mle (for MLE estimation). The BiCopEst() function uses u1 and u2 in the
# formula for the Gumbel-copula and uses MLE to estimate the only unknown
# parameter, θ.

u1 <- pnorm((rets$AAPL-M2$mu[1])/sqrt(M2$sigma[1,1]))
u2 <- pnorm((rets$TSLA-M2$mu[2])/sqrt(M2$sigma[2,2]))

test <- BiCopEst(u1, u2, family = 1, method = "mle")


# M4: Stock returns with t-distribution and Gaussian copula. The first step is
# to solve for the ε of each stock, which in our case means standardizing net
# returns like in M3 (i.e., we use the parameters from M2). One then needs to
# find the degrees of freedom of each ε. This is done in the next two lines. The
# fitdistr() function uses the respective ε as input for the standardized
# t-distribution, given by ft,1(x; ν) = Γ[(ν+1)/2]/[sqrt(vπ)Γ(v/2)]*(1+(x^2)/v)^-[(v+1)/2], 
# and uses MLE to find the only unkown parameter, v. After obtaining estimates
# the two v, we can find u1 and u2. The approach similar to the one in M3, but
# because we have a t-distribution, we use pt() to get the t-distribution and we
# input both ε and the respective v. Finally, we use u1 and u2 as input for
# BiCopEst() again, but we set family=1 to compute the parameter (rho) of a
# Gaussian copula.

eps1 <- (rets$AAPL-M2$mu[1])/sqrt(M2$sigma[1,1])
eps2 <- (rets$TSLA-M2$mu[2])/sqrt(M2$sigma[2,2])

v1 <- fitdistr(eps1, densfun = "t")$estimate[3]
v2 <- fitdistr(eps2, densfun = "t")$estimate[3]

u1 <- pt(eps1, v1)
u2 <- pt(eps2, v2)

M4 <- BiCopEst(u1, u2, family = 1, method = "mle")


# Clean up environment.

rm("eps1")
rm("eps2")
rm("u1")
rm("u2")
rm("v1")
rm("v2")

################################
###      (ii) VaR M1         ###
################################

# Create portfolio using M1. The Return.portfolio() function computes the
# portfolio returns using a time series object of asset returns. These returns
# are given by the first argument, which in our case contains the empirical net
# return distributions of AAPL and TSLA. The second argument specifies how these
# returns should be weighted when computing portfolio returns. Third one
# indicates that the weighing should take place every week, or in our case every
# observation. Without doing this, weights would be applied in the first row
# only, which is not what we want. Finally, we specify that the (weighted)
# geometric mean should be used. While the latter makes no noticable difference,
# it seems like the more correct way to compute portfolio returns.

pf_M1 <- Return.portfolio(M1, 
                         weights = c(0.3, 0.7), 
                         rebalance_on = "week", 
                         geometric = TRUE)


# Create 10k simulations using portfolio 1. Since this is just a time series of
# returns, we directly sample from the pf_1 object created above. It has less
# than 10k observations, so we set replace=TRUE.

rets_M1 <- sample(pf_1, 10000, replace = TRUE)


# Create empty tibble to store results and vector containing different
# confidence levels.

results_1 <- tibble("alpha" = rep(NA_real_, 3),
                    "VaR" = rep(NA_real_, 3),
                    "ES" = rep(NA_real_, 3))

conf <- c(0.9, 0.95, 0.99)


# Computing VaR and ES. For every of the three confidence levels, the loop
# computes VaR and ES of the simulated portfolio return series. Both measures
# and the corresponding confidence level are stored in the tibble created
# before.

for (i in 1:3) {
  
  results_1$alpha[i] <- conf[i]
  results_1$VaR[i] <- VaR(rets_M1, p = conf[i])
  results_1$ES[i] <- ES(rets_M1, p = conf[i])
  
}


# Clean up environment.

rm("i")
rm("conf")
rm("rets_M1")
rm("results_1")


################################
###    (iii) VaR M2 & M4     ###
################################

# Simulate 10k net returns for M2. The mvrnorm() function generates random
# values for a multivariate Gaussian distribution, given multiple parameters.
# The first one is the number of simulations, namely 10k. The next two are the
# means and the covariance matrix, which have been derived earlier using MLE. To
# preserve the mean and covariance structure, we set empirical = TRUE.

values <- mvrnorm(10000, mu = M2$mu, Sigma = M2$sigma, empirical = TRUE) %>% 
          as.data.frame() 


# Create xts object with simulated returns from above. While not required here,
# the simulated returns of M2 are also used in question v), where they need to
# be in the form of a time series starting on 1 Jan 2013. Therefore, we create a
# corresponding vector of weekly dates, which is then used together with the
# simulated returns to create a xts object.

dates <- seq(as.Date("2013-01-01"), by = "weeks", length.out = 10000)
rets_M2 <- xts(values, dates)


# Compute weighted portfolio for M2. We use the same function and additional
# arguments as for M1.

pf_M2 <- Return.portfolio(rets_M2, 
                         weights = c(0.3, 0.7), 
                         rebalance_on = "week", 
                         geometric = TRUE)


# Simulate 10k net returns for M4. The mvdc() function computes parameters for
# the multivariate distribution of a copula, and takes three arguments. The
# first one specifies the copula type, which is Gaussian, to which we can pass
# the rho obtained earlier using MLE, and the dimension, which in our case is
# the number of stocks. The second argument specifies that the distributions of
# both marginals, which we set to normal. Third, we enter two lists, each of
# which contains the desired mean and standard deviation for one of the
# marginals. The values used here are again the MLE estimates from M2, which
# were used initially to compute the correlation parameter of the copula. The
# parameters obtained by this code will contain not just the marginal
# parameters, but also the correlation parameter.

dist_M4 <- mvdc(normalCopula(param = M4$par, dim = 2), 
                margins = c("norm","norm"),
                paramMargins = list(list(mean = M2$mu[1], sd = sqrt(M2$sigma[1,1])), 
                                    list(mean = M2$mu[2], sd = sqrt(M2$sigma[2,2]))))


# Simulate 10k net returns for M4. The rMvdc() function draws the specified
# number of sample returns, the values of which are are chosen to match the
# properties we specified above.

values <- rMvdc(10000, dist_M4) %>% 
                as.data.frame()

colnames(values) <- c("AAPL", "TSLA")


# Create xts object with simulated returns from M4. Like the returns from M2,
# this data will be needed later on in time series format. The code is the same
# used before.

dates <- seq(as.Date("2013-01-01"), by = "weeks", length.out = 10000)
rets_M4 <- xts(values, dates)


# Compute weighted portfolio for M4. We use the same function and additional
# arguments as for M1 and M2.

pf_M4 <- Return.portfolio(rets_M4, 
                          weights = c(0.3, 0.7), 
                          rebalance_on = "week", 
                          geometric = TRUE)


# Estimate 1-week VaR and ES for both models at confidence levels of 90%, 95%
# and 99%. The looping approach is same used for M1.

conf <- c(0.9, 0.95, 0.99)

results_M2 <- tibble("alpha" = rep(NA_real_, 3),
                     "VaR" = rep(NA_real_, 3),
                     "ES" = rep(NA_real_, 3))

results_M4 <- tibble("alpha" = rep(NA_real_, 3),
                    "VaR" = rep(NA_real_, 3),
                    "ES" = rep(NA_real_, 3))

for (i in 1:3) {
  
  results_M2$alpha[i] <- conf[i]
  results_M2$VaR[i] <- VaR(pf_M2, conf[i])
  results_M2$ES[i] <- ES(pf_M2, conf[i])
  
  results_M4$alpha[i] <- conf[i]
  results_M4$VaR[i] <- VaR(pf_M4, conf[i])
  results_M4$ES[i] <- ES(pf_M4, conf[i])
  
}


rm("i")
rm("conf")
rm("dist_M4")
rm("values")
rm("dates")


################################
###        (iv) VaR N        ###
################################

# Compute VaR and ES for last 100, 200, 300 and 400 observations of empirical
# distribution. No confidence level was specified, so we did computed it for
# every level between 80%-99%. Since no specific data was mentioned either, we
# computed both the portfolio and the separate stocks. First, we compute the
# portfolio returns using the same approach as before.

pf_EMP <- Return.portfolio(rets, 
                           weights = c(0.3, 0.7), 
                           rebalance_on = "week", 
                           geometric = TRUE)


# Initialize tibbles for results to be stored.

results_EMP <- tibble("VaR" = rep(NA_real_, 80),
                      "ES" = rep(NA_real_, 80),
                      "N" = rep(c(100, 200, 300, 400), 20),
                      "alpha" = rep(NA_real_, 80))

results_AAPL <- tibble("VaR" = rep(NA_real_, 80),
                       "ES" = rep(NA_real_, 80),
                       "N" = rep(c(100, 200, 300, 400), 20),
                       "alpha" = rep(NA_real_, 80))

results_TSLA <- tibble("VaR" = rep(NA_real_, 80),
                       "ES" = rep(NA_real_, 80),
                       "N" = rep(c(100, 200, 300, 400), 20),
                       "alpha" = rep(NA_real_, 80))


# Use loop to compute risk measures for different return series (portfolio,
# AAPL, TSLA), confidence levels (100-i%), and subsets of each series (100*v).

for(i in 1:20) {
  
  results_EMP$alpha[c(1:4)+(i-1)*4] <- (1-(i/100))
  
  results_AAPL$alpha[c(1:4)+(i-1)*4] <- (1-(i/100))
  
  results_TSLA$alpha[c(1:4)+(i-1)*4] <- (1-(i/100))
  
  for(v in 1:4) {
    
    results_EMP$VaR[v+(i-1)*4] <- VaR(last(pf_EMP, 100*v), (1-(i/100)))
    results_EMP$ES[v+(i-1)*4] <- ES(last(pf_EMP, 100*v), (1-(i/100)))
    
    results_AAPL $VaR[v+(i-1)*4] <- VaR(last(rets[, "AAPL"], 100*v), (1-(i/100)))
    results_AAPL $ES[v+(i-1)*4] <- ES(last(rets[, "AAPL"], 100*v), (1-(i/100)))
    
    results_TSLA$VaR[v+(i-1)*4] <- VaR(last(rets[, "TSLA"], 100*v), (1-(i/100)))
    results_TSLA$ES[v+(i-1)*4] <- ES(last(rets[, "TSLA"], 100*v), (1-(i/100)))
    
  }
  
}



# Plot results.

ggplot(results_EMP, aes(x = alpha, group = N)) + 
  geom_line(aes(y = VaR, color = N), linetype = "dashed") + 
  geom_line(aes(y = ES, color = N))



################################
###  (v) VaR Rolling Window  ###
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



###########################################################################
###########################################################################
###                                                                     ###
###                                PLOTS                                ###
###                                                                     ###
###########################################################################
###########################################################################




################################
###    (iii) VaR M2 & M4     ###
################################


p <- ggplot(rets_M4, aes(AAPL, TSLA)) + 
          geom_point() + 
          theme(axis.text.x = element_text(size = 14), 
                axis.text.y = element_text(size = 14)) +
          xlab("Apple") + ylab("Tesla") + 
          ggtitle(expression(paste(rho, "=0.5")))

ggExtra::ggMarginal(p, type = "histogram")

pairs.panels(rets_M4)




get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

density <- get_density(rets_M4$AAPL, rets_M4$TSLA)

ggplot(rets_M4, aes(AAPL, TSLA)) + 
  geom_point() +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14)) +
  labs(x = "Apple", 
       y = "Tesla",
       title = "Net returns, Tesla vs. Apple",
       subtitle = paste("10'000 simulations from Gaussion copula", rho, "=0.48"))


