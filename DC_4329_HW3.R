library(tseries)
library(forecast)

setwd("~/Downloads/4329_Materials/HW3/data")

bond_data <- read.csv("data.csv")
arrests <- read.csv("arrests.csv")
hitters <- read.csv("hitters.csv")


## bond data 
  # column r is the quarterly t-bill rate,
  # column y the quarterly log GDP 
  # column pi the quarterly inflation rate.

# plot time series
plot(ts(bond_data))

# plot autocorrelation functions
acf(bond_data$r)
acf(bond_data$y)
acf(bond_data$pi)

# ADF test on 3 time series
adf_stats <- apply(ts(bond_data), MARGIN = 2, adf.test)

# Phillips-Perron Unit Root Test on 3 time series
pp_stats <- apply(ts(bond_data), MARGIN = 2, pp.test)

# (e) Describe the signs of non stationarity seen in the time series and ACF plots.
## Significant level of correlation on t-bill rate and quarterly log of GDP
## Lowest oscilating value of autocorrelation by lags can be observed in inflation rate
## thus indicating that r and y are non-stationary, and pi is potentially statiorary
  

# (f) Use the ADF tests to decide which of the series are nonstationary. Do the tests
# corroborate the conclusions of the time series and ACF plots?
adf_stats$r
pp_stats$r
# adf stats indicate that quarterly t-bill rate is non-stationary, H0 accepted
# this result corraborates with the ACF interpretation

adf_stats$y
pp_stats$y
# adf stats indicate that quarterly log GDP is non-stationary, H0 accepted
# this result corraborates with the ACF interpretation


adf_stats$pi
pp_stats$pi
kpss.test(ts(bond_data$pi))
# in case of inflation based on ADF test, null hypothesis is accepted as non-stationary 
# but due to low p-value, additional tests PP reject the H0 and indicate that inflation is stationary
# while KPSS rejects the null hypothesis of inflation being stationary. 
## hence the unit root tests are somewhat contradictory indicating the possibility that 
## the inflation rates are stationary with long-term memory.


# (g) Run the ADF test on the differenced series and plot them.
# Are they stationary?
d_bond_data = diff(ts(bond_data))
plot(d_bond_data)

adf_stats <- apply(d_bond_data, MARGIN = 2, adf.test)
## based on the ADF test all 3 time series are stationary
pp_stats <- apply(d_bond_data, MARGIN = 2, pp.test)

## this corraborates with pp test
kpss.test(ts(d_bond_data[,3]))


# (h) Do you see evidence of autocorrelation in the differenced series?
#  If so, describe these correlations.
acf(d_bond_data)

###############################################################
# Fit an automatic ARIMA model to the r and y series
arima_r <- auto.arima(d_bond_data[,1], max.P=0, max.Q=0,ic="aic")
arima_y <- auto.arima(d_bond_data[,2], max.P=0, max.Q=0,ic="aic")

# plot check 
plot(arima_r$fitted)
lines(d_bond_data[,1],col="green")

#plot residuals
plot(arima_r$residuals)

acf(arima_r$residuals) # autocorrelation persists for residuals
Box.test(arima_r$residuals, lag = 5, type="Ljung")

summary(arima_r)  #best model is ARIMA(0,0,3) = MA(3)


# plot check 
plot(arima_y$fitted)
lines(d_bond_data[,2],col="green")

#plot residuals
plot(arima_y$residuals)
acf(arima_y$residuals) 
Box.test(arima_y$residuals, lag = 5, type="Ljung") 
# null hypothesis accepted, no autocorrelation

summary(arima_y)   #best model is ARIMA(1,0,1) = ARMA(1,1)
###############################################################

###############################################################
# Fit an automatic ARIMA model to the r and y series
arima_r <- auto.arima(d_bond_data[,1], max.P=0, max.Q=0,ic="bic")
arima_y <- auto.arima(d_bond_data[,2], max.P=0, max.Q=0,ic="bic")

# plot check 
plot(arima_r$fitted)
lines(d_bond_data[,1],col="green")

#plot residuals
plot(arima_r$residuals)
summary(arima_r)  #best model is ARIMA(0,0,3) = MA(3)


# plot check 
plot(arima_y$fitted)
lines(d_bond_data[,2],col="green")

#plot residuals
plot(arima_y$residuals)
summary(arima_y)   #best model is ARIMA(1,0,1) = AR(1), best model changed
###############################################################

# GARCH effects
resid2 = (residuals(arima_r) - mean(residuals(arima_r)))^2
plot(resid2)
acf(resid2)
Box.test(resid2, lag = 5, type="Ljung")

resid2_y = (residuals(arima_y) - mean(residuals(arima_y)))^2
plot(resid2_y)
acf(resid2_y)
Box.test(resid2_y, lag = 5, type="Ljung")

# therefore there is GARCH effects becasue the mean centre sqared residual series has serial
# autocorrelation as shown by LB test.

#################################################################

var_model <- ar(d_bond_data,order.max=4, aic=T)
var_res <- residuals(var_model) 
var_res <- na.omit(var_res) # omit NA
plot(var_res)

# residual  of the model are correlated
cor(var_res)
cor(var_res) - cor(d_bond_data) 
# residuals are correlated to the extent that underlying data is correlated

acf(var_res[,1]) 
Box.test(var_res[,1], lag = 5, type="Ljung") # no autocorrelation for r
acf(var_res[,2])
Box.test(var_res[,2], lag = 5, type="Ljung") # no autocorrelation for y
acf(var_res[,3])
Box.test(var_res[,3], lag = 5, type="Ljung") # no autocorrelation for pi

# c) What are the coecient estimates (e.g.,Φ 1 , ..., Φ p )?
F_parameters <- var_model$ar

# estimated correlation matrix for the residuals
cov(var_res)

#Do you believe that the model fits well?

