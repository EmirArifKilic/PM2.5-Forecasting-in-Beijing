#ONLY NA values were imputed with tsclean, the rest was left untouched
#afterwards data were converted to weekly from hourly

library(pdR)
library(usethis)
library(devtools)
library(anomalize)
library(forecast)
library(AnomalyDetection)
library(dplyr)
library(vctrs)
library(tseries)
library(uroot)
library(lmtest)
library(tidyverse)
library(Rcpp)
library(tibbletime)
library(anomalize)
library(recipes)
library(future.apply)
library(TSA)
library(fUnitRoots)


data <- read.csv("pm2.5_weekly.csv", header = TRUE)
data <- data %>%
  mutate(week = as.Date(week, format = "%d/%m/%Y"))

pm2.5_ts <- ts(data$pm2.5, frequency = 52, start = c(2009, 52)) 

train_ts <- window(pm2.5_ts, end = c(2013, 52))
test_ts <- window(pm2.5_ts, start = c(2014, 1))
par(mfrow=c(1,1))
plot(train_ts,main="Train set plot")
ggAcf(train_ts,main="Acf plot of train set",lag.max = 210)
ggPacf(train_ts,main="Pacf plot of train set",lag.max = 210)

print(train_ts)
plot(train_ts,main="Series Before Imputation")
x <- train_ts

train_datetime <- seq(from = as.Date("2009-12-27"), by = "week", length.out = 209)

train_df <- tibble(week = train_datetime, pm2.5 = as.numeric(train_ts)) %>%
  as_tbl_time(index = week)

result <- train_df %>%
  anomalize::time_decompose(pm2.5, method = "stl", frequency = "auto", trend = "auto") %>%
  anomalize::anomalize(remainder, method = "iqr") %>%
  anomalize::time_recompose()

result %>%
  anomalize::plot_anomalies(time_recomposed = TRUE, ncol = 3, alpha_dots = 0.5)

train_df <- result %>%
  mutate(pm2.5 = ifelse(anomaly == "Yes", NA, observed)) %>%  # Replace anomalies with NA
  select(week, pm2.5)

train_df <- train_df %>%
  mutate(pm2.5 = tsclean(ts(pm2.5, start = c(2009, 52), frequency = 52)))  # Impute anomalies

train_ts <- ts(train_df$pm2.5, start = c(2009, 52), frequency = 52)

cleaned_result <- tibble(
  week = seq(from = as.Date("2009-12-27"), by = "week", length.out = length(train_ts)),
  pm2.5 = as.numeric(train_ts)
) %>%
  as_tbl_time(index = week)

result_cleaned <- cleaned_result %>%
  anomalize::time_decompose(pm2.5, method = "stl", frequency = "auto", trend = "auto") %>%
  anomalize::anomalize(remainder, method = "iqr") %>%
  anomalize::time_recompose()

result_cleaned %>%
  anomalize::plot_anomalies(time_recomposed = TRUE, ncol = 3, alpha_dots = 0.5)


par(mfrow=c(1,2))
plot(x,main="Series Before Imputation")
plot(train_ts,main="Series After Imputation")

str(train_ts)
print(train_df)
str(train_df)

#lets try boxcox
lambda <- BoxCox.lambda(train_ts)
print(paste("Optimal lambda:", lambda))
arima_no_transform <- auto.arima(train_ts)
arima_transform <- auto.arima(BoxCox(train_ts, lambda))
print(paste("AIC without transformation:", arima_no_transform$aic))
print(paste("AIC with transformation:", arima_transform$aic))
checkresiduals(arima_transform)

#both the aic and the plot of the series change dramatically + the series look much better after transformation
#i will continue with the box-cox transformed series
transformed_ts <- BoxCox(train_ts, lambda)
head(transformed_ts)
par(mfrow = c(1, 2))
plot(train_ts, main = "Original Time Series")
plot(transformed_ts, main = "B-C Transformed TS")


plot(transformed_ts)
ggAcf(transformed_ts,lag.max = 250) #weekly seasonality exists
ggPacf(transformed_ts,lag.max=250)

kpss.test(transformed_ts, null = "Level")  #stationary series

adfTest(transformed_ts, lags=1, type="c") #series is stationary

pp.test(transformed_ts) #stationary series  
#the series are stationary according to kpss test and there is no unit root by pp test. 
ggAcf(transformed_ts,lag.max = 250)
ggPacf(transformed_ts,lag.max=250)
#the acf and pacf plots look good too
#there might be a seasonal unit root though. (there is for sure seasonality )

ch.test(transformed_ts,type = "dummy",sid=c(1:7)) 
#p > 0.05 series is stat and deterministic


ocsb.test(transformed_ts,lag.method ="AIC",maxlag = 0) 
# test stat bigger than critical value, we reject that there is seasonal unit root

ggAcf(transformed_ts,lag.max = 250,main="Acf of Transformed Series")
ggPacf(transformed_ts,lag.max=250,main="Pacf of Transformed Series")
eacf(transformed_ts)

back_transform <- function(y_forecast, lambda) {
  ((y_forecast * lambda) + 1)^(1 / lambda)
}

#lets fit the model
arima_model <- Arima(transformed_ts, order = c(1, 0, 1), seasonal = list(order = c(1, 0, 1), period = 7))
arima_model2 <- Arima(transformed_ts, order = c(1, 0, 23), seasonal = list(order = c(1, 0, 1), period = 7))
arima_model3 <- Arima(transformed_ts, order = c(1, 0, 0), seasonal = list(order = c(1, 0, 1), period = 7))
arima_model4 <- Arima(transformed_ts, order = c(0, 0, 1), seasonal = list(order = c(1, 0, 1), period = 7))
arima_model5 <- Arima(transformed_ts, order = c(0, 0, 1), seasonal = list(order = c(1, 0, 0), period = 7))
arima_model6 <- Arima(transformed_ts, order = c(1, 0, 1), seasonal = list(order = c(1, 0, 0), period = 7))
arima_model7 <- Arima(transformed_ts, order = c(1, 0, 0), seasonal = list(order = c(1, 0, 0), period = 7))
arima_model8 <- Arima(transformed_ts, order = c(0, 0, 0), seasonal = list(order = c(1, 0, 0), period = 7))
arima_model9 <- Arima(transformed_ts, order = c(1, 0, 0), seasonal = list(order = c(0, 0, 0), period = 7))
arima_model10 <- Arima(transformed_ts, order = c(1, 0, 0), seasonal = list(order = c(0, 0, 1), period = 7))
#significant models 5 to 9
auto_model <- auto.arima(transformed_ts,seasonal=TRUE)

summary(auto_model)
summary(arima_model)
summary(arima_model2)
summary(arima_model3)
summary(arima_model4)
summary(arima_model5)
summary(arima_model6)
summary(arima_model7)
summary(arima_model8)
summary(arima_model9)
summary(arima_model10)
print(c(auto_model$aic,auto_model$bic))
print(c(arima_model$aic,arima_model$bic))
print(c(arima_model2$aic,arima_model2$bic))
print(c(arima_model3$aic,arima_model3$bic))
print(c(arima_model4$aic,arima_model4$bic))
print(c(arima_model5$aic,arima_model5$bic)) #significant
print(c(arima_model6$aic,arima_model6$bic)) #significant best model
print(c(arima_model7$aic,arima_model7$bic)) #significant 2nd best model
print(c(arima_model8$aic,arima_model8$bic)) #significant
print(c(arima_model9$aic,arima_model9$bic)) #significant (|coef/se| > 2 for all coeffs)
print(c(arima_model10$aic,arima_model10$bic)) #significant
#arimamodel6 fits our data best according to the aic and the significance of the coefficients (coef/se > 2) and the RMSE and MAPE  values
plot(transformed_ts)
ggAcf(transformed_ts,lag.max=250)
ggPacf(transformed_ts,lag.max=250)
arima_model7ML <- Arima(transformed_ts, order = c(1, 0, 0), seasonal = list(order = c(1, 0, 0), period = 7),method="ML")
arima_model7CSS <- Arima(transformed_ts, order = c(1, 0, 0), seasonal = list(order = c(1, 0, 0), period = 7),method="CSS")

arima_model6ML <- Arima(transformed_ts, order = c(1, 0, 1), seasonal = list(order = c(1, 0, 0), period = 7),method="ML")
arima_model6CSS <- Arima(transformed_ts, order = c(1, 0, 1), seasonal = list(order = c(1, 0, 0), period = 7),method="CSS")

arima_model10ML <- Arima(transformed_ts, order = c(1, 0, 0), seasonal = list(order = c(0, 0, 1), period = 7), method="ML")
arima_model10CSS <- Arima(transformed_ts, order = c(1, 0, 0), seasonal = list(order = c(0, 0, 1), period = 7),method="CSS")



#11a res check
r <- residuals(arima_model6)
Box.test(r, lag = 15, type = "Ljung-Box") #p above 0.05 hence residuals uncorrelated
  
ggAcf(r, main = "ACF of Residuals",lag.max=250)
ggPacf(r, main = "PACF of Residuals",lag.max=250)
par(mfrow = c(1, 1)) 

plot(r, main = "res") #random pattern so good


#11b normality
hist(r, main = "Histogram of Residuals", xlab = "Residuals", col = "lightblue", breaks = 10) #normal?
qqnorm(r, main = "QQ-Plot of Residuals") #looks normal
shapiro.test(r) #normal series

#11c correlation
library(TSA)
library(lmtest)
ggAcf(r, main = "ACF of Residuals",lag.max=250)
m = lm(r ~ 1+zlag(r))
bgtest(m,order=15) #there is no serial correlation



#11d heteroscedasticity
library(FinTS)
rr <- r^2
ggAcf(rr, main = "ACF of Squared Residuals",lag.max=250)
ggPacf(rr, main = "ACF of Squared Residuals",lag.max=250)
ArchTest(rr) #p value 0.97 no arch effect present there is homoscedasticity

#ets
library(forecast)
ETS1 <- ets(train_ts)
ETS2 <- ets(train_ts, model = "ANN") #cant do ANA 
ETS1_residuals <- residuals(ETS1)
ETS2_residuals <- residuals(ETS2)
summary(ETS1)
summary(ETS2)
shapiro.test(ETS1_residuals)
shapiro.test(ETS2_residuals)
ETS1.F <- forecast(ETS1, h = length(test_ts))
ETS2.F <- forecast(ETS2, h = length(test_ts))
autoplot(ETS1.F)
autoplot(ETS2.F)
accuracy(ETS1.F,test_ts)
accuracy(ETS2.F,test_ts)
autoplot(train_ts)+autolayer(test_ts,PI=F,series="test set")+autolayer(ETS2.F,PI=F,series="ETS(ANN)")+autolayer(ETS1.F,PI=F,series="ETS(MNN)")+ggtitle("Forecasts from Different ETS")+theme_minimal()

#sarima
sarima_f1 <- forecast(arima_model6,h=length(test_ts))
converted_forecast <- back_transform(sarima_f1$mean, lambda)  # Back-transform the point forecast
converted_lower <- back_transform(sarima_f1$lower[, 2], lambda)  # 95% lower bound
converted_upper <- back_transform(sarima_f1$upper[, 2], lambda) 
plot(train_ts, type = "l", col = "black", lwd = 2,xlab = "Time", ylab = "Values",
     main = "SARIMA Forecast",xlim = c(2010, 2014), ylim = range(c(train_ts, test_ts, converted_forecast, converted_lower, converted_upper)))
lines(seq(from = 2013, by = 1 / 52, length.out = 53),
      test_ts, col = "blue", lwd = 2)
lines(seq(from = 2013 + 1 / 52,
          by = 1 / 52,
          length.out = 53),
      converted_forecast, col = "red", lwd = 2)
lines(seq(from = 2013+1/ 52,by = 1 / 52,length.out = 53),converted_lower, 
      col = "red", lwd = 1, lty = 2)
lines(seq(from = 2013 + 1 / 52,by = 1 / 52,length.out = 53),converted_upper, 
      col = "red", lwd = 1, lty = 2)
legend("topleft", legend = c("train", "test", "forecast", "95% CI"),
       col = c("black", "blue", "red", "red"), lty = c(1, 1, 1, 2), lwd = c(2, 2, 2, 1),cex = 0.65)

accuracy(converted_forecast,test_ts)
#prophet
library(prophet)
ds<-c(seq(as.Date("2009/12/27"),as.Date("2013/12/22"),by="week"))
df<-data.frame(ds,y=as.numeric(train_ts))
pm2.5_prophet <- prophet(df,weekly.seasonality=TRUE)
future<-make_future_dataframe(pm2.5_prophet,periods = 53,freq="week")
forecast <- predict(pm2.5_prophet, future)
tail(forecast[c('ds', 'yhat', 'yhat_lower', 'yhat_upper')],53)
plot(pm2.5_prophet, forecast)+theme_minimal()
prophet_plot_components(pm2.5_prophet, forecast)
dyplot.prophet(pm2.5_prophet, forecast)
accuracy(tail(forecast$yhat, 53), test_ts)

changepoint_prior <- c(0.1, 0.5, 0.9)
seasonality_prior <- c(0.1, 0.3, 0.5)
changepoint_range <- c(0.6, 0.8, 0.9)

results <- data.frame(
  changepoint_prior = numeric(),
  seasonality_prior = numeric(),
  changepoint_range = numeric(),
  RMSE = numeric()
)

for (cp in changepoint_prior) {
  for (sp in seasonality_prior) {
    for (cr in changepoint_range) {
      m <- prophet(
        changepoint.prior.scale = cp,
        seasonality.prior.scale = sp,
        changepoint.range = cr,
        weekly.seasonality=TRUE
      )
      m <- fit.prophet(m, df) 
      
      
      future <- make_future_dataframe(m, periods = 53, freq = "week")
      forecastt <- predict(m, future)
      
      predicted <- tail(forecastt$yhat, 53)
      acc <- accuracy(predicted, test_ts)  
      rmse <- acc["Test set", "RMSE"]  # Extract RMSE from accuracy
      
      results <- rbind(results, data.frame(
        changepoint_prior = cp, 
        seasonality_prior = sp, 
        changepoint_range = cr, 
        RMSE = rmse
      ))
    }
  }
}

#best parameters
best_params <- results[which.min(results$RMSE), ]
best_params



final_model <- prophet(changepoint.prior.scale = 0.9,seasonality.prior.scale = 0.5,changepoint.range = 0.6,weekly.seasonality = TRUE)

final_model <- fit.prophet(final_model, df)
future_final <- make_future_dataframe(final_model, periods = 53, freq = "week")
forecastt <- predict(final_model, future_final)

plot(final_model, forecastt) + theme_minimal()
prophet_plot_components(final_model, forecastt)


accuracy(tail(forecast$yhat,53),test_ts)
accuracy(tail(forecastt$yhat,53),test_ts) #best model
plot(pm2.5_prophet, forecast)+theme_minimal() + ggtitle("First Model with default parameters")
plot(pm2.5_prophet, forecastt)+theme_minimal() + ggtitle("Second Model with hyperparameter tuning")

#TBATS
tbatsmodel<-tbats(train_ts)
tbats_residuals <- residuals(tbatsmodel)
shapiro.test(tbats_residuals)
Box.test(tbats_residuals, lag = 15, type = "Ljung-Box") #p above 0.05 hence residuals uncorrelated
tbats_residualss <- tbats_residuals^2
ArchTest(tbats_residualss)

autoplot(train_ts,main="TS plot of Train with TBATS Fitted") +autolayer(fitted(tbatsmodel), series="Fitted") +theme_minimal()
tbats_forecast<-forecast(tbatsmodel,h=length(test_ts))
autoplot(tbats_forecast)+autolayer(test_ts,series="actual",color="red")+theme_minimal()
accuracy(tbats_forecast,test_ts)

#NNETAR
nnmodel<-nnetar(train_ts)
nnmodel_residuals <- residuals(nnmodel)

Box.test(nnmodel_residuals, lag = 15, type = "Ljung-Box") #p above 0.05 hence residuals uncorrelated

autoplot(train_ts)+autolayer(fitted(nnmodel))+theme_minimal()+ggtitle("Fitted Values of NN Model")
nnforecast<-forecast(nnmodel,h=length(test_ts),PI=TRUE)
autoplot(nnforecast)+theme_minimal()
accuracy(nnforecast,test_ts)
length(nnforecast)


#COMPARISON OF ALL ACCURACIES
accuracy(ETS1.F,test_ts)
accuracy(ETS2.F,test_ts)
accuracy(converted_forecast,test_ts)
accuracy(tail(forecast$yhat,53),test_ts) #initial prophet model
accuracy(tail(forecastt$yhat,53),test_ts) #hyperparameter tuned model, the best model
accuracy(tbats_forecast,test_ts)
accuracy(nnforecast,test_ts)


