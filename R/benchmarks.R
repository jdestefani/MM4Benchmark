#' Fucntion producing the forecasts according to the benchmark methods of the M4 competition.
#' Adapted from https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R
#' @param input - Input time series
#' @param fh - Forecasting horizon
#'
#' @export
#' @return List containing the forecast and the computational times of the different methods (Naive, Naive Seasonal, Naive 2, Simple Exponential smoothing ,Holt Winters, Holt Winters Damped, Theta, Comb)
#'
M4_benchmarks <- function(input, fh){
  #Used to estimate the statistical benchmarks of the M4 competition

  ptm <- proc.time()
  #decomposition_results <- multiplicativeSeasonalityDecomposition(input)
  ppy <- stats::frequency(input)
  ST <- F
  if (ppy>1){ ST <- SeasonalityTest(input,ppy) }
  if (ST==T){
    Dec <- stats::decompose(input,type="multiplicative")
    des_input <- input/Dec$seasonal
    SIout <- utils::head(rep(Dec$seasonal[(length(Dec$seasonal)-ppy+1):length(Dec$seasonal)], fh), fh)
  }else{
    des_input <- input
    SIout <- rep(1, fh)
  }
  time_decompose <- proc.time() - ptm

  #des_input <- decomposition_results$des_input
  #SIout <- decomposition_results$SIout

  ptm <- proc.time()
  f1 <- forecast::naive(input, h=fh)$mean #Naive
  time_naive <- proc.time() - ptm

  ptm <- proc.time()
  f2 <- naive_seasonal(input, fh=fh) #Seasonal Naive
  time_naive_seasonal <- proc.time() - ptm

  ptm <- proc.time()
  f3 <- forecast::naive(des_input, h=fh)$mean*SIout #Naive2
  time_naive_2 <- proc.time() - ptm

  ptm <- proc.time()
  f4 <- forecast::ses(des_input, h=fh)$mean*SIout #Ses
  time_ses <- proc.time() - ptm

  ptm <- proc.time()
  f5 <- forecast::holt(des_input, h=fh, damped=F)$mean*SIout #Holt
  time_holt <- proc.time() - ptm

  ptm <- proc.time()
  f6 <- forecast::holt(des_input, h=fh, damped=T)$mean*SIout #Damped
  time_holt_damped <- proc.time() - ptm

  ptm <- proc.time()
  f7 <- Theta.classic(input=des_input, fh=fh)$mean*SIout #Theta
  time_theta <- proc.time() - ptm

  f8 <- (f4+f5+f6)/3 #Comb

  return(list(Naive=f1,NaiveSeasonal=f2,Naive2=f3,SimpleES=f4,HoltWinters=f5,HoltWintersDamped=f6,Theta=f7,Comb=f8,
              TimeNaive=time_naive[3],TimeNaiveSeasonal=time_naive_seasonal[3],TimeNaive2=time_naive_2[3]+time_decompose[3],
              TimeSimpleES=time_ses[3]+time_decompose[3],TimeHoltWinters=time_holt[3]+time_decompose[3],TimeHoltWintersDamped=time_holt_damped[3]+time_decompose[3],
              TimeTheta=time_theta[3]+time_decompose[3],TimeComb=time_ses[3]+time_holt[3]+time_holt_damped[3]+time_decompose[3]))
}

#' Quick and dirty multivariate implementation of the M4 benchmarks
#' https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R
#'
#' @param Xtr - Training data
#' @param H - Forecasting horizon
#'
#' @export
#' @return List of forecasts and excecution time for the different methods (Naive, Naive Seasonal, Naive 2, Simple Exponential smoothing ,Holt Winters, Holt Winters Damped, Theta, Comb)
#'         Forecast are in a HxN(number of variables) matrix format, whereas the time is in a numeric format
multivariate_M4benchmarks <- function(Xtr,H){

  forecast_list <- apply(Xtr, 2, M4_benchmarks,fh=H)

  naive_matrix <- Reduce(cbind,lapply(forecast_list, function(x){matrix(x$Naive)}))
  naive_seasonal_matrix <- Reduce(cbind,lapply(forecast_list, function(x){matrix(x$NaiveSeasonal)}))
  naive_2_matrix <- Reduce(cbind,lapply(forecast_list, function(x){matrix(x$Naive2)}))
  simple_es_matrix <- Reduce(cbind,lapply(forecast_list, function(x){matrix(x$SimpleES)}))
  holt_winters_matrix <- Reduce(cbind,lapply(forecast_list, function(x){matrix(x$HoltWinters)}))
  holt_winters_damped_matrix <- Reduce(cbind,lapply(forecast_list, function(x){matrix(x$HoltWintersDamped)}))
  theta_matrix <- Reduce(cbind,lapply(forecast_list, function(x){matrix(x$Theta)}))
  comb_matrix <- Reduce(cbind,lapply(forecast_list, function(x){matrix(x$Comb)}))

  naive_time <- Reduce('+',lapply(forecast_list, function(x){x$TimeNaive}))
  naive_seasonal_time <- Reduce('+',lapply(forecast_list, function(x){x$TimeNaiveSeasonal}))
  naive_2_time <- Reduce('+',lapply(forecast_list, function(x){x$TimeNaive2}))
  simple_es_time <- Reduce('+',lapply(forecast_list, function(x){x$TimeSimpleES}))
  holt_winters_time <- Reduce('+',lapply(forecast_list, function(x){x$TimeHoltWinters}))
  holt_winters_damped_time <- Reduce('+',lapply(forecast_list, function(x){x$TimeHoltWintersDamped}))
  theta_time <- Reduce('+',lapply(forecast_list, function(x){x$TimeTheta}))
  comb_time <- Reduce('+',lapply(forecast_list, function(x){x$TimeComb}))

  return(list(Naive=naive_matrix,NaiveSeasonal=naive_seasonal_matrix,Naive2=naive_2_matrix,SimpleES=simple_es_matrix,
              HoltWinters=holt_winters_matrix,HoltWintersDamped=holt_winters_damped_matrix,Theta=theta_matrix,Comb=comb_matrix,
              TimeNaive=naive_time,TimeNaiveSeasonal=naive_seasonal_time,TimeNaive2=naive_seasonal_time,
              TimeSimpleES=simple_es_time,TimeHoltWinters=holt_winters_time,TimeHoltWintersDamped=holt_winters_damped_time,
              TimeTheta=theta_time,TimeComb=comb_time))

}

#' Auxiliary function for computing a statistical test for seasonality
#' From https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R
#' @param input - Input time series
#' @param ppy - Seasonality period
#'
#' @return test_seasonal <- Boolean value describing whether the time series is seasonal or not
SeasonalityTest <- function(input, ppy){
  #Used to determine whether a time series is seasonal
  tcrit <- 1.645
  if (length(input)<3*ppy){
    test_seasonal <- FALSE
  }else{
    xacf <- stats::acf(input, plot = FALSE)$acf[-1, 1, 1]
    clim <- tcrit/sqrt(length(input)) * sqrt(cumsum(c(1, 2 * xacf^2)))
    test_seasonal <- ( abs(xacf[ppy]) > clim[ppy] )

    if (is.na(test_seasonal)==TRUE){ test_seasonal <- FALSE }
  }

  return(test_seasonal)
}

#' Auxiliary function for computing the multiplicative seasonality decomposition
#' From https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R
#' @param input - Input time series
#' @param fh - Forecasting horizon
#'
#' @return des_input <- Deseasonalized input time series
#'         SIout <- Additional deseasonalization parameters
multiplicativeSeasonalityDecomposition <- function(input,fh){
  ppy <- stats::frequency(input)
  ST <- F
  if (ppy>1){ ST <- SeasonalityTest(input,ppy) }
  if (ST==T){
    Dec <- stats::decompose(input,type="multiplicative")
    des_input <- input/Dec$seasonal
    SIout <- utils::head(rep(Dec$seasonal[(length(Dec$seasonal)-ppy+1):length(Dec$seasonal)], fh), fh)
  }else{
    des_input <- input
    SIout <- rep(1, fh)
  }
  return(list(des_input=des_input,SIout=SIout))
}

#'
#' Auxiliary function for the seasonal naive forecasting method.
#' From https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R
#'
#' @param input - Input time series (numeric vector)
#' @param fh - Forecasting horizon (numeric)
#' @param level - Numeric vector (length 2) containing the upper and lower bound for interval forecasting
#'
#' @import forecast
#' @export
#' @return h-step forecast for the seasonal naive method
naive_seasonal <- function(input, fh, level=c(80,95)){
  input_check_univariate(input,fh)
  input_check_level(level)

  #Used to estimate Seasonal Naive
  frcy <- stats::frequency(input)
  frcst <- forecast::naive(input, h=fh,level=level)$mean
  if (frcy>1){
    frcst <- utils::head(rep(as.numeric(utils::tail(input,frcy)), fh), fh) + frcst - frcst
  }
  return(frcst)
}


#' Auxiliary function for the Theta classic forecasting method.
#' From https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R
#' @param input - Input time series
#' @param fh - Forecasting horizon
#'
#' @import forecast
#' @export
#' @return h-step forecast for the theta classic method
Theta.classic <- function(input, fh){
  input_check_univariate(input,fh)
  #Used to estimate Theta classic

  #Set parameters
  wses <- wlrl<-0.5
  theta <- 2
  #Estimate theta line (0)
  observations <- length(input)
  xt <- c(1:observations)
  xf <- c((observations+1):(observations+fh))
  train <- data.frame(input=input, xt=xt)
  test <- data.frame(xt = xf)

  estimate <- stats::lm(input ~ poly(xt, 1, raw=TRUE))
  thetaline0In <- as.numeric(stats::predict(estimate))
  thetaline0Out <- as.numeric(stats::predict(estimate,test))

  #Estimate theta line (2)
  thetalineT <- theta*input+(1-theta)*thetaline0In
  sesmodel <- forecast::ses(thetalineT, h=fh)
  thetaline2In <- sesmodel$fitted
  thetaline2Out <- sesmodel$mean

  #Theta forecasts
  forecastsIn <- (thetaline2In*wses)+(thetaline0In*wlrl)
  forecastsOut <- (thetaline2Out*wses)+(thetaline0Out*wlrl)

  #Zero forecasts become positive
  for (i in 1:length(forecastsOut)){
    if (forecastsOut[i]<0){ forecastsOut[i]<-0 }
  }

  output=list(fitted = forecastsIn, mean = forecastsOut,
              fitted0 = thetaline0In, mean0 = thetaline0Out,
              fitted2 = thetaline2In, mean2 = thetaline2Out)

  return(output)

}

#' Auxiliary function for the naive forecasting method.
#' From https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R
#' @param input - Input time series
#' @param h - Forecasting horizon
#' @param level - Numeric vector (length 2) containing the upper and lower bound for interval forecasting
#'
#' @import forecast
#' @export
#' @return h-step forecast for the theta classic method

naiveBenchmark <- function(input,h,level=c(80,95)){
  input_check_univariate(input,h)
  input_check_level(level)

  ptm <- proc.time()
  model <- forecast::naive(input, h=h, level=level)
  forecasts <- model$mean #Naive
  time<- proc.time() - ptm
  return(list(model=model,forecasts=forecasts,time=time[3]))
}

#' Auxiliary function for the seasonal naive forecasting method.
#' From https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R
#' @param input - Input time series
#' @param h - Forecasting horizon
#' @param level - Numeric vector (length 2) containing the upper and lower bound for interval forecasting
#'
#' @export
#' @return h-step forecast for the theta classic method
naiveSeasonalBenchmark <- function(input,h,level=c(80,95)){
  input_check_univariate(input,h)
  input_check_level(level)

  ptm <- proc.time()
  forecasts <- naive_seasonal(input, fh=h, level=level) #Seasonal Naive
  time<- proc.time() - ptm
  return(list(model=forecasts[1],forecasts=forecasts,time=time[3]))
}

#' Auxiliary function for the naive 2 forecasting method.
#' From https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R
#' @param input - Input time series
#' @param h - Forecasting horizon
#' @param level - Numeric vector (length 2) containing the upper and lower bound for interval forecasting
#'
#' @import forecast
#' @export
#' @return h-step forecast for the theta classic method
naive2Benchmark <- function(input,h,level=c(80,95)){
  input_check_univariate(input,h)
  input_check_level(level)

  ptm <- proc.time()
  decomposition_results <- multiplicativeSeasonalityDecomposition(input,h)
  des_input <- decomposition_results$des_input
  SIout <- decomposition_results$SIout

  model <- forecast::naive(des_input, h=h, level=level)
  forecasts <- model$mean*SIout #Naive 2
  time<- proc.time() - ptm
  return(list(model=model,forecasts=forecasts,time=time[3]))
}

#' Auxiliary function for the exponential smoothing forecasting method.
#' From https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R
#' @param input - Input time series
#' @param h - Forecasting horizon
#' @param level - Numeric vector (length 2) containing the upper and lower bound for interval forecasting
#'
#' @import forecast
#' @export
#' @return h-step forecast for the theta classic method
ESBenchmark <- function(input,h,level=c(80,95)){
  input_check_univariate(input,h)
  input_check_level(level)

  ptm <- proc.time()
  decomposition_results <- multiplicativeSeasonalityDecomposition(input,h)
  des_input <- decomposition_results$des_input
  SIout <- decomposition_results$SIout

  model <- forecast::ses(des_input, h=h,level = level)
  forecasts <- model$mean*SIout # Exponential smoothing
  time<- proc.time() - ptm
  return(list(model=model,forecasts=forecasts,time=time[3]))
}

#' Auxiliary function for the exponential smoothing forecasting method.
#' From https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R
#' @param input - Input time series
#' @param h - Forecasting horizon
#' @param level - Numeric vector (length 2) containing the upper and lower bound for interval forecasting
#'
#' @import forecast
#' @export
#' @return h-step forecast for the Holt-Winters classic method
HoltWintersBenchmark <- function(input,h,level=c(80,95)){
  input_check_univariate(input,h)
  input_check_level(level)

  ptm <- proc.time()
  decomposition_results <- multiplicativeSeasonalityDecomposition(input,h)
  des_input <- decomposition_results$des_input
  SIout <- decomposition_results$SIout

  model <- forecast::holt(des_input, h=h, damped=F, level = level)
  forecasts <- model$mean*SIout #Holt
  time<- proc.time() - ptm
  return(list(model=model,forecasts=forecasts,time=time[3]))
}

#' Auxiliary function for the Holt-Winters damped forecasting method.
#' From https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R
#' @param input - Input time series
#' @param h - Forecasting horizon
#' @param level - Numeric vector (length 2) containing the upper and lower bound for interval forecasting
#'
#' @import forecast
#' @export
#' @return h-step forecast for the theta classic method
HoltWintersDampedBenchmark <- function(input,h,level=c(80,95)){
  input_check_univariate(input,h)
  input_check_level(level)

  ptm <- proc.time()
  decomposition_results <- multiplicativeSeasonalityDecomposition(input,h)
  des_input <- decomposition_results$des_input
  SIout <- decomposition_results$SIout

  model <- forecast::holt(des_input, h=h, damped=T,level=level)
  forecasts <- model$mean*SIout #Damped
  time<- proc.time() - ptm
  return(list(model=model,forecasts=forecasts,time=time[3]))
}

#' Auxiliary function for the theta forecasting method.
#' From https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R
#' @param input - Input time series
#' @param h - Forecasting horizon
#'
#' @export
#' @return h-step forecast for the theta classic method
thetaBenchmark <- function(input,h){
  input_check_univariate(input,h)

  ptm <- proc.time()
  decomposition_results <- multiplicativeSeasonalityDecomposition(input,h)
  des_input <- decomposition_results$des_input
  SIout <- decomposition_results$SIout

  model <- Theta.classic(input=des_input, fh=h)
  forecasts <- model$mean*SIout #Theta
  time<- proc.time() - ptm
  return(list(model=model,forecasts=forecasts,time=time[3]))
}

#' Auxiliary function for the combined forecasting method.
#' From https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R
#' @param input - Input time series
#' @param h - Forecasting horizon
#' @param level - Numeric vector (length 2) containing the upper and lower bound for interval forecasting
#'
#' @import forecast
#' @export
#' @return h-step forecast for the theta classic method
combinedBenchmark <- function(input,h,level=c(80,95)){
  input_check_univariate(input,h)
  input_check_level(level)

  ptm <- proc.time()
  decomposition_results <- multiplicativeSeasonalityDecomposition(input,h)
  des_input <- decomposition_results$des_input
  SIout <- decomposition_results$SIout

  models <- list(forecast::ses(des_input, h=h, level=level),
                 forecast::holt(des_input, h=h, damped=F, level=level),
                 forecast::holt(des_input, h=h, damped=T, level=level))
  forecasts <- (models[[1]]$mean*SIout + models[[2]]$mean*SIout + models[[3]]$mean*SIout)/3

  time<- proc.time() - ptm
  return(list(model=models,forecasts=forecasts,time=time[3]))
}
