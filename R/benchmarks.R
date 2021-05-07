#' M4_benchmarks
#'
#' Function producing the univariate forecasts employing all the benchmark methods of the M4 competition.
#' Adapted from \url{https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R}
#' @param input - Input time series (numeric vector)
#' @param fh - Forecasting horizon (numeric scalar)
#'
#' @export
#' @return List containing:
#'         \itemize{
#'         \item{\code{Naive}: }{h-step forecast for the naive method (numeric vector - length h)}
#'         \item{\code{NaiveSeasonal}: }{h-step forecast for the seasonal naive method (numeric vector - length h)}
#'         \item{\code{SimpleES}: }{h-step forecast for the simple exponential smoothing method (numeric vector - length h)}
#'         \item{\code{HoltWinters}: }{h-step forecast for the Holt-Winters method (numeric vector - length h)}
#'         \item{\code{HoltWintersDamped}: }{h-step forecast for Holt-Winters damped method (numeric vector - length h)}
#'         \item{\code{Theta}: }{h-step forecast for the theta method (numeric vector - length h)}
#'         \item{\code{Combined}: }{h-step forecast for the combined method (numeric vector - length h)}
#'         \item{\code{TimeNaive}: }{Computational time for the naive method (numeric scalar)}
#'         \item{\code{TimeNaiveSeasonal}: }{Computational time for the seasonal naive method (numeric scalar)}
#'         \item{\code{TimeSimpleES}: }{Computational time for the simple exponential smoothing method (numeric scalar)}
#'         \item{\code{TimeHoltWinters}: }{Computational time for the Holt-Winters method (numeric scalar)}
#'         \item{\code{TimeHoltWintersDamped}: }{Computational time for Holt-Winters damped method (numeric scalar)}
#'         \item{\code{TimeTheta}: }{Computational time for the theta method (numeric scalar)}
#'         \item{\code{TimeCombined}: }{Computational time for the combined method (numeric scalar)}
#'         }
#' @examples
#' x <- AirPassengers
#' splitting_point <- round(2*length(x)/3)
#' x_train <- x[1:splitting_point]
#' h <- 5
#' univariate_results <- M4_benchmarks(x_train,h)
M4_benchmarks <- function(input, fh){
  #Used to estimate the statistical benchmarks of the M4 competition

  ptm <- proc.time()
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

#' multivariate_M4benchmarks
#'
#'  Function producing the multivariate forecasts employing all the benchmark methods of the M4 competition.
#'  The multivariate implementation is made by launching in parallel the different univariate M4 benchmarks
#'  \url{https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R} on all the multivariate series
#'
#' @param Xtr - Training data (numeric matrix with N columns)
#' @param h - Forecasting horizon
#'
#' @export
#' @return List containing:
#'         \itemize{
#'         \item{\code{Naive}: }{h-step forecast for the naive method (numeric matrix - dimensions h x N)}
#'         \item{\code{NaiveSeasonal}: }{h-step forecast for the seasonal naive method (numeric matrix - dimensions h x N)}
#'         \item{\code{SimpleES}: }{h-step forecast for the simple exponential smoothing method (numeric matrix - dimensions h x N)}
#'         \item{\code{HoltWinters}: }{h-step forecast for the Holt-Winters method (numeric matrix - dimensions h x N)}
#'         \item{\code{HoltWintersDamped}: }{h-step forecast for Holt-Winters damped method (numeric matrix - dimensions h x N)}
#'         \item{\code{Theta}: }{h-step forecast for the theta method (numeric matrix - dimensions h x N)}
#'         \item{\code{Combined}: }{h-step forecast for the combined method (numeric matrix - dimensions h x N)}
#'         \item{\code{TimeNaive}: }{Computational time for the naive method (numeric scalar)}
#'         \item{\code{TimeNaiveSeasonal}: }{Computational time for the seasonal naive method (numeric scalar)}
#'         \item{\code{TimeSimpleES}: }{Computational time for the simple exponential smoothing method (numeric scalar)}
#'         \item{\code{TimeHoltWinters}: }{Computational time for the Holt-Winters method (numeric scalar)}
#'         \item{\code{TimeHoltWintersDamped}: }{Computational time for Holt-Winters damped method (numeric scalar)}
#'         \item{\code{TimeTheta}: }{Computational time for the theta method (numeric scalar)}
#'         \item{\code{TimeCombined}: }{Computational time for the combined method (numeric scalar)}
#'         }
#' @examples
#' X <- EuStockMarkets
#' splitting_point <- round(2*dim(X)[1]/3)
#' X_train <- X[1:splitting_point,]
#' h <- 5
#' multivariate_results <- multivariate_M4benchmarks(X_train,h)
multivariate_M4benchmarks <- function(Xtr,h){

  forecast_list <- apply(Xtr, 2, M4_benchmarks,fh=h)

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

#' SeasonalityTest
#'
#' Auxiliary function for computing a statistical test for seasonality
#' From \url{https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R}
#' @param input - Input time series
#' @param ppy - Seasonality period
#'
#' @return Boolean value describing whether the time series is seasonal or not
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

#' multiplicativeSeasonalityDecomposition
#'
#' Auxiliary function for computing the multiplicative seasonality decomposition
#' From \url{https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R}
#' @param input - Input time series
#' @param fh - Forecasting horizon
#'
#' @return \itemize{
#'         \item{\code{des_input}: }{Deseasonalized input time series}
#'         \item{\code{SIOut}: }{Additional deseasonalization parameters}
#'         }
#'
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

#' naive_seasonal
#'
#' Auxiliary function for the seasonal naive forecasting method.
#' From \url{https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R}
#'
#' @param input - Input time series (numeric vector)
#' @param fh - Forecasting horizon (numeric)
#' @param level - Numeric vector (length 2) containing the upper and lower bound for interval forecasting
#'
#' @import forecast
#' @return h-step forecast for the seasonal naive method (numeric vector - length h)
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


#' Theta.classic
#'
#' Auxiliary function for the Theta classic forecasting method.
#' From \url{https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R}
#' @param input - Input time series (numeric vector)
#' @param fh - Forecasting horizon (numeric scalar)
#'
#' @import forecast
#'
#' @return \itemize{
#'         \item{\code{fitted}: }{Fitted values for the ensemble method averaging \eqn{\theta_0} and \eqn{theta_2}}
#'         \item{\code{fitted0}: }{Fitted values for \eqn{\theta_0}}
#'         \item{\code{fitted2}: }{Fitted values for \eqn{\theta_2}}
#'         \item{\code{mean}: }{Forecast values for the ensemble method averaging \eqn{\theta_0} and \eqn{\theta_2}}
#'         \item{\code{mean0}: }{Forecast values for \eqn{\theta_0}}
#'         \item{\code{mean2}: }{Forecast values for \eqn{\theta_0}}
#'         }
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

#' naiveBenchmark
#'
#' Auxiliary function for the naive forecasting method.
#' From \url{https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R}
#' @param input - Input time series (numeric vector)
#' @param h - Forecasting horizon (numeric scalar)
#' @param level - Numeric vector (length 2) containing the upper and lower bound for interval forecasting
#'
#' @import forecast
#' @export
#' @return h-step forecast for the naive forecasting method (numeric vector - length h)
#' @examples
#' x <- AirPassengers
#' splitting_point <- round(2*length(x)/3)
#' x_train <- x[1:splitting_point]
#' h <- 5
#' x_hat <- naiveBenchmark(x_train,h)
naiveBenchmark <- function(input,h,level=c(80,95)){
  input_check_univariate(input,h)
  input_check_level(level)

  ptm <- proc.time()
  model <- forecast::naive(input, h=h, level=level)
  forecasts <- model$mean #Naive
  time<- proc.time() - ptm
  return(list(model=model,forecasts=forecasts,time=time[3]))
}

#' naiveSeasonalBenchmark
#'
#' Auxiliary function for the seasonal naive forecasting method.
#' From \url{https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R}
#' @param input - Input time series (numeric vector)
#' @param h - Forecasting horizon (numeric scalar)
#' @param level - Numeric vector (length 2) containing the upper and lower bound for interval forecasting
#'
#' @export
#' @return h-step forecast for the naive seasonal forecasting method (numeric vector - length h)
#' @examples
#' x <- AirPassengers
#' splitting_point <- round(2*length(x)/3)
#' x_train <- x[1:splitting_point]
#' h <- 5
#' x_hat <- naiveSeasonalBenchmark(x_train,h)
naiveSeasonalBenchmark <- function(input,h,level=c(80,95)){
  input_check_univariate(input,h)
  input_check_level(level)

  ptm <- proc.time()
  forecasts <- naive_seasonal(input, fh=h, level=level) #Seasonal Naive
  time<- proc.time() - ptm
  return(list(model=forecasts[1],forecasts=forecasts,time=time[3]))
}

#' naive2Benchmark
#'
#' Auxiliary function for the naive 2 forecasting method.
#' From \url{https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R}
#' @param input - Input time series (numeric vector)
#' @param h - Forecasting horizon (numeric scalar)
#' @param level - Numeric vector (length 2) containing the upper and lower bound for interval forecasting
#'
#' @import forecast
#' @export
#' @return h-step forecast for the naive 2 forecasting method (numeric vector - length h)
#' @examples
#' x <- AirPassengers
#' splitting_point <- round(2*length(x)/3)
#' x_train <- x[1:splitting_point]
#' h <- 5
#' x_hat <- naive2Benchmark(x_train,h)
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

#' ESBenchmark
#'
#' Auxiliary function for the simple exponential smoothing forecasting method.
#' From \url{https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R}
#' @param input - Input time series (numeric vector)
#' @param h - Forecasting horizon (numeric scalar)
#' @param level - Numeric vector (length 2) containing the upper and lower bound for interval forecasting
#'
#' @import forecast
#' @export
#' @return h-step forecast for the simple exponential smoothing forecasting method (numeric vector - length h)
#' @examples
#' x <- AirPassengers
#' splitting_point <- round(2*length(x)/3)
#' x_train <- x[1:splitting_point]
#' h <- 5
#' x_hat <- ESBenchmark(x_train,h)
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

#' HoltWintersBenchmark
#'
#' Auxiliary function for the exponential smoothing forecasting method.
#' From \url{https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R}
#' @param input - Input time series (numeric vector)
#' @param h - Forecasting horizon (numeric scalar)
#' @param level - Numeric vector (length 2) containing the upper and lower bound for interval forecasting
#'
#' @import forecast
#' @export
#' @return h-step forecast for the Holt-Winters classic method (numeric vector - length h)
#' @examples
#' x <- AirPassengers
#' splitting_point <- round(2*length(x)/3)
#' x_train <- x[1:splitting_point]
#' h <- 5
#' x_hat <- HoltWintersBenchmark(x_train,h)
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

#' HoltWintersDampedBenchmark
#'
#' Auxiliary function for the Holt-Winters damped forecasting method.
#' From \url{https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R}
#' @param input - Input time series (numeric vector)
#' @param h - Forecasting horizon (numeric scalar)
#' @param level - Numeric vector (length 2) containing the upper and lower bound for interval forecasting
#'
#' @import forecast
#' @export
#' @return h-step forecast for the theta classic method (numeric vector - length h)
#' @examples
#' x <- AirPassengers
#' splitting_point <- round(2*length(x)/3)
#' x_train <- x[1:splitting_point]
#' h <- 5
#' x_hat <- HoltWintersDampedBenchmark(x_train,h)
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

#' thetaBenchmark
#'
#' Auxiliary function for the Theta forecasting method.
#' From \url{https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R}
#' @param input - Input time series (numeric vector)
#' @param h - Forecasting horizon (numeric scalar)
#'
#' @export
#' @return h-step forecast for the Theta method (numeric vector - length h)
#' @examples
#' x <- AirPassengers
#' splitting_point <- round(2*length(x)/3)
#' x_train <- x[1:splitting_point]
#' h <- 5
#' x_hat <- thetaBenchmark(x_train,h)
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

#' combinedBenchmark
#'
#' Auxiliary function for the combined forecasting method.
#' From \url{https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R}
#' @param input - Input time series (numeric vector)
#' @param h - Forecasting horizon (numeric scalar)
#' @param level - Numeric vector (length 2) containing the upper and lower bound for interval forecasting
#'
#' @import forecast
#' @export
#' @return h-step forecast for the combined forecasting method (numeric vector - length h)
#' @examples
#' x <- AirPassengers
#' splitting_point <- round(2*length(x)/3)
#' x_train <- x[1:splitting_point]
#' h <- 5
#' x_hat <- combinedBenchmark(x_train,h)
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
