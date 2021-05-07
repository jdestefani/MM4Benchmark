# Constant string literals
error_X_type_string <- "[ERROR] - Parameter X is not a matrix"
error_x_type_string <- "[ERROR] - Parameter x is not numeric"
error_x_not_vector_string <- "[ERROR] - Parameter x is not a numeric vector"
error_h_type_string <- "[ERROR] - Parameter h is not numeric"
error_h_not_scalar_string <- "[ERROR] - Parameter x is not a numeric scalar"
error_level_type_string <- "[ERROR] - Parameter level is notnumeric"
error_level_length_string <- "[ERROR] - Parameter level is not of length 2"
error_level_order_string <- "[ERROR] - Parameter level's values should be ordered"

#' input_check_multivariate
#'
#' Verifies that the input values to the multivariate forecasting function have the correct types and dimensions
#'
#' @param X - Target matrix (hxk) for k variables
#' @param h - Forecasting horizon (numeric scalar)
#'
#' @return - NULL
input_check_multivariate <- function(X,h){
  if(!is.matrix(X)){
    stop(error_X_type_string)
  }

  if(!is.numeric(h)){
    stop(error_h_type_string)
  }

  if(length(h) !=1){
    stop(error_h_not_scalar_string)
  }


}

#' input_check_univariate
#'
#' Verifies that the input values to the univariate forecasting function have the correct types and dimensions
#'
#' @param x - Input univariate series (numeric vector)
#' @param h - Forecasting horizon (numeric scalar)
#'
#' @return - NULL
input_check_univariate <- function(x,h){
  if(!is.numeric(x)){
    stop(error_x_type_string)
  }

  if(!is.null(dim(x))){
    stop(error_x_not_vector_string)
  }

  if(!is.numeric(h)){
    stop(error_h_type_string)
  }

  if(length(h) !=1){
    stop(error_h_not_scalar_string)
  }

}

#' input_check_level
#'
#' Verifies that the level input to the univariate forecasting function has the correct types and shape
#'
#' @param level - Numeric vector (length 2) containing the upper and lower bound for interval forecasting
#'
#' @return - NULL
input_check_level <- function(level){
  if(!is.numeric(level)){
    stop(error_level_type_string)
  }

  if(length(level) != 2){
    stop(error_level_length_string)
  }

  if(level[1] > level[2]){
    stop(error_level_order_string)
  }

}

