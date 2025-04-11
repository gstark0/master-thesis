# sem_utils.R
# Contains utility functions for model fitting and predictions

calculate_mse <- function(predictions, actual) {
  return(mean((actual - predictions)^2))
}

# Update a lavaan model with new parameter values
get_new_fit <- function(new_params, fit) {
    # Create a copy of the original fit
    updated_fit <- fit

    # Update the parameter estimates in ParTable
    free_idx <- which(updated_fit@ParTable$free > 0)
    updated_fit@ParTable$est[free_idx] <- new_params

    # Update the internal model parameters
    updated_fit@Model <- lavaan:::lav_model_set_parameters(updated_fit@Model, new_params)
    
    # Ensure implied matrices are updated
    tryCatch({
        implied <- lavaan:::lav_model_implied(updated_fit@Model)
        updated_fit@implied <- implied
    }, error = function(e) {
        # If the implied calculation fails, just return the model with updated parameters
    })

    return(updated_fit)
}

# Predict outcomes (Y) based on predictors (X) using SEM model parameters
predict_sem <- function(fit, data, x_names, y_names) {
    implied_values <- lav_model_implied(fit@Model)
    Sigma <- implied_values$cov[[1]]

    # Formula variables
    Sigma_XX <- Sigma[x_names, x_names]
    Sigma_YX <- Sigma[y_names, x_names]

    Sigma_XX_inv <- solve(Sigma_XX)  
    S_YX_S_XX_inv <- Sigma_YX %*% Sigma_XX_inv  # Σ_YX * Σ_XX^(-1)

    # Model-implied means
    mu_X <- implied_values$mean[[1]][x_names]  
    mu_Y <- implied_values$mean[[1]][y_names]

    centered_X <- as.matrix(sweep(data[x_names], 2, mu_X))

    # Transpose the centered_X before multiplication
    predictions <- t(mu_Y + S_YX_S_XX_inv %*% t(centered_X))
    
    return(predictions)
}

predicty.lavaan = function(object, newdata, xnames, ynames){
  # predict function for computing predicted values for set of response variables
  # given a set of predictor variables
  # based on the joint distribution estimated by a SEM model
  # INPUT
  # object: lavaan output object obtained from sem()
  # newdata: data frame with new values for the predictors
  # xnames: variables designated as predictors
  # ynames: variables designated as response variables
  
  #
  Sxx = fitted(object)$cov[xnames , xnames]
  Sxy = fitted(object)$cov[xnames , ynames]
  mx = fitted(object)$mean[xnames]
  my = fitted(object)$mean[ynames]
  
  #
  Xtest = as.matrix(newdata[, xnames])
  Xtest = scale(Xtest, center = mx, scale = FALSE)
  yhat = matrix(my, nrow = nrow(Xtest), ncol = length(ynames), byrow = TRUE) + Xtest %*% solve(Sxx) %*% Sxy
  
  # return
  return(yhat)
}