# objective_functions.R

library(mvtnorm)

# Calculate mean squared error for SEM predictions
mse_sem <- function(fit, data, x_names, y_names) {
    # Make predictions
    predictions <- predict_sem_nonlinear(fit, data, x_names, y_names)

    # Calculate MSE
    actual_Y <- data[, y_names, drop = FALSE][, 1]
    mse <- mean((actual_Y - predictions)^2)

    return(mse)
}

# Calculate negative log-likelihood for SEM model
likelihood_sem <- function(new_params, fit, data, x_names, y_names) {
    # Make a new model with new parameters
    updated_fit <- get_new_fit(new_params, fit)
    implied_values <- lav_model_implied(updated_fit@Model)

    Sigma <- implied_values$cov[[1]]

    # Formula variables
    Sigma_XX <- Sigma[x_names, x_names]
    Sigma_YX <- Sigma[y_names, x_names]
    Sigma_YY <- Sigma[y_names, y_names]

    if (any(is.na(Sigma)) || any(is.infinite(Sigma))) {
    return(1e6)  # Large penalty instead of Inf
    }
    #Sigma_XX_inv <- solve(Sigma_XX)  
    #Sigma_XX_inv <- solve(Sigma_XX + diag(1e-6, nrow(Sigma_XX)))
    Sigma_XX_inv <- solve(Sigma_XX + diag(1e-6, nrow(Sigma_XX)))
    S_YX_S_XX_inv <- Sigma_YX %*% Sigma_XX_inv  # Σ_YX * Σ_XX^(-1)

    # Model-implied means
    mu_X <- implied_values$mean[[1]][x_names]  
    mu_Y <- implied_values$mean[[1]][y_names]

    #print(implied_values$mean)

    pred_cov <- Sigma_YY - S_YX_S_XX_inv %*% Sigma_YX
    pred_cov <- (pred_cov + t(pred_cov)) / 2 # Fix to make it symmetric
    pred_cov <- pred_cov + diag(1e-6, nrow(pred_cov))

    if (any(eigen(pred_cov, symmetric = TRUE, only.values = TRUE)$values <= 1e-6)) {
    return(1e6)
    }
    # print(pred_cov)

    minus2_log_likelihood <- 0
    for (i in 1:nrow(data)) {
    row_pred_mean <- mu_Y + S_YX_S_XX_inv %*% (t(data[i, x_names, drop = FALSE]) - mu_X)

    row_ll <- dmvnorm(data[i, y_names], mean = row_pred_mean, sigma = pred_cov, log = TRUE)
    minus2_log_likelihood <- minus2_log_likelihood + row_ll
    }
    return(-2 * minus2_log_likelihood)
}

mse_objective <- function(new_params, fit, data, x_names, y_names) {
    updated_fit <- get_new_fit(new_params, fit)
    predictions <- predict_lm(updated_fit, data, x_names, y_names)
    return(mean((data$x4 - predictions)^2))
}

prediction_likelihood_sem <- function(new_params, fit, data, x_names, y_names) {
    # Make a new model with new parameters
    updated_fit <- get_new_fit(new_params, fit)
    
    # Get predictions using SEM approach
    predictions <- predict_sem(updated_fit, data, x_names, y_names)
    
    # Calculate prediction-focused negative log-likelihood
    pred_error <- data[, y_names, drop = FALSE][, 1] - predictions
    n <- nrow(data)
    
    # Minimize negative log-likelihood of prediction errors
    sigma2 <- mean(pred_error^2)  # MLE of error variance
    ll <- -n/2*log(2*pi*sigma2) - sum(pred_error^2)/(2*sigma2)
    
    return(-ll)  # Return negative log-likelihood
}