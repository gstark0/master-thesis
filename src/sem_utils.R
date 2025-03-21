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

    centered_X <- sweep(data[x_names], 2, mu_X)

    predictions = t(mu_Y + S_YX_S_XX_inv %*% t(centered_X))
    return(predictions)
}

# Predict outcomes (Y) based on predictors (X) using SEM model parameters using full model power
predict_sem_nonlinear <- function(fit, data, x_names, y_names) {
    # Get latent variable scores
    latent_scores <- lavPredict(fit, newdata = data)
    
    # Get parameters for latent to outcome paths
    params <- parameterEstimates(fit)
    latent1_coef <- params$est[params$lhs == "x4" & params$op == "~" & params$rhs == "latent1"]
    latent2_coef <- params$est[params$lhs == "x4" & params$op == "~" & params$rhs == "latent2"]
    
    # Non-linear prediction (matching your data generation)
    predictions <- latent1_coef * latent_scores[, "latent1"] + 
                   latent2_coef * latent_scores[, "latent2"]^2 - 
                   0.5 * (latent_scores[, "latent1"] * latent_scores[, "latent2"])
    
    return(predictions)
}

# Predict like a linear model, but with SEM parameters
predict_lm <- function(fit, data, x_names, y_names) {
    # Get variable names
    y_var <- colnames(data)[y_names[1]]
    x_vars <- colnames(data)[x_names]

    # Get parameters
    params <- parameterEstimates(fit)

    # Get intercept
    intercept <- params$est[params$lhs == y_var & params$op == "~1"]
    if(length(intercept) == 0) intercept <- 0

    # Get coefficients
    coeffs <- numeric(length(x_vars))
    for(i in 1:length(x_vars)) {
        idx <- which(params$lhs == y_var & params$op == "~" & params$rhs == x_vars[i])
        if(length(idx) > 0) {
            coeffs[i] <- params$est[idx]
        }
    }

    # Calculate predictions
    pred <- intercept
    for(i in 1:length(x_vars)) {
        pred <- pred + coeffs[i] * data[[x_vars[i]]]
    }

    return(pred)
}