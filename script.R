library(lavaan)
library(mvtnorm)

# Define the SEM model
model <- '
  # Define latent variables
  latent1 =~ x1 + x2 + x3
  latent2 =~ x4 + x5 + x6

  # Define covariance between latent variables
  latent1 ~~ latent2

  # Define means of latent variables
  latent1 ~ 1
  latent2 ~ 1
'

# Create a sample dataset (not used for fitting, just for variable names)
set.seed(123)
data <- data.frame(matrix(rnorm(300), nrow=50, ncol=6))
colnames(data) <- c("x1", "x2", "x3", "x4", "x5", "x6")

# Create a lavaan model object without fitting the model
fit <- sem(model, data = data, do.fit = FALSE)

x_names = c(1, 2, 3)
y_names = c(4, 5, 6)

predict_sem <- function(implied_values, data, x_names, y_names) {
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

get_new_fit <- function(new_params, fit, data) {
  par_table <- parTable(fit)
  par_table$est[par_table$free > 0] <- new_params
  
  updated_fit <- sem(model = model, data = data, start = new_params, do.fit = FALSE)
  return (updated_fit)
}

mse_sem <- function(new_params, fit, data, x_names, y_names) {
  # Make a new model with new parameters
  updated_fit <- get_new_fit(new_params, fit, data)
  implied_values <- lav_model_implied(updated_fit@Model)
  
  new_predictions = predict_sem(implied_values, data, x_names, y_names)
  
  mse <- mean(rowMeans(((actual_Y - predictions)^2)))
  return(mse)
}

loss_values <- c()
likelihood_sem <- function(new_params, fit, data, x_names, y_names) {
  # Make a new model with new parameters
  updated_fit <- get_new_fit(new_params, fit, data)
  implied_values <- lav_model_implied(updated_fit@Model)
  
  Sigma <- implied_values$cov[[1]]
  
  # Formula variables
  Sigma_XX <- Sigma[x_names, x_names]
  Sigma_YX <- Sigma[y_names, x_names]
  Sigma_YY <- Sigma[y_names, y_names]
  
  if (det(Sigma_XX) <= 1e-10) {
    return(1e6)  # Large penalty instead of Inf
  }
  #Sigma_XX_inv <- solve(Sigma_XX)  
  Sigma_XX_inv <- solve(Sigma_XX + diag(1e-6, nrow(Sigma_XX)))  
  S_YX_S_XX_inv <- Sigma_YX %*% Sigma_XX_inv  # Σ_YX * Σ_XX^(-1)
  
  # Model-implied means
  mu_X <- implied_values$mean[[1]][x_names]  
  mu_Y <- implied_values$mean[[1]][y_names]
  
  pred_cov <- Sigma_YY - S_YX_S_XX_inv %*% Sigma_YX
  pred_cov <- (pred_cov + t(pred_cov)) / 2 # Fix to make it symmetric
  
  if (any(eigen(pred_cov)$values <= 0)) {
    return(1e6)  # Large penalty instead of Inf
  }
  # print(pred_cov)
  
  minus2_log_likelihood <- 0
  for (i in 1:nrow(data)) {
    row_pred_mean <- mu_Y + S_YX_S_XX_inv %*% (t(data[i, x_names, drop = FALSE]) - mu_X)
    
    row_ll <- dmvnorm(data[i, y_names], mean = row_pred_mean, sigma = pred_cov, log = TRUE)
    minus2_log_likelihood <- minus2_log_likelihood + row_ll
  }
  
  loss_value <- -2 * minus2_log_likelihood
  loss_values <<- c(loss_values, loss_value)
  return(-2 * minus2_log_likelihood)
}

free_params <- parTable(fit)$est[parTable(fit)$free > 0]
opt_result <- optim(
  par = free_params, 
  fn = likelihood_sem, 
  fit = fit, 
  data = data, 
  x_names = x_names, 
  y_names = y_names,
  method = "Nelder-Mead"
)
opt_result

