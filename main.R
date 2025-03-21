library(lavaan)
source("src/sem_utils.R")
source("src/objective_functions.R")
source("src/data_generation.R")
set.seed(123)

# Define the SEM model
model <- '
    # Define latent variables
    # latent1 =~ 1*x1 + x2 + x3
    # define correlations

    # Define means of latent variables
    x1 ~ 1
    x2 ~ 1
    x3 ~ 1
    x4 ~ 1

    x4 ~ x1 + x2 + x3
'

# Define variable indices for predictors and outcome
x_names <- c(1, 2, 3)
y_names <- c(4)

# Generate synthetic data
data <- generate_data(n = 100)

# Create a lavaan model object without fitting
zero_fit <- sem(model, data = data, do.fit = FALSE)

# Extract initial free parameters
free_params <- parTable(zero_fit)$est[parTable(zero_fit)$free > 0]

# Optimize the model using maximum likelihood
opt_result <- optim(
    par = free_params,
    fn = likelihood_sem,
    fit = zero_fit,
    data = data,
    x_names = x_names,
    y_names = y_names,
    method = "BFGS",
    #control = list(maxit = 5000, reltol = 1e-10)
)

# Models
optimized_fit <- get_new_fit(opt_result$par, zero_fit)
lavaan_fit <- sem(model, data = data)
lm_model <- lm(x4 ~ x1 + x2 + x3, data = data)

# Compare MSE on test data
test_data = generate_data(n = 100)
mse_lm <- calculate_mse(
    predict(lm_model, test_data),
    test_data$x4
)
mse_lavaan <- mse_sem(
    lavaan_fit,
    test_data,
    x_names,
    y_names
)
mse_optimized <- mse_sem(
    optimized_fit,
    test_data,
    x_names,
    y_names
)