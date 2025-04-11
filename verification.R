library(lavaan)
source("src/sem_utils.R")
source("src/objective_functions.R")
source("src/data_generation.R")
set.seed(123)

x_names <- c(1, 2)
y_names <- c(3)

set.seed(123)

# Define the SEM model
get_model_specification <- function(scenario = "ideal") {
    if (scenario == "adventage") {
        return('
            # Intercepts for all variables
            x1 ~ 1
            x2 ~ 1
            x3 ~ 1
            x4 ~ 1
            x5 ~ 1
            x6 ~ 1
            x7 ~ 1
            x8 ~ 1
            x9 ~ 1
            x10 ~ 1
            
            # Define latent variables
            latent1 =~ x1 + x2 + x3
            latent2 =~ x4 + x5 + x6
            latent3 =~ x7 + x8 + x9
            
            # Define structural paths to outcome
            x10 ~ latent1 + latent2 + latent3
        ')
    }
  else if (scenario %in% c("ideal", "skewed", "high_noise")) {
    return('
            x1 ~ 1
            x2 ~ 1
            x3 ~ 1
            latent1 =~ x1 + x2 + x3
    
            # Structural model with latent variable
            x4 ~ latent1
            
            # Allow for residual correlations to capture non-modeled complexity
            x1 ~~ x2
            x2 ~~ x3
            x1 ~~ x3
        ')

  } else if (scenario == "misspecified") {
    return('
            x1 ~ 1
            x2 ~ 1
            x3 ~ 1
            latent1 =~ x1 + x2 + x3
            x4 ~ latent1 + x1
    ')  # intentionally redundant

  } else if (scenario == "regression") {
    return('
        x1 ~ 1
        x2 ~ 1
        x3 ~ x1 + x2
    ')

  } else if (scenario == "latent_signal_only") {
    return('
        x1 ~ 1
        x2 ~ 1
        x3 ~ 1
        latent1 =~ x1 + x2 + x3
        x4 ~ latent1
    ')

  } else {
    stop("Unknown scenario")
  }
}


# Scenario testing function
test_scenario <- function(ntraining, ntesting, scenario='ideal') {
    #data <- generate_sem_data(n=ntraining, scenario)
    data <- generate_sem_data(ntraining, scenario = scenario)
    model = get_model_specification(scenario)
    zero_fit <- sem(model, data = data, do.fit = FALSE)
    lavaan_fit <- sem(model, data = data, meanstructure = FALSE)
    #free_params <- parTable(zero_fit)$est[parTable(zero_fit)$free > 0]
    free_params <- parTable(lavaan_fit)$est[parTable(lavaan_fit)$free > 0]

    opt_result <- optim(
        par = free_params,
        fn = mse_objective,
        fit = lavaan_fit,
        data = data,
        x_names = x_names,
        y_names = y_names,
        method = "BFGS",
        control=c(maxit = 2000, reltol = 1e-10),
        hessian = TRUE,
    )

    # ------------------ HESSIAN ------------------
    hessian_matrix <- opt_result$hessian
    eigen_results <- eigen(hessian_matrix)
    eigen_values <- eigen_results$values
    cat("Smallest eigenvalue:", min(eigen_values), "\n")
    cat("Largest eigenvalue:", max(eigen_values), "\n")
    # ----------------------------------------------

    test_data = generate_sem_data(ntesting, scenario = scenario)

    # Models
    optimized_fit <- get_new_fit(opt_result$par, zero_fit)
    lavaan_fit <- sem(model, data = data, meanstructure = FALSE)
    # This needs a redefinition for different scenarios !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    lm_model <- lm(x3 ~ x1 + x2, data = data)

    # Compare MSE on test data
    mse_lm <- calculate_mse(
        predict(lm_model, test_data),
        test_data$x3
    )

    actual_Y <- test_data[, y_names, drop = FALSE][, 1]
    # This needs a redefinition for different scenarios !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    mse_lavaan_paper = mean((actual_Y - predicty.lavaan(lavaan_fit, test_data, c("x1", "x2"), c("x3")))^2)
    mse_optimized_paper = mean((actual_Y - predicty.lavaan(optimized_fit, test_data, c("x1", "x2"), c("x3")))^2)

    return(
        list(
            mse_lm = mse_lm,
            mse_lavaan_paper = mse_lavaan_paper,
            mse_optimized_paper = mse_optimized_paper,
            min_eigenvalue = min(eigen_values),
            max_eigenvalue = max(eigen_values),
            lavaan_fit = lavaan_fit,
            optimized_fit = optimized_fit
        )
    )
}

# Test different scenarios
scenarios <- c('regression')
results <- lapply(scenarios, function(scenario) {
    test_scenario(100000, 1000, scenario)
})
print(results)
# print paramaters
#print(results[[1]]$lavaan_fit@ParTable)
#print(results[[1]]$optimized_fit@ParTable)