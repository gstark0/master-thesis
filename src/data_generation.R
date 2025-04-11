# data_generation.R

generate_sem_data <- function(n = 300, scenario = "ideal") {
  if (scenario == "ideal") {
    latent1 <- rnorm(n)
    x1 <- 0.8 * latent1 + rnorm(n, sd = 0.5)
    x2 <- 0.7 * latent1 + rnorm(n, sd = 0.5)
    x3 <- 0.6 * latent1 + rnorm(n, sd = 0.5)
    x4 <- 0.9 * latent1 + rnorm(n, sd = 0.5)

  } else if (scenario == "misspecified") {
    latent1 <- rnorm(n)
    x1 <- 0.8 * latent1 + rnorm(n, sd = 0.5)
    x2 <- 0.7 * latent1 + rnorm(n, sd = 0.5)
    x3 <- 0.6 * latent1 + rnorm(n, sd = 0.5)
    x4 <- 0.5 * latent1 + 0.5 * x1 + rnorm(n, sd = 0.5)

  } else if (scenario == "skewed") {
    latent1 <- scale(exp(rnorm(n)))
    x1 <- 0.8 * latent1 + rnorm(n, sd = 0.5)
    x2 <- 0.7 * latent1 + rnorm(n, sd = 0.5)
    x3 <- 0.6 * latent1 + rnorm(n, sd = 0.5)
    x4 <- 0.9 * latent1 + rnorm(n, sd = 0.5)

  } else if (scenario == "regression") {
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    x3 <- x1 + x2 + rnorm(n, sd = 1)
    return (data.frame(x1, x2, x3))

  } else if (scenario == "high_noise") {
    latent1 <- rnorm(n)
    x1 <- 0.8 * latent1 + rnorm(n, sd = 1.5)
    x2 <- 0.7 * latent1 + rnorm(n, sd = 1.5)
    x3 <- 0.6 * latent1 + rnorm(n, sd = 1.5)
    x4 <- 0.9 * latent1 + rnorm(n, sd = 1.5)

  } else if (scenario == "latent_signal_only") {
    latent1 <- rnorm(n)

    # Very noisy predictors
    x1 <- 0.8 * latent1 + rnorm(n, sd = 1.5)
    x2 <- 0.7 * latent1 + rnorm(n, sd = 1.5)
    x3 <- 0.6 * latent1 + rnorm(n, sd = 1.5)

    x4 <- 0.9 * latent1 + rnorm(n, sd = 0.3)
    } else if (scenario == "adventage") {

    latent_cor <- matrix(c(
      1.0, 0.5, 0.3,
      0.5, 1.0, 0.4,
      0.3, 0.4, 1.0
    ), nrow = 3, ncol = 3)
    
    # Generate latent scores from multivariate normal distribution
    latent_scores <- MASS::mvrnorm(n, mu = c(0, 0, 0), Sigma = latent_cor)
    latent1 <- latent_scores[, 1]
    latent2 <- latent_scores[, 2]
    latent3 <- latent_scores[, 3]
    
    # Generate indicators for latent1
    x1 <- 0.8 * latent1 + rnorm(n, sd = 0.5)
    x2 <- 0.7 * latent1 + rnorm(n, sd = 0.6)
    x3 <- 0.6 * latent1 + rnorm(n, sd = 0.7)
    
    # Generate indicators for latent2
    x4 <- 0.8 * latent2 + rnorm(n, sd = 0.5)
    x5 <- 0.7 * latent2 + rnorm(n, sd = 0.6)
    x6 <- 0.6 * latent2 + rnorm(n, sd = 0.7)
    
    # Generate indicators for latent3
    x7 <- 0.8 * latent3 + rnorm(n, sd = 0.5)
    x8 <- 0.7 * latent3 + rnorm(n, sd = 0.6)
    x9 <- 0.6 * latent3 + rnorm(n, sd = 0.7)
    
    # Create outcome variable that depends on all latent variables
    x10 <- 0.7 * latent1 + 0.3 * latent2 - 0.2 * latent3 + rnorm(n, sd = 0.4)
    
    # Return data frame with all variables
    return(data.frame(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10))

  } else {
    stop("Unknown scenario.")
  }
  return(data.frame(x1, x2, x3, x4))
}