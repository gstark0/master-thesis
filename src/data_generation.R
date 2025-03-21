# data_generation.R

generate_data <- function(n) {
    # Latent variables
    latent1 <- rnorm(n)
    latent2 <- 0.5 * latent1 + rnorm(n)  # latent2 is correlated with latent1

    # Indicators
    x1 <- 0.8 * latent1 + rnorm(n, sd = 0.5)
    x2 <- 0.7 * latent1 + rnorm(n, sd = 0.5)
    x3 <- 0.6 * latent1 + rnorm(n, sd = 0.5)
    x4 <- 0.5 * x1 + 0.3 * x2 + 0.2 * x3 + rnorm(n, sd = 0.4)

    # Create the dataset
    data <- data.frame(x1, x2, x3, x4)

    data <- as.data.frame(data)

    return(data)
}

generate_data_nonlinear <- function(n) {
    # Latent variables
    latent1 <- rnorm(n)
    latent2 <- 0.5 * latent1 + rnorm(n)

    # Indicators
    x1 <- 0.8 * latent1 + rnorm(n, sd = 0.3)
    x2 <- 0.7 * latent1 + rnorm(n, sd = 0.3)
    x3 <- 0.6 * latent1 + rnorm(n, sd = 0.3)

    # non-linear relationships
    x4 <- 0.2 * latent1 + 0.6 * latent2^2 - 0.5 * (latent1 * latent2) + 0.3 * sin(latent1) + rnorm(n, sd = 0.3)

    # Create the dataset
    data <- data.frame(x1, x2, x3, x4)
    return(data)
}

# A data where SEM should show its potential with latent variables :P
generate_sem_advantage_data <- function(n) {
  # Create latent variables
  latent1 <- rnorm(n)
  latent2 <- 0.5 * latent1 + rnorm(n, sd = 0.7)
  
  # Introduce heavy measurement error to observed variables
  # This makes the latent structure more important
  x1 <- 0.7 * latent1 + rnorm(n, sd = 1.2)  # High measurement error
  x2 <- 0.6 * latent1 + rnorm(n, sd = 1.3)
  x3 <- 0.5 * latent1 + rnorm(n, sd = 1.4)
  
  # Make outcome directly dependent on latent variables only
  # This creates an indirect relationship that linear regression won't capture well
  x4 <- 0.8 * latent1 - 0.4 * latent2 + 0.3 * latent1 * latent2 + rnorm(n, sd = 0.5)
  
  data <- data.frame(x1, x2, x3, x4)
  return(data)
}