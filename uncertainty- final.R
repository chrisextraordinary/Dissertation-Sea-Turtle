# Function for alpha and beta
calculate_beta_parameters <- function(mu, sigma) {
  var <- sigma^2
  common_factor <- mu * (1 - mu) / var - 1
  alpha <- mu * common_factor
  beta <- (1 - mu) * common_factor
  return(c(alpha = alpha, beta = beta))
}

calculate_lambda_distribution <- function(mu, sigma, n_iterations = 10000, stage, stage_duration) {
  # Get alpha and beta
  params <- calculate_beta_parameters(mu, sigma)
  
  # Initialize storage vectors
  new_vals <- numeric(n_iterations)
  lambda_values <- numeric(n_iterations)
  
  for (i in 1:n_iterations) {
    
    new_L <- L  # Ensure L is your original matrix
    
    # Generate a new value from the beta distribution
    new_survival_rate <- rbeta(1, params['alpha'], params['beta'])
    
    # Assuming calc_P function exists, replace with your specific calculation
    new_P <- calc_P(new_survival_rate, stage_duration)
    new_L[stage, stage] <- new_P
    
    # Compute the new lambda
    eigen_vals <- eigen(new_L)$values
    new_lambda <- max(Re(eigen_vals))
    
    # Store the sampled value and lambda
    new_vals[i] <- new_P
    lambda_values[i] <- new_lambda
  }
  
  return(list(
    mean = mean(lambda_values),
    min = min(lambda_values),
    max = max(lambda_values)
  ))
}

calculate_lambda_distribution_G <- function(mu, sigma, n_iterations = 10000, stage, stage_duration) {
  # Get alpha and beta
  params <- calculate_beta_parameters(mu, sigma)
  
  # Initialize storage vectors
  new_vals_G <- numeric(n_iterations)
  lambda_values_G <- numeric(n_iterations)
  
  for (i in 1:n_iterations) {
    
    new_L <- L  # Ensure L is your original matrix
    
    # Generate a new value from the beta distribution
    new_survival_rate <- rbeta(1, params['alpha'], params['beta'])
    
    # Assuming calc_G function exists, replace with your specific calculation
    new_G <- calc_G(new_survival_rate, stage_duration)
    new_L[stage, stage-1] <- new_G
    
    # Compute the new lambda
    eigen_vals <- eigen(new_L)$values
    new_lambda <- max(Re(eigen_vals))
    
    # Store the sampled value and lambda
    new_vals_G[i] <- new_G
    lambda_values_G[i] <- new_lambda
  }
  
  return(list(
    mean = mean(lambda_values_G),
    min = min(lambda_values_G),
    max = max(lambda_values_G)
  ))
}

# Uncertainty of lambda in changes of adult survival rate 
mu_adult_sur <- mean(c(0.79, 0.83, 0.8091))
sigma_adult_sur <- sd(c(0.79, 0.83, 0.8091))
results_adult_sur_P <- calculate_lambda_distribution(mu_adult_sur, sigma_adult_sur, stage = 4, stage_duration = 19)
# Uncertainty for G4 (it is an exception cuz it needs to time 1/4.4 the remigration interval)
params_adult_sur <- calculate_beta_parameters(mu_adult_sur, sigma_adult_sur)
new_vals_adult_sur_G <- numeric(10000)
lambda_adult_sur_G <- numeric(10000)

for (i in 1:10000) {
  new_L_G4 <- L  # Ensure L is your original matrix
  
  # Generate a new value from the beta distribution
  new_adult_sur <- rbeta(1,  params_adult_sur[[1]],  params_adult_sur[[2]])
  
  # Update G4 with new survival
  new_G4 <- calc_G(new_adult_sur, 19)*1/4.4
  new_L_G4[4, 3] <- new_G4
  
  # Compute the new lambda
  eigen_vals_G4 <- eigen(new_L_G4)$values
  new_lambda_adult_sur_G4 <- max(Re(eigen_vals_G4))
  
  # Store the sampled value and lambda
  new_vals_adult_sur_G[i] <- new_G4
  lambda_adult_sur_G[i] <- new_lambda_adult_sur_G4
}

mean(lambda_adult_sur_G)
min(lambda_adult_sur_G)
max(lambda_adult_sur_G)

results_adult_sur_G <- c(mean(lambda_adult_sur_G), min(lambda_adult_sur_G), max(lambda_adult_sur_G))

juv_ner_sur <- c(0.674, 0.893, 0.6758, 0.7425)
# Compute mean and standard deviation for juv_ner_sur
mu_juv_ner_sur <- mean(juv_ner_sur)
sigma_juv_ner_sur <- sd(juv_ner_sur)
results_juv_ner_P <- calculate_lambda_distribution(mu_juv_ner_sur, sigma_juv_ner_sur, stage = 3, stage_duration = 20)
results_juv_ner_G <- calculate_lambda_distribution_G(mu_juv_ner_sur, sigma_juv_ner_sur, stage = 3, stage_duration = 20)

# for juvenile oceanic, data mainly from bjorndal 2003
juv_ocean_sur <- c(0.720, 0.911, .643, .656, .894, .608, .725)
# Compute mean and standard deviation for juv_ocean_sur
mu_juv_ocean_sur <- mean(juv_ocean_sur)
sigma_juv_ocean_sur <- sd(juv_ocean_sur)
results_juv_ocean_P <- calculate_lambda_distribution(mu_juv_ocean_sur, sigma_juv_ocean_sur, stage = 2, stage_duration = 8.2)
results_juv_ocean_G <- calculate_lambda_distribution_G(mu_juv_ocean_sur, sigma_juv_ocean_sur, stage = 2, stage_duration = 8.2)

df <- data.frame(Stage = c("Juvenile Oceanic", "Juvenile Neritic", "Adult Female", "Juvenile Oceanic", "Juvenile Neritic", "Adult Female"),
                 Parameter = c("P", "P", "P", "G", "G", "G"),
                 Max_Lambda = c(results_juv_ocean_P$max, results_juv_ner_P$max, results_adult_sur_P$max, 
                                results_juv_ocean_G$max, results_juv_ner_G$max, results_adult_sur_G[3]),
                 Mean_Lambda = c(results_juv_ocean_P$mean, results_juv_ner_P$mean, results_adult_sur_P$mean,
                                 results_juv_ocean_G$mean, results_juv_ner_G$mean, results_adult_sur_G[1]),
                 Min_Lambda = c(results_juv_ocean_P$min, results_juv_ner_P$min, results_adult_sur_P$min,
                                results_juv_ocean_G$min, results_juv_ner_G$min, results_adult_sur_G[2]))

ggplot(df, aes(x = as.factor(Stage), y = Mean_Lambda, colour = Parameter)) +
  geom_linerange(aes(ymin = Min_Lambda, ymax = Max_Lambda), position = position_dodge(width = 0.5), linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.5), size = 1) +
  labs(title = "Uncertainty Analysis",
       x = "Stages",
       y = "Range of Lambda") +
  theme_bw() +
  scale_y_continuous(limits = c(0.75, 1.25), expand = c(0, 0)) +
  theme(axis.text.x = element_text(),
        plot.title = element_text(hjust = 0.5))  # Corrected the specification for centering the title

