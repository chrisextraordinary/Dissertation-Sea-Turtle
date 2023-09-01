set.seed(123)
# Load libraries
library(ggplot2)
library(popbio)
library(dplyr)
library(tibble)
library(tidyr)

# Get the initial nesting female and total female population
# randomly drew 1,000 values of c from a normal distribution with mean = 1.7 and SD = 0.4
# Given parameters
ni <- 1345  # sum nest count for year i (2022)
n_2021 <- 1047
mean_c <- 1.7 # mean of clutch frequency
sd_c <- 0.4 # sd of clutch frequency

mean_r <- 4.4 # mean value of remigration interval
sd_r <- 2.4

# Drawing 1000 values of c from the normal distribution
set.seed(123)  # For reproducibility
c_samples <- rnorm(10000, mean_c, sd_c)

# make sure r is always > 0
r_samples <- rep(NA, 10000)
for (i in 1:10000) {
  sample <- rnorm(1, mean_r, sd_r)
  while (sample <= 0) {
    sample <- rnorm(1, mean_r, sd_r)
  }
  r_samples[i] <- sample
}

# Calculating ANFy for each sample
ANF_samples <- ni / c_samples # thats the nesting female number

# Calculating the mean and 95% confidence interval of ANFy
mean_ANF <- mean(ANF_samples)
conf_limits_ANF <- quantile(ANF_samples, c(0.025, 0.975))

print(paste("Mean of ANF:", mean_ANF))
print(paste("95% CI of ANF: [", conf_limits_ANF[1], ",", conf_limits_ANF[2], "]"))

# Calculating Tfy for each ANFy sample
Tf_samples <- ANF_samples * r_samples # thats the total female living in this area

# Calculating the mean and 95% confidence interval of Tfy
mean_Tf <- mean(Tf_samples)
conf_limits_Tf <- quantile(Tf_samples, c(0.025, 0.975))

print(paste("Mean of TF:", mean_Tf))
print(paste("95% CI of TF: [", conf_limits_Tf[1], ",", conf_limits_Tf[2], "]"))


ecs <- 101 # estimated clutch size
es <- .484 # emergence success
ini_post <- floor(ecs*ni*es)
ini_juv_oce <- floor(ini_post * 0.484) 
ini_juv_ner <- floor(ini_juv_oce * 0.212)
ini_adult <- floor(mean_Tf - mean_ANF)
ini_nesting <- floor(mean_ANF)

# The population matirx
# Load necessary libraries
library(tibble)

# Create the dataframe
loggerhead_data <- tibble(
  Stage = c("Post-Hatchlings", "Juvenile Oceanic", "Juvenile Neritic", "Adult", "Nesting Females"),
  Stage_Duration = c(NA, 8.2, 20, 19, NA),
  Survival_Rate = c(NA, 0.725, 0.893, 0.812, 0) 
)

loggerhead_data <- as.data.frame(loggerhead_data)

# Define the functions for the probability of surviving and remaining in the same stage (P)
# and the probability of surviving and growing into the next stage (G).
calc_P <- function(SR, SD) { # where SR is the survival rate and SD is the stage duration
  if (is.na(SR) | is.na(SD)) {
    return(NA)
  } else {
    return((1 - SR^(SD - 1)) / (1 - SR^SD) * SR)
  }
}

calc_G <- function(SR, SD) {
  if (is.na(SR) | is.na(SD)) {
    return(NA)
  } else {
    return((SR^SD) * (1 - SR) / (1 - SR^SD))
  }
}

# Apply the functions to the dataframe
loggerhead_data$P <- mapply(calc_P, loggerhead_data$Survival_Rate, loggerhead_data$Stage_Duration)
loggerhead_data$G <- mapply(calc_G, loggerhead_data$Survival_Rate, loggerhead_data$Stage_Duration)

# View the updated dataframe
print(loggerhead_data)

# Define stages and demographic parameters
stages <- c("Post Hatchling", "Juvenile Oceanic", "Juvenile Neritic", "Adult Female", "Nesting Female")

# survival probability of each stage, updated, might need to refine
P <- c(0, loggerhead_data$P[[2]], loggerhead_data$P[[3]], loggerhead_data$P[[4]], 0)  #0.808*(1-1/breeding)

# Growth probabilitiy of each stage, updated, mightneed to refine
G <- c(0.484, loggerhead_data$G[[2]], loggerhead_data$G[[3]], loggerhead_data$G[[4]])

# Fecundity = clutch frequency * clutch size * hatchling success * sex ratio (females only)
f <- round(1.7*101*0.873*0.5*0.516, 0)

# female surviving - breeding 
# Define the function
L <- matrix(c(0, 0, 0, f, 0, # P[1] is 0
              G[1], P[2], 0, 0, 0,
              0, G[2], P[3], 0, 0,
              0, 0, G[3], P[4]*(1-1/4.4), 0.73, #0.73 is the survival rate of nesting females, from 'Survival and remigration probabilities for loggerhead turtles (Caretta caretta) nesting in the eastern Gulf of Mexico'
              0, 0, 0, G[4]*1/4.4, 0), # P5 is 0
              nrow = 5, byrow = TRUE,
              dimnames = list(c("Post Hatchling", "Juvenile Oceanic", "Juvenile Neritic", 
                                "Adult", "Nesting Female"),
                              c("Post Hatchling", "Juvenile Oceanic", "Juvenile Neritic", 
                                "Adult", "Nesting Female")))

# Calculate eigenvalues
eigen_vals <- eigen(L)$values

# Find the dominant eigenvalue (lambda)
lambda_L <- max(Re(eigen_vals))

# Function to calculate stable stage distribution and reproductive values
df_SSD_RRV <- function(L) {
  
  # Stable Stage Distribution
  SSD <- round(stable.stage(L), 3)
  
  # Reproductive Value
  r <- round(reproductive.value(L), 3)
  
  # Create dataframe for the reproductive values
  df <- data.frame(Reproductive_Value = r, Stable_Stage_Distribution = SSD)
  
  # Convert row names to a column
  df <- df %>% rownames_to_column(var = "Stage")
  
  # Convert reproductuve value to relative reproductive value (RRV)
  df$RRV <- df$Reproductive_Value / df$Reproductive_Value[which.max(df$Reproductive_Value)]
  df$RRV <- round(df$RRV, 2)
  # Calculate the dominant eigenvalue of L, which is the population growth rate (lambda)
  lambda <- max(Re(eigen(L)$values))
  
  # Return a list containing the dataframe and the growth rate
  return(list(df = df, lambda = lambda))
}

df_SSD_RRV(L)

# Sensitivity Analysis
s <- sensitivity(L, zero = TRUE)

# Create vectors of Ps, Gs, and Fs sensitivities in each stage.
Ps <- rep(0, 5)
Gs <- rep(0, 5)
Fs <- rep(0, 5)
for (j in 1:5) {
  Ps[j] <- s[j,j]
  if (j != 5) {
    Gs[j] <- s[j+1,j]
  } else(Gs[j] <- 0)
  Fs[j] <- s[1,j]
}

sensitivity_data <- data.frame(Stage = stages, Fs = Fs, Ps = Ps, Gs = Gs)

dodge <- position_dodge(width=0.3)

sensitivity_data_long <- data.frame(
  Stage = rep(stages, each = 3),
  Factor = rep(c("Fs", "Ps", "Gs"), times = 5), 
  Sensitivity = c(Fs, Ps, Gs) 
)

# Convert sensitivity data from wide format to long format
sensitivity_long <- sensitivity_data %>%
  pivot_longer(cols = c("Fs", "Ps", "Gs"), names_to = "Parameter", values_to = "Sensitivity")

library(stringr)
ggplot(sensitivity_long, aes(x = Stage, y = Sensitivity, fill = Parameter)) +
  geom_bar(stat = 'identity', width = 0.5, position = "dodge") +
  labs(title = "Sensitivity of Different Parameters by Life Stages",
       x = "Life Stage",
       y = "Sensitivity",
       color = "Parameter") +
  scale_color_manual(values = c("Fs" = "green", "Ps" = "blue", "Gs" = "red")) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1,
                                   margin = margin(t = 5, r = 0, b = 5, l = 0),
                                   family = "serif", size = 10,
                                   lineheight = 0.8), plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme_bw()

# The elasticity of lambda to changes in F, G, and P
elas <- elasticity(L)

Pi <- rep(0, 5)
Gi <- rep(0, 5)
Fi <- rep(0, 5)

for (j in 1:5) {
  Pi[j] <- elas[j,j]
  if (j != 5) {
    Gi[j] <- elas[j+1,j]
  } else (Gi[j] <- 0)
  Fi[j] <- elas[1,j]
}

# Create a line plot
plot(Pi, type = "l", col = "blue",
     xlab = "Life Stage Number", ylab = "Elasticity of lambda", ylim = c(0, 0.7),
     main = "Elasticity of lambda vs. Life Stage Number")
lines(Gi, type = "l", col = "red")
lines(Fi, type = "l", col = "green")
legend("topright", legend = c("Pi", "Gi", "Fi"), col = c("blue", "red", "green"), lty = 1)
# The elasticities of lambda in Gi and Fi are rather low. 

# increase in population survival rate by percentage
new_survival <- c(1.01, 1.05, 1.10)

new_lam <- vector("list", length(new_survival))  # Initialize a list to store eigenvalues

for (stage in c(2, 3, 4)) {
  for (x in 1:length(new_survival)) {
    temp_L <- L  # Create a temporary copy of L for each iteration
    
    # Update the survival probability for the specified stage
    temp_L[stage, stage] <- L[stage, stage] * new_survival[x]  # increase survival probabilities in the specified stage by x%
    
    new_eigen_vals <- Re(eigen(temp_L)$values)
    new_lam[[x]] <- c(new_lam[[x]], as.numeric(new_eigen_vals[1]))  # Store the eigenvalue
  }
}

output <- matrix(NA, nrow = 3, ncol = 3)
for (w in c(1:3)) {
  output[w, ] <- as.vector(new_lam[[w]])[which(as.vector(new_lam[[w]]) >0)]
}

lower <- 0.80

# proportion changes on lambda - elasticity analysis 
par(mar=c(5, 5, 4, 4))
barplot(output-lower, beside=TRUE, ylim=c(lower, 1.12), xlab = "Life Stages", ylab="Population Growth Rate",
        col=c("#007765","#6AA341","gray90"), offset= lower)
abline(h=lower, lwd=2)
abline(h=1.0, lwd=2, col="blue")
abline(h= lambda_L, lwd=2, col="orange")
text(x=c(4.5, 8.3, 12.2), y=0.78,
     labels=c("Juvenile Oceanic","Juvenile Neritic","Adult Female"),
     pos=2, xpd=TRUE, cex=0.9)
title(main = "Elasticity of Lambda to Changes in Survival Probibility")
legend("topright", inset=c(0.01, 0.01), c("+1%","+5%","+10%"),
       fill=c("#007765","#6AA341","gray90"), cex=0.7)

stored_matrices <- list()

for (x in 1:length(new_survival)) {
    temp <- L  # Create a temporary copy of L for each iteration
    
    # Update the survival probability for the juvenile neritic
    temp[3, 3] <- temp[3, 3] * new_survival[x]  # increase survival probabilities in the specified stage by x%
    stored_matrices[[x]] <- temp
}

# Function to simulate and plot growth
growth_plot <- function(L, time_period = 100, plot = TRUE) {
  # Define stages
  stages <- c("Post Hatchling", "Juvenile Oceanic", "Juvenile Neritic", "Adult Female", "Nesting Female")
  
  # Initialize matrix to hold population data
  N <- matrix(NA, nrow = 5, ncol = time_period)
  N[, 1] <- c(ini_post, ini_juv_oce, ini_juv_ner, ini_adult, ini_nesting)
  
  # Simulate population over time
  for (i in 2:time_period) {
    N[, i] <- L %*% N[, i-1]
  }
  
  # Convert matrix to data frame for plotting
  proj_df <- data.frame(N = c(N),
                        Time = rep(1:time_period, each = 5),
                        Stage = factor(rep(stages, time_period)))
  
  # Plot data 
    ggplot(proj_df, aes(x = Time, y = N, group = Stage, color = Stage)) +
      geom_line() +
      geom_point() +
      theme_bw() +
      facet_wrap(.~Stage, nrow = 2, scales = "free") +
      scale_fill_brewer(palette = "Set1") +
      scale_y_continuous(limits = c(0, NA))
}

growth_plot(L)

# population growth plots for 3 new matrixs
for (i in 1:3) {
  print(growth_plot(stored_matrices[[i]]))
}

# Initialize a data frame to hold the results
results <- list()

# loop over 3 new matrixs and get their df_SSD-RRV results
for(i in 1:3) {
  results[[i]] <- df_SSD_RRV(stored_matrices[[i]])
}

# Print the results
print(results)

growth_plot(stored_matrices[[3]], time_period = 200)
df_SSD_RRV(stored_matrices[[3]])
