library(ggplot2)
library(matrixcalc)
library(MASS)
library(mvtnorm)
library(LearnBayes)
library(nlme)

continuous_updating_by_n <- function(update_n_times, init_sample_size, new_sample_size){
  var.cov.mat <- matrix(c(1, 0.6, 0.6, 1), nrow = 2, byrow = TRUE)
  
  random_data <- mvrnorm(n = init_sample_size, mu = c(2, 3), Sigma = var.cov.mat)
  random_data <- as.data.frame(random_data)
  
  f_init <- lm(V2 ~ V1, data = random_data)
  sigma0 <- vcov(f_init)
  mu0 <- coef(f_init)
  
  errors_intercept <- c()
  errors_coeff <- c()
  
  for(i in 1:update_n_times){
    new_data_length <- new_sample_size
    new_random_data <- mvrnorm(n = new_data_length, mu = c(4, 5), Sigma = var.cov.mat)
    new_random_data <- as.data.frame(new_random_data)
    
    #Append the data and graph
    combined_random_data_x <- append(unlist(random_data["V1"], use.names = FALSE), values = unlist(new_random_data["V1"], use.names = FALSE))
    combined_random_data_y <- append(unlist(random_data["V2"], use.names = FALSE), values = unlist(new_random_data["V2"], use.names = FALSE))
    combined_random_data <- data.frame(combined_random_data_x, combined_random_data_y)
    
    f_new_data <- lm(V2 ~ V1, data = new_random_data)
    sigma_new_data <- vcov(f_new_data)
    mu_new_data <- coef(f_new_data)
    
    inverse_data_length <- 1
    #sum_of_new_data <- sum(new_random_data["V1"])
    sum_of_new_data <- c(sum(new_random_data["V1"]), sum(new_random_data["V2"]))
    sigma_n <- sigma0%*%matrix.inverse((sigma0 + inverse_data_length*sigma_new_data))%*%(inverse_data_length*sigma_new_data)
    
    mu_n <- sigma0%*%matrix.inverse(sigma0 + inverse_data_length*sigma_new_data)%*%(inverse_data_length*mu_new_data) + inverse_data_length*sigma_new_data%*%matrix.inverse(sigma0 + inverse_data_length*sigma_new_data)%*%mu0
    
    f_final <- lm(combined_random_data_y ~ combined_random_data_x, data = combined_random_data)
    sigma_actual <- vcov(f_final)
    mu_actual <- coef(f_final)
    
    errors_intercept <- c(errors_intercept, 100*abs(mu_actual[1]-mu_n[1])/mu_actual[1])
    errors_coeff <- c(errors_coeff, 100*abs(mu_actual[2]-mu_n[2])/mu_actual[2])
    
    mu0 <- mu_n
    sigma0 <- sigma_n
    
    random_data <- rbind(random_data, new_random_data)
    print(i)
  }
  data.frame(errors_intercept = errors_intercept, errors_coeff = errors_coeff, n = 1:update_n_times)
}
res11 <- continuous_updating_by_n(500, 1000, ceiling(sqrt(1000)))
ggplot(data = res11, aes(x = n, y = errors_intercept)) + 
  geom_point(aes(y = errors_intercept, col = 'red')) +
  stat_smooth(method = 'loess', color = 'blue') + 
  ggtitle("n=2000, initial_size=1000, batch_size=32") + 
  ylab("Error for Intercept (%)")
ggplot(data = res11, aes(x = n, y = errors_coeff)) + 
  geom_point(aes(y = errors_coeff, col = 'red')) +
  stat_smooth(method = 'loess', color = 'blue') + 
  ggtitle("n=2000, initial_size=1000, batch_size=32") + 
  ylab("Error for Coefficient (%)")