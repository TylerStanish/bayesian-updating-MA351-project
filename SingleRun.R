library(ggplot2)
library(matrixcalc)
library(MASS)
library(mvtnorm)
library(LearnBayes)
library(nlme)

run <- function(init_sample_size, new_sample_size){
  # adjust the 0.6 to whatever you want the covariance to be between the variables
  var.cov.mat <- matrix(c(1, 0.6, 0.6, 1), nrow = 2, byrow = TRUE)
  # variance-covariance matrices must be at least positive semi-definite
  is.positive.semi.definite(var.cov.mat)
  # is.positive.definite(var.cov.mat)
  
  # Thanks to this: https://stat.ethz.ch/pipermail/r-help/2007-April/128925.html
  random_data <- mvrnorm(n = init_sample_size, mu = c(2, 3), Sigma = var.cov.mat)
  random_data <- as.data.frame(random_data)
  ggplot(random_data, aes(V1, V2)) + geom_point() + geom_smooth(method = lm, se = FALSE)
  f_init <- lm(V2 ~ V1, data = random_data)
  sigma0 <- vcov(f_init)
  mu0 <- coef(f_init)
  
  new_data_length <- new_sample_size
  new_random_data <- mvrnorm(n = new_data_length, mu = c(4, 5), Sigma = var.cov.mat)
  new_random_data <- as.data.frame(new_random_data)
  ggplot(new_random_data, aes(V1, V2)) + geom_point() + geom_smooth(method = lm, se = FALSE)
  
  #Append the data and graph
  combined_random_data_x <- append(unlist(random_data["V1"], use.names = FALSE), values = unlist(new_random_data["V1"], use.names = FALSE))
  combined_random_data_y <- append(unlist(random_data["V2"], use.names = FALSE), values = unlist(new_random_data["V2"], use.names = FALSE))
  combined_random_data <- data.frame(combined_random_data_x, combined_random_data_y)
  ggplot(combined_random_data, aes(combined_random_data_x, combined_random_data_y)) + geom_point() + geom_smooth(method = lm, se = FALSE)
  
  f_new_data <- lm(V2 ~ V1, data = new_random_data)
  sigma_new_data <- vcov(f_new_data)
  mu_new_data <- coef(f_new_data)
  
  inverse_data_length <- 1
  sum_of_new_data <- c(sum(new_random_data["V1"]), sum(new_random_data["V2"]))
  sigma_n <- sigma0%*%matrix.inverse((sigma0 + inverse_data_length*sigma_new_data))%*%(inverse_data_length*sigma_new_data)
  
  mu_n <- sigma0%*%matrix.inverse(sigma0 + inverse_data_length*sigma_new_data)%*%(inverse_data_length*mu_new_data) + inverse_data_length*sigma_new_data%*%matrix.inverse(sigma0 + inverse_data_length*sigma_new_data)%*%mu0
  
  f_final <- lm(combined_random_data_y ~ combined_random_data_x, data = combined_random_data)
  sigma_actual <- vcov(f_final)
  mu_actual <- coef(f_final)
  
  mu0
  mu_new_data
  
  mu_n
  mu_actual
  
  c(mu_n[1], mu_n[2], mu_actual[1], mu_actual[2])
}
res <- run(1000, 100)
percent_err <- c(100*abs(res[1] - res[3])/res[3], 100*abs(res[2] - res[4])/res[4])
percent_err

# for a fixed initial sample training size, say 1000, graph percent error against incoming sample size
percent_err_intercept <- c()
n_count <- c()
percent_err_coeff <- c()
for(i in seq(100, 100000, 500)){
  res <- run(i, ceiling(sqrt(i)))
  percent_err <- c(100*abs(res[1] - res[3])/res[3], 100*abs(res[2] - res[4])/res[4])
  percent_err_intercept <- c(percent_err_intercept, res[1])
  percent_err_coeff <- c(percent_err_coeff, res[2])
  n_count <- c(n_count, i)
  print(i)
}
plot(n_count, percent_err_intercept)
plot(n_count, percent_err_coeff)
ggplot(data = data.frame(n_count = n_count, percent_err_intercept = percent_err_intercept), aes(x = n_count, y = percent_err_intercept)) + 
  geom_point(aes(y = percent_err_intercept, col = 'red')) +
  stat_smooth(method = 'loess', color = 'blue') + 
  ggtitle(" ")
ggplot(data = data.frame(n_count = n_count, percent_err_coeff = percent_err_coeff), aes(x = n_count, y = percent_err_coeff)) + 
  geom_point(aes(y = percent_err_coeff, col = 'red')) +
  stat_smooth(method = 'loess', color = 'blue') + 
  ggtitle(" ")

error_df <- data.frame(n_count=n_count, percent_err_intercept=percent_err_intercept, percent_err_coeff=percent_err_coeff)
ggplot(error_df, aes(n_count, percent_err_coeff)) + 
  geom_point() + 
  stat_smooth(method = 'loess', col = 'red') +
  xlab('Initial training size') + 
  ylab('Percent Error') + 
  ggtitle('Coefficient estimate')
ggplot(error_df, aes(n_count, percent_err_intercept)) +
  geom_point() +
  stat_smooth(method = 'loess', col = 'red') + 
  xlab('Initial training size') + 
  ylab('Percent Error (intercept estimate)') +
  ggtitle('Intercept estimate')