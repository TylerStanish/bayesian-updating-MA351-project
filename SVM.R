library(ggplot2)
library(matrixcalc)
library(MASS)
library(mvtnorm)
library(LearnBayes)
library(nlme)

beer_data <- read.table("recipeData.csv", header = TRUE)
# this function will always have the output because I split the data in order. Everything before the cutoff point is taken
# from the beginning of the data set, and everything after is taken from the end
svm <- function(cutoff_point){
  beer_data$FG <- normalize.vector(beer_data$FG)
  f_init <- glm(formula = FG ~ BoilSize, family = gaussian, data = beer_data[1:(nrow(beer_data)-cutoff_point),])
  sigma0 <- vcov(f_init)
  mu0 <- coef(f_init)
  
  f_new_data <- glm(formula = FG ~ BoilSize, family = gaussian, data = beer_data[(nrow(beer_data)-cutoff_point):nrow(beer_data),])
  sigma_new_data <- vcov(f_new_data)
  mu_new_data <- coef(f_new_data)
  
  inverse_data_length <- 1
  sigma_n <- sigma0%*%matrix.inverse((sigma0 + inverse_data_length*sigma_new_data))%*%(inverse_data_length*sigma_new_data)
  
  mu_n <- sigma0%*%matrix.inverse(sigma0 + inverse_data_length*sigma_new_data)%*%(inverse_data_length*mu_new_data) + inverse_data_length*sigma_new_data%*%matrix.inverse(sigma0 + inverse_data_length*sigma_new_data)%*%mu0
  
  f_final <- glm(formula = FG ~ BoilSize, family = gaussian, data = beer_data)
  sigma_actual <- vcov(f_final)
  mu_actual <- coef(f_final)
  
  mu0
  mu_new_data
  
  mu_n
  mu_actual
  c(mu_n[1], mu_n[2], mu_actual[1], mu_actual[2])
}
# I choose the square root here and for the demonstrations for the OLS case
res_nonlinear <- svm(ceiling(sqrt(nrow(beer_data))))
percent_err <- c(100*abs(res_nonlinear[1] - res_nonlinear[3])/res_nonlinear[3], 100*abs(res_nonlinear[2] - res_nonlinear[4])/res_nonlinear[4])
percent_err