rm(list=ls())
source("Required_functions.R")

load("temperature_daily.Rda")
load("wind_daily.Rda")
load("solar_daily.Rda")

###Predetermined number of basis functions and roughness parameter
smoothing_params_Y = c(150, -2.333333) #solar
smoothing_params_X1 = c(147, -4.777778) #temperature
smoothing_params_X2 = c(62, 0.1111111) #wind
smoothing_params_X = list(smoothing_params_X1, smoothing_params_X2)

######################################################
### N       : sample size
### d_time  : discrete time points
### n       : number of design points 

N = 70
n = 365
d_time = c(1:365)
######################################################

###Data 
Y = as.data.frame(my.solar)
predictor_list = list(as.data.frame(my.temp), as.data.frame(my.wind))

###Number of predictors
np = length(predictor_list)

###Smooth the data using the optimal parameters
S_Y = smooth_all(data = Y, param_vec = smoothing_params_Y, d_time = d_time)
S_X = vector("list",)
for(i in 1:np){
  S_X[[i]] = smooth_all(data = predictor_list[[i]], param_vec = smoothing_params_X[[i]], d_time = d_time)
}

###Basis functions
Phi_y = Bsplines_FDA(d_time = d_time, nbf = smoothing_params_Y[1], norder=4)
Phi_x = vector("list",)
for(i in 1:np){
  Phi_x[[i]] = Bsplines_FDA(d_time = d_time, nbf = smoothing_params_X[[i]][1], norder=4)
}

###Observed response functions
solar_fun = S_Y%*%t(Phi_y)
#ts.plot(t(solar_fun))

train_stations = c(1:50) #training stations
test_stations = c(51:70) #test stations

###Inner product matrices
Inn_prod = vector("list",)
for(i in 1:np){
  m_basis = create.bspline.basis(rangeval = range(d_time), nbasis = smoothing_params_X[[i]][1], norder=4)
  Inn_prod[[i]] = inprod(m_basis, m_basis)
}


###Center the response and predictors for training stations
#Response
CY = matrix(NA, nrow = length(train_stations), ncol = ncol(S_Y))
for(i in 1:length(train_stations))
  CY[i,] = S_Y[train_stations,][i,] - apply(S_Y[train_stations,], 2, mean)

#Predictors
CX = vector("list",)
for(j in 1:np){
  cen_x = matrix(NA, nrow = length(train_stations), ncol = ncol(S_X[[j]]))
  for(i in 1:length(train_stations))
    cen_x[i,] = S_X[[j]][train_stations,][i,] - apply(S_X[[j]][train_stations,], 2, mean)
  CX[[j]] = cen_x
}

###Center the predictors for test stations
CX_test = vector("list",)
for(j in 1:np){
  cen_x_test = matrix(NA, nrow = length(test_stations), ncol = ncol(S_X[[j]]))
  for(i in 1:length(test_stations))
    cen_x_test[i,] = S_X[[j]][test_stations,][i,] - apply(S_X[[j]][test_stations,], 2, mean)
  CX_test[[j]] = cen_x_test
}

###Z matrix for the regression
#Training sample
Z = NULL
for(i in 1:np){
  Z = cbind(Z, CX[[i]] %*% Inn_prod[[i]])
}

#Test sample
Z_test = NULL
for(i in 1:np){
  Z_test = cbind(Z_test, CX_test[[i]] %*% Inn_prod[[i]])
}

###SIMPLS
model_simpls = pls.regression(Xtrain = Z, Ytrain = CY, Xtest = Z,
                            ncomp = 5, unit.weights=FALSE)
Bhat_simpls = model_simpls$B


###NIPALS
model_nipals = plsreg2(predictors = Z, responses = CY, comps = 5,crosval = TRUE)
Bhat_nipals = model_nipals$reg.coefs[1:(smoothing_params_X1[1]+smoothing_params_X2[1]),]


###Obtain fitted coefficients
fitted_coefs_simpls = Z_test%*%Bhat_simpls
fitted_coefs_nipals = Z_test%*%Bhat_nipals

for(i in 1: length(test_stations)){
  fitted_coefs_simpls[i,] = fitted_coefs_simpls[i,] + apply(S_Y[train_stations,], 2, mean)
  fitted_coefs_nipals[i,] = fitted_coefs_nipals[i,] + apply(S_Y[train_stations,], 2, mean)
}

###Obtain fitted functions
fitted_functions_simpls = fitted_coefs_simpls%*%t(Phi_y)
fitted_functions_nipals = fitted_coefs_nipals%*%t(Phi_y)

ts.plot(t(solar_fun[test_stations,])) # Time series plot of the observed solar radiation functions in the test sample
ts.plot(t(fitted_functions_simpls)) # Time series plot of the fitted solar radiation functions in the test sample (SIMPLS)
ts.plot(t(fitted_functions_nipals)) # Time series plot of the fitted solar radiation functions in the test sample (NIPALS)


sum((fitted_functions_simpls - solar_fun[test_stations,])^2)/length(test_stations) #AMSE_p for SIMPLS
sum((fitted_functions_nipals - solar_fun[test_stations,])^2)/length(test_stations) #AMSE_p for NIPALS

