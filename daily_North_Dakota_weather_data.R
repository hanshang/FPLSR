rm(list=ls())
source("auxiliary_functions.R")

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
m_basis_y = create.bspline.basis(rangeval = range(d_time), nbasis = smoothing_params_Y[1], norder=4)
Inn_prod_y = inprod(m_basis_y, m_basis_y)

Inn_prod_x = vector("list",)
for(i in 1:np){
  m_basis_x = create.bspline.basis(rangeval = range(d_time), nbasis = smoothing_params_X[[i]][1], norder=4)
  Inn_prod_x[[i]] = inprod(m_basis_x, m_basis_x)
}

###Square roots of the inner product matrices
Inn_prod_y_sqrt = sqrtm(Inn_prod_y)
Inn_prod_x_sqrt = vector("list",)

for(i in 1:np)
  Inn_prod_x_sqrt[[i]] = sqrtm(Inn_prod_x[[i]])

###Arguments for the PLS regression
#Response
Reg_y = Reg_y_out = S_Y %*% Inn_prod_y_sqrt

#Predictors
Reg_mat = vector("list",)

for(i in 1:np){
  Reg_mat[[i]] = S_X[[i]] %*% Inn_prod_x_sqrt[[i]]
}

###Divide the data into training and test samples
#Response
Reg_y_train = Reg_y[train_stations,]

#Predictors
Reg_mat_train = vector("list",)
Reg_mat_test = vector("list",)

for(i in 1:np){
  Reg_mat_train[[i]] = Reg_mat[[i]][train_stations,]
  Reg_mat_test[[i]] = Reg_mat[[i]][test_stations,]
}

Reg_mat_train = do.call(cbind, Reg_mat_train)
Reg_mat_test = do.call(cbind, Reg_mat_test)

###Mean of variables
mean_y = apply(Reg_y_train, 2, mean)

mean_x = apply(Reg_mat_train, 2, mean)
mean_x_test = apply(Reg_mat_test, 2, mean)


###Center the response and predictors for training stations
#Response
CY = matrix(NA, nrow = nrow(Reg_y_train), ncol = ncol(Reg_y_train))
for(i in 1:nrow(Reg_y_train))
  CY[i,] = Reg_y_train[i,] - mean_y

#Predictors
CX = matrix(NA, nrow = nrow(Reg_mat_train), ncol = ncol(Reg_mat_train))
for(i in 1:nrow(Reg_mat_train))
  CX[i,] = Reg_mat_train[i,] - mean_x

###Center the predictors for test stations
CX_test = matrix(NA, nrow = nrow(Reg_mat_test), ncol = ncol(Reg_mat_test))
for(i in 1:nrow(Reg_mat_test))
  CX_test[i,] = Reg_mat_test[i,] - mean_x_test

###SIMPLS
model_simpls = pls.regression(Xtrain = CX, Ytrain = CY, Xtest = CX,
                            ncomp = 5, unit.weights=FALSE)
Bhat_simpls = model_simpls$B


###NIPALS
model_nipals = plsreg2(predictors = CX, responses = CY, comps = 5,crosval = TRUE)
Bhat_nipals = model_nipals$reg.coefs[1:(smoothing_params_X1[1]+smoothing_params_X2[1]),]


###Obtain fitted coefficients
fitted_coefs_simpls = CX_test%*%Bhat_simpls
fitted_coefs_nipals = CX_test%*%Bhat_nipals

for(i in 1: length(test_stations)){
  fitted_coefs_simpls[i,] = fitted_coefs_simpls[i,] + mean_y
  fitted_coefs_nipals[i,] = fitted_coefs_nipals[i,] + mean_y
}

###Obtain fitted functions
fitted_functions_simpls = (fitted_coefs_simpls %*% solve(Inn_prod_y_sqrt))%*%t(Phi_y)
fitted_functions_nipals = (fitted_coefs_nipals %*% solve(Inn_prod_y_sqrt))%*%t(Phi_y)

ts.plot(t(solar_fun[test_stations,])) # Time series plot of the observed solar radiation functions in the test sample
ts.plot(t(fitted_functions_simpls)) # Time series plot of the fitted solar radiation functions in the test sample (SIMPLS)
ts.plot(t(fitted_functions_nipals)) # Time series plot of the fitted solar radiation functions in the test sample (NIPALS)


sum((fitted_functions_simpls - solar_fun[test_stations,])^2)/length(test_stations) #AMSE_p for SIMPLS
sum((fitted_functions_nipals - solar_fun[test_stations,])^2)/length(test_stations) #AMSE_p for NIPALS

