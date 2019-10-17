#################################### 
###         Functions            ###
####################################
install.packages(c("fda", "plsgenomics", "plsdepot"))
library(fda)
library(plsgenomics)
library(plsdepot)

### Generate Matrix of Bsplines basis functions using the FDA package
Bsplines_FDA = function(d_time, nbf, norder=4){
  require(fda)
  basis = create.bspline.basis(rangeval = range(d_time), nbasis = nbf, norder)
  Phi = eval.basis(evalarg = d_time, basisobj = basis)
  return(Phi)
}

###Smooth all the FD by param_vec
smooth_all = function(data, param_vec, d_time){
  
  ### Generate Matrix of Bsplines basis functions using the FDA package
  Bsplines_FDA = function(d_time, nbf, norder=4){
    require(fda)
    basis = create.bspline.basis(rangeval = range(d_time), nbasis = nbf, norder)
    Phi = eval.basis(evalarg = d_time, basisobj = basis)
    return(Phi)
  }
  
  ###Penalized Maximum Likelihood function
  PML = function(Phi, n, lambda, data){
    D = matrix(0,(n-2),n)
    D[1, ] = c(1,-2,1,rep(0,(n-3)))
    for (i in 1:(n-4)) {
      D[(i+1), ] = c(rep(0,i),1,-2,1,rep(0,(n-3)-i))
    }
    D[(n-2), ] = c(rep(0,(n-3)),1,-2,1)
    K = t(D)%*%D
    
    lamda = 10^(lambda)
    sigma = 2
    sigma1 = 1
    
    while((sigma-sigma1)^2 > 1e-7){
      Phi_inv = try(
        solve(t(Phi)%*%Phi+length(data)*(lamda)*(sigma)*K,diag(ncol(K))),silent = T)
      weight = (Phi_inv)%*%t(Phi)%*%matrix(data, ncol=1)
      sigma1 = sigma
      sigma1 = as.vector(sigma1)
      sigma = (1/length(data))*t(matrix(data, ncol=1)-Phi%*%weight)%*%(matrix(data, ncol=1)-Phi%*%weight)
      sigma = as.vector(sigma)       
    }
    list(lamda=lamda,sigma=sigma,K=K,weight=weight)
  }
  
  smoothed_data = matrix(0,nrow(data), param_vec[1])
  for(i in 1:nrow(data)) {
    my_data = as.matrix(data[i,])
    Phi = Bsplines_FDA(d_time = d_time, nbf = param_vec[1], norder=4)
    smoothed_data[i,] = PML(Phi = Phi, n = param_vec[1], lambda = param_vec[2], data = my_data)$weight
  }
  return(smoothed_data)
}
