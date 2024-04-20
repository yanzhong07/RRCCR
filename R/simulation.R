set_A = function(x,y){
  a1 = c(2.5, -2, -1,1, 0.5) 
  a2 = c(0, 1, -2,0.5, -1) 
  a3 = c(1, 1.5, -0.5, -1, 2) 
  
  group = NA
  a20 = cbind(a1,a2,a3)
  
  if(x < 0.5 & y < 0.5){
    a_use = a20[1,]
    group = 1
  }else if(x < 0.5 & y > 0.5){
    a_use = a20[2,]
    group = 2
    
  }else if(x > 0.5 & y < 0.5){
    a_use = a20[3,]
    group = 3
    
  }else{
    a_use = a20[4,]
    group = 4
    
  }
  if( sqrt(((x-0.5)^2 + (y-0.5)^2)) < 0.25){
    a_use = a20[5,]
    group = 5
    
  }
  return(list(a_use = a_use, group = group))
}


build_A = function(loc){
  n = dim(loc)[1]
  A = matrix(0, ncol = 3, nrow = n)
  group = rep(0,n)
  for(i in 1:n){
    tmp = set_A(loc[i,1], loc[i,2])
    group[i] = tmp$group
    A[i,] = tmp$a_use
  }
  
  a1 = c(2.5, -2, -1,1, 0.5) 
  a2 = c(0, 1, -2,0.5, -1) 
  a3 = c(1, 1.5, -0.5, -1, 2) 
  
  # a1 %*% a2
  # a1 %*% a3
  # a2 %*% a3
  
  a_coef = cbind(a1,a2,a3)
  return(list(A = A, a_coef = a_coef, group = group))
}

#cor_selection = 1,2,3, three selections for multicolinearity.
#
build_simulation = function(n = 1000, p = 10, rho = 0.5, cor_selection = 1, sigma.error = 0.25, seed = 123){
  set.seed(seed)
  library(mvtnorm)
  library(Matrix)
  r = 3
  locations = matrix(runif(n*2),nrow = n, ncol = 2) # generate locations from [0,1]* [0,1] space
  
  ## construct B
  dsample = matrix(rnorm(p*r),ncol = 3)
  Pm = svd(dsample)$u
  
  tmp = build_A(locations)
  A = tmp$A * sqrt(p/10)
  group = tmp$group
  a_coef = tmp$a_coef
  B = A %*% t(Pm)
  
  # construct X
  dd = as.matrix(dist(locations))
  covm = exp(-dd/rho)
  Z = NA
  for(i in 1:p){
    tmp = mvtnorm::rmvnorm(1,mean = rep(0,n), sigma = covm )
    Z = cbind(Z, c(tmp))
  }
  Z = as.matrix(Z[,-1])
  
  if(cor_selection == 1){
    Sigma = matrix(0,p,p)
    for(i in 1:p){
      for(j in 1:p){
        Sigma[i,j] = 0.3^(abs(i-j))  # this is an important parameter
      }
    }
    choS =chol(Sigma)
    X = Z %*% choS 
  }else if(cor_selection == 2){
    Sigma = matrix(0,p,p)
    for(i in 1:p){
      for(j in 1:p){
        Sigma[i,j] = 0.8^(abs(i-j))  # this is an important parameter
      }
    }
    choS =chol(Sigma)
    X = Z %*% choS 
  }else{
    
    Sigma1 = matrix(0,p,p)
    Sigma2 = matrix(0,p,p)
    for(i in 1:p){
      for(j in 1:p){
        Sigma1[i,j] = 0.8^(abs(i-j))  # this is an important parameter
        Sigma2[i,j] = 0.3^(abs(i-j))  # this is an important parameter
      }
    }
    
    choS1 =chol(Sigma1)
    choS2 =chol(Sigma2)
    X = Z  
    X[group %in% c(1,2,3),] = Z[group %in% c(1,2,3),] %*% choS2
    X[group %in% c(4,5),] = Z[group %in% c(4,5),] %*% choS1
  }
  
  y = rowSums(X * B) + rnorm(n,0,sigma.error)
  return(list(y = y, X = X,locations = locations, B = B, A = A, Pm = Pm, a_coef = a_coef, group = group, cor_selection = cor_selection, phi = 0.5, sigma.error = sigma.error))
}

