library(Rcpp)
library(RcppArmadillo)
library(Matrix)
library(RSpectra)
library(igraph)
library(sparsegl)

sourceCpp("./src/manifold_gradient_descent.cpp")

# This code is build only for the case that E is the edge set of a tree graph.
# H is the incidence matrix for the edge set.
TreeRRCCR= function(X, y, H = NULL, r, Bprior, lambda_list = exp(-5), tol = 10^(-7), max.iter = 100, A_initial = NULL, P_initial = NULL){
  n = length(y)
  p = dim(X)[2]
  
  ## step 1: build Htilde matrix
  Htilde <- rbind(H, 1/n)
  iHtilde <- as.matrix(solve(Htilde))
  
  # initializating
  if(is.null(A_initial)){
    At = matrix(rnorm(n*r), nrow = n, ncol = r)
  }else{
    stopifnot(dim(A_initial)[1] == n, dim(A_initial)[2] == (r))
    At = A_initial
  }
  if(is.null(P_initial)){
    Pt = matrix(rnorm(p*r), nrow = p, ncol = r)
    Pt = qr(Pt)
    Pt = qr.Q(Pt)
  }else{
    stopifnot(dim(P_initial)[1] == p, dim(P_initial)[2] == r)
    Pt = P_initial
  }
  At_1 = At
  Pt_1 = Pt
  
  lambda_list = sort(lambda_list, decreasing = T)
  nlambda = length(lambda_list)
  B = list()

  Bpriortilde = Htilde %*% Bprior
  weight_list = (1/sqrt(rowSums(Bpriortilde[-n,]^2)))
  weight_list  = weight_list/min(weight_list) * 0.1
  weight_list[weight_list > 5] = 5
  
  print("Begin Modeling.")
  for(ilambda in 1:nlambda){
    diff = Inf
    count = 0
    t0 = Sys.time()
    while(diff >tol & count < max.iter){
      # fix Pt, update At and b0: proximal gradient descent
      Xstar = X %*% Pt
      Xtilde = matrix(0,nrow = n,n*r)
      for(i in 1:r){
        Xtilde[,c(1:n)*r-r+i] = iHtilde * Xstar[,i]
      }
      groups <- rep(1:n, each = r)
      fit <- sparsegl(Xtilde, y, group = groups, asparse = 0, standardize = FALSE, pf = c(weight_list,0), intercept = F,lambda = lambda_list[ilambda])
      At = matrix(fit$beta, ncol = r, byrow = T)
      At_back = as.matrix(iHtilde %*% At) # Atilde back to original A
      Pt = updateD(y=y, X = X, A = At_back, b0 = 0, D_initial = Pt_1, tol = tol, MaxIt = max.iter, subMaxIt = 10, epsilon = 0.5)

      count = count + 1
      diff = sqrt(sum((At %*% t(Pt) - At_1 %*% t(Pt_1))^2) / n / p)
      At_1 = At
      Pt_1 = Pt
    }
    t1 = Sys.time()
    print(paste0("Finish",ilambda,". Computing Time: ", difftime(t1,t0,units = "mins")))
    B[[ilambda]] = At_back %*% t(Pt)
  }
  print("Model Finish.")
  
  return(list(B = B, lambda_list = lambda_list))
}



