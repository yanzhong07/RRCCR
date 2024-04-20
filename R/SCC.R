library(Matrix)
library(RSpectra)
library(igraph)

## build Knn network.
buildKNN = function(dd, k){
  set.seed(123)
  n = dim(dd)[1]
  stopifnot(dim(dd)[1] != dim(dd[2]))
  diag(dd) = 0
  edge = cbind(NA, NA, NA)
  for (i in 1:n) {
    edge_index = order(dd[i, ])[2:(k + 1)]
    dd[i,]
    
    edge = rbind(edge, cbind(i, edge_index, dd[i, edge_index]))
  }
  edge = edge[-1, ]
  edge = rbind(edge, edge[, c(2, 1, 3)])
  edge = edge[!duplicated(edge[, 1:2]), ]
  edge = edge[edge[, 1] < edge[, 2], ]
  colnames(edge) = c("i", "j", "distance")
  Ge = make_graph(c(t(edge[, 1:2])), directed = F)
  return(list(graph = Ge, edge = edge))
}

# construct H matrix for a minimum spanning tree
construct_H_mst = function(Graph0){
  n = length(V(Graph0$graph))
  Ge = Graph0$graph
  edge = Graph0$edge
  g_mst <- mst(Ge, weights = edge[,3])
  sample_order = as.vector(dfs(g_mst, root = 1)$order)
  
  zhong = g_mst[1:n,1:n]
  
  zhong = Matrix::t(zhong) + zhong
  edge_mst = Matrix::summary(zhong)
  edge_mst  = edge_mst[edge_mst[,1] < edge_mst[,2],]
  
  ## construct H matrix for mst
  H = Matrix(0, ncol =n,nrow = n-1)
  for(i in 1:dim(edge_mst)[1]){
    H[i,edge_mst[i,1]] = 1
    H[i,edge_mst[i,2]] = -1
  }  
  return(list(H = H, Ge = Ge, edge = edge, sample_order = sample_order))
}

# gic to select tuning parameters for lasso
gic.glmnet = function (x, y, crit = c("gic", "mse", "aic"), ...) {
  if (is.matrix(x) == FALSE) {
    x = as.matrix(x)
  }
  if (is.vector(y) == FALSE) {
    y = as.vector(y)
  }
  crit = match.arg(crit)
  n = length(y)
  model = glmnet(x = x, y = y, ...)
  coef = coef(model)
  lambda = model$lambda
  df = model$df
  yhat = cbind(1, x) %*% coef
  residuals = as.matrix(y - yhat)
  mse = colMeans(residuals^2)
  if(crit == "mse"){
    mse[df >= 2*n] = Inf
  }
  sse = colSums(residuals^2)
  nvar = df + 1
  gic = n * log(mse) + nvar * log(n) #* log(log(n))
  aic = n * log(mse) + 2 * log(n)
  sst = (n - 1) * var(y)
  #r2 = 1 - (sse/sst)
  #adjr2 = (1 - (1 - r2) * (n - 1)/(nrow(x) - nvar - 1))
  crit = switch(crit, gic = gic, mse = mse, aic = aic)
  selected = best.model = which(crit == min(crit))
  ic = c(gic = gic[selected], mse = mse[selected])
  result = list(coefficients = coef[, selected], ic = ic, lambda = lambda[selected], 
                nvar = nvar[selected], glmnet = model, residuals = residuals[,  selected], fitted.values = yhat[, selected], ic.range = crit, 
                call = match.call())
  class(result) = "gic.glmnet"
  return(result)
}


SCC = function(X, y, loc, K = 5, intercept = F, crit = "gic"){
  library(glmnet)
  dd = dist(loc)
  while(T){
    Graph0 = buildKNN(as.matrix(dd), k = K)
    if(is.connected(Graph0$graph)){
      break
    }
    K = K+1
  }

  n <- length(y)
  p <- dim(X)[2]
  HGe = construct_H_mst(Graph0)
  H = HGe$H
  #Matrix(HGe$H, sparse = T) # make H to sparse matrix
  Htilde <- rbind(H, 1/n)
  inHtilde <- solve(Htilde)
  inHtilde <- as.matrix(inHtilde)
  Xnew1 = NA
  for(i in 1:p){
    Xnew1 = cbind(Xnew1, diag(c(X[,i])) %*% inHtilde)
  }
  Xnew1 = Xnew1[,-1]
  penfac = rep(1,p*n)
  penfac[seq(n,p*n, by = n)] = 0
  
  Cmodel <- gic.glmnet(x = Xnew1, y = y, crit = crit, intercept = intercept, alpha = 1, family = "gaussian", standardize = F, penalty.factor = penfac, standardize.response = F)
  C = matrix(Cmodel$coefficients[-1], nrow = n)
  B = inHtilde %*% C
  a = Cmodel$coefficients[1]
  
  return(list(B = B, a = a, H =H)) 
}
