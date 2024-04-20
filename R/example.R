
source("./R/simulation.R")
source("./R/SCC.R")
source("./R/RRCCR.R")

################################################################3
# This example code is build only for the case that E is the edge set of a tree graph.



# cor_selection: 1. low level of multicolineariy 2. high level of multicolineariy 3. Local multicolineariy
sim1 = build_simulation(n = 1000, p = 10, cor_selection = 3, seed = 123)
X = sim1$X
y = sim1$y
loc = sim1$locations
dd = dist(sim1$locations)

Btrue = sim1$B
groupT = sim1$group

# SCC to generate graph and initial estimator
t0 = Sys.time()
compare_scc = SCC(X, y, loc = loc, K = 5, intercept = F, crit = "gic")
Huse = compare_scc$H  # E for the following study.
Bscc = as.matrix(compare_scc$B)


# Two step SCC to built the initialization
X_svd_full = svd(X)
kuse = max(which.max(cumsum(X_svd_full$d^2)/sum(X_svd_full$d^2) > 0.9), 7)
X_svd = svds(X, k = kuse)
compare_2s = SCC(X = X_svd$u %*% diag(X_svd$d), y, loc = loc, K = 5, intercept = F, crit = "gic")
B2s = t(X_svd$v %*% t(compare_2s$B))


## The proposed method
B2s_svd = svds(B2s,3)
A2s_use = B2s_svd$u %*% diag(B2s_svd$d)
D2s_use = B2s_svd$v

my1 = TreeRRCCR(X, y, H = Huse, r = 3, Bprior = Bscc, lambda_list = exp(-5), tol = 10^(-7), max.iter = 100, A_initial = A2s_use, P_initial = D2s_use)
## take around 5-10 minutes. The computing time may vary depending on the computer used.
Bmy = my1$B[[1]]

## comparsion

# Estimation Error
sqrt(mean((Bscc - Btrue)^2))
sqrt(mean((B2s - Btrue)^2))
sqrt(mean((Bmy - Btrue)^2))

cor(Bscc[,1],Btrue[,1])
cor(B2s[,1], Btrue[,1])
cor(Bmy[,1],Btrue[,1])

