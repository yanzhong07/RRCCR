
\name{TreeRRCCR}
\alias{TreeRRCCR}

\title{TreeRRCCR}
\usage{

TreeRRCCR(X, y, H = NULL, r, Bprior, lambda_list = exp(-5), lambda_list_initial = lambda_list, tol = 10^(-7), max.iter = 100, A_initial = NULL, P_initial = NULL)

}
\arguments{
\item{X}{An n*p matrix for n observations with p covariates. Set the first column of X as 1 if you want to include the varying intercept term.}

\item{y}{The response variable for n observations.}

\item{H}{The incidence matrix for a tree graph }

\item{r}{The number of rank used.}

\item{Bprior}{The initial estimator of B for calculating the adaptive weight.}

\item{lambda_list}{A list of penalty parameters.}

\item{tol}{Tolerance for convergence.}

\item{max.iter}{Max number of iterations.}

\item{A_initial, P_initial}{Initial value of A and P.}

\item{alpha}{The allocation of adaptive distance. 0 means that only use the orignal auxillary information}

}
\value{
Return a list of the estimators of B with different lambda, and a vector of lambda used.
}
\description{
This function is designed for the reduced-rank clustered coefficient regression with the edge set of a tree graph as E.
}
