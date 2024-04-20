#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(cpp11)]]
//#include <Rcpp.h>
#include <vector>
#include <math.h>
using namespace Rcpp;
using namespace std;


// X: p*n B:k*n Y:n*1 A:k*r theta: p*r
// calculate residual of Y


arma::vec Residual(const arma::vec *y,const arma::mat *X,const arma::mat *A,const double *b0,const arma::mat *D){
  unsigned int i=0;
  arma::vec diff(X->n_rows,arma::fill::zeros);
  for(i=0;i<X->n_rows;++i){
    diff.row(i) = (X->row(i))*(*D)*(A->row(i).t());
  } 
  diff = diff + *b0;
  diff=(*y - diff);
  return diff;
}

// calculate the gradient of A

arma::mat gradD(const arma::vec *y,const arma::mat *X,const arma::mat *A,const double *b0,const arma::mat *D){
  unsigned int i;
  arma::vec diff=Residual(y,X,A,b0,D);
  arma::mat grad(D->n_rows,D->n_cols,arma::fill::zeros);   //K*r
  arma::mat grad_temp(D->n_rows,D->n_rows,arma::fill::zeros);  // K *p
  //arma::mat grad_tempSum(D->n_rows,D->n_rows,arma::fill::zeros);
  //arma::mat firstPart(A->n_rows,1,arma::fill::zeros);      //K
  
  for(i=0;i<y->n_elem;++i){
    grad_temp = diff(i) * (X->row(i).t()) * (A->row(i));
    grad -= grad_temp;
    // grad=grad+(-2.0)*(B->col(i))*(Y->col(i)-arma::kron((X->col(i)).t(),I_q)*(*Theta)*A->t()*B->col(i)).t()*(arma::kron((X->col(i)).t(),I_q))*(*Theta);
  }
    grad=grad / y->n_elem;  // not -2
  return(grad);
}

// calculate the projection of A
arma::mat projectionD(arma::mat *gradF,const arma::mat *D){
  arma::mat grad_normal=D->t()*(*gradF);
  grad_normal=0.5*(*D)*(grad_normal+grad_normal.t());
  arma::mat grad=*gradF-grad_normal;
  return(grad);
}

// calculate QR 
arma::mat retractD(const arma::mat *direc,const double* stepSize,const arma::mat *D){
  arma::mat Dt=*D - *stepSize*(*direc); // the sign of this is negative.
  arma::mat retract_Q,retract_R;
  arma::qr_econ(retract_Q,retract_R,Dt);
  if(retract_R(0,0)<0){
    retract_Q=-retract_Q;
  }
  return(retract_Q);
}


// square error of regression
double totalloss(const arma::vec *y,const arma::mat *X,const arma::mat *A,const double *b0,const arma::mat *D){
  double obj=0.0;
  int n = y->n_elem;
  arma::vec diff=Residual(y,X,A,b0,D);
  obj = arma::dot(diff,diff) /2 / n;
  return obj;
}

//[[Rcpp::export]]
arma::mat updateD(const arma::vec& y, const arma::mat& X, const arma::mat& A, const double& b0, const arma::mat& D_initial,
                        double tol = 0.0000001,
                        double MaxIt = 100,
                        double subMaxIt = 10,
                        double epsilon = 0.5){
  // Initialization
  //int n = X.n_rows;
  int r = D_initial.n_cols;
  int p = D_initial.n_rows;
  arma::mat D(p,r,arma::fill::zeros);
  D = D_initial;

  // controls
  int count=0,subcount=0;
  arma::mat *D_temp=new arma::mat;
  arma::mat *steepD=new arma::mat;
  double TS,TS_temp,eAscent,TS_last,TS_Des, tol_now, tol_now2;
  double stepSize=1;
  TS = totalloss(&y,&X,&A,&b0,&D); // goal function
  
  arma::mat D_last;
  vector<double> save_TS;
  vector<arma::mat> save_D;
  *D_temp=D;
  bool flag=true;
  //#Rcout<<"Before Updating D, the objetive value is "<<TS<<endl;
  // Steepest Ascent and Armijo back tracking
  
  do{
    // control of inner Armijo back tracking
    count++;
    TS_Des=-1.0;
    TS_last=TS; // save TS and D for the last time
    D_last=D;
    // gradient.
    *steepD=gradD(&y,&X,&A,&b0,&D);
    *steepD=projectionD(steepD,&D);
    eAscent=stepSize*epsilon*arma::dot(*steepD,*steepD); //a kind of tolerance
    
    // this innerIter is used to save the power of step size "beta=0.5".
    // Find the best step of gradient
    subcount=0;
    tol_now=-1.0;
    do{
      subcount++;
      stepSize=0.5*stepSize;  //step size: beta = 0.5
      eAscent=0.5*eAscent;
      *D_temp=retractD(steepD,&stepSize,&D);
      TS_temp=totalloss(&y,&X,&A,&b0,D_temp);
      save_TS.push_back(TS_temp);
      save_D.push_back(*D_temp);
     // Rcout<<"The objetive value now is "<<TS_temp<<endl;
     //   Rcout<<"The objetive value now is "<<TS_temp<<endl;
     //   Rcout<<"TS_last-TS_temp now is "<<TS_last-TS_temp<<endl;
     //   Rcout<<"eAscent is "<<eAscent<<endl;

    }while (subcount<subMaxIt && (TS_last-TS_temp)<eAscent ); //out of Armijo
    
    TS_Des=TS_last-TS_temp;
    if(TS_Des<0.0) {
      TS=*min_element(save_TS.begin(),save_TS.end());
      D=save_D[distance(save_TS.begin(),min_element(save_TS.begin(),save_TS.end()))];
    } else{
      D=*D_temp;
      TS=totalloss(&y,&X,&A,&b0,D_temp);
    }
    // clean the container
    vector<double>().swap(save_TS);
    vector<arma::mat>().swap(save_D);
    
    tol_now=arma::norm((D-D_last),"fro");
    tol_now2=abs(TS_last-TS);
    if(tol_now < tol || tol_now2 < tol) flag=false;
    
  }while (count<MaxIt && flag); // Out of iteration
  delete D_temp;
  delete steepD;
  //# Rcout<<"After Updating D, the objetive value is "<<TS<<endl;
 // Rcout<<"After Updating D, the iteration is "<<count<<endl;
 // Rcout<<"After Updating D, the sub iteration is "<<subcount<<endl;
  return(D);
}


