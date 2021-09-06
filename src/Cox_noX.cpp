#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <ctime>
#include "BS.h"
#include "compare.h"

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// U: score function; 
// I: information matrix; 
// gamma: hazard ratio
// time: event time (already sorted); 
// t: unique failure time; 
// delta: event status;
// W: entry time; 
// S: vaccination time; 
// knots: change points of hazard ratio;
// constantVE: indicator of constant VE
double getUI_noX(arma::vec& U, arma::mat& I, 
                 const arma::vec& gamma, const arma::vec& time, 
                 const arma::vec& t, const arma::ivec& delta, 
                 const arma::vec& W, const arma::vec& S, 
                 const arma::vec& knots, bool constantVE) {

  int k = t.n_elem, n = delta.n_elem, p2 = gamma.n_elem, di;
  double ti, xbeta, exbeta, ll = 0.0, SR0, SD0, tempS0;
  arma::vec SR1, SD1, tempS1;
  arma::rowvec tmpe;
  arma::mat SR2, SD2, tempS2;

  for (int i=k-1; i>=0; --i) {
    ti = t(i);
    di = 0;
    // SR: sum in the risk set
    // SD: sum in the failure set
    SR0 = SD0 = 0.0; 
    SR1.zeros(p2); 
    SD1.zeros(p2);
    SR2.zeros(p2, p2); 
    SD2.zeros(p2, p2);
    for (int l = n-1; l >= 0; --l) {
      if (!ge(time(l), ti)) break;
      if (ge(ti, W(l))) {
        if (S(l)<ti) {
          tmpe = BS(ti-S(l), knots, constantVE);
          xbeta = arma::as_scalar(tmpe*gamma);
          exbeta = exp(xbeta);
          SR0 += exbeta;
          SR1 += exbeta*tmpe.t();
          SR2 += exbeta*tmpe.t()*tmpe;

          if (eq(time(l),ti) && delta(l) == 1) {
            di++;
            SD0 += exbeta;
            SD1 += exbeta*tmpe.t();
            SD2 += exbeta*tmpe.t()*tmpe;
            U += tmpe.t();
            ll += xbeta;
          }

        } else {
          tmpe.zeros(p2);
          exbeta = 1.0;
          SR0 += exbeta;

          if (eq(time(l),ti) && delta(l) == 1) {
            di++;
            SD0 += exbeta;
          }
        }
      }
    } // end l
    
    for (int j = 0; j < di; ++j) {
      tempS0 = SR0 - j*SD0/di;
      tempS1 = SR1 - j*SD1/di;
      tempS2 = SR2 - j*SD2/di;
      U -= tempS1/tempS0;
      I += tempS2/tempS0 - tempS1*tempS1.t()/pow(tempS0, 2);
      ll -= log(tempS0);
    }
  } // end i
  return ll;
} 

arma::mat getw_noX(const arma::vec& gamma, 
                   const arma::vec& time, const arma::vec& t, 
                   const arma::ivec& delta, const arma::vec& W, 
                   const arma::vec& S, const arma::vec& knots, bool constantVE) {

  // compute w matrix for robust variance estimator
  int k = t.n_elem, n = delta.n_elem, p2 = gamma.n_elem;
  int di;
  double exbeta, S0, SR0, SD0;
  arma::vec SR1, SD1, S1;
  arma::rowvec tmpe;
  arma::mat w;
  
  w.zeros(n, p2);
  
  for (int q = k-1; q >= 0; --q) {

    di = 0;
    SR0 = 0.0; 
    SD0 = 0.0;
    SR1.zeros(p2); 
    SD1.zeros(p2); 

    for (int l = n-1; l >= 0; --l) {
      if (!ge(time(l), t(q))) break;
      if (ge(t(q), W(l))) {
        if (S(l)<t(q)) {
          tmpe = BS(t(q)-S(l), knots, constantVE);
          exbeta = exp(arma::as_scalar(tmpe*gamma));
          SR0 += exbeta;
          SR1 += exbeta*tmpe.t();
          if (eq(time(l),t(q)) && delta(l) == 1) {
            di++;
            SD0 += exbeta;
            SD1 += exbeta*tmpe.t();
          }
        } else {
          tmpe.zeros(p2);
          exbeta = 1.0;
          SR0 += exbeta;
          if (eq(time(l),t(q)) && delta(l)==1) {
            di++;
            SD0 += exbeta;
          }
        }
        
      }
    }
  
    for (int l = n-1; l >= 0; --l) {
      if(!ge(time(l), t(q))) break;
      if(ge(t(q), W(l))) {
        if (S(l)<t(q)) {
          tmpe = BS(t(q)-S(l), knots, constantVE);
          exbeta = exp(arma::as_scalar(tmpe*gamma));
        } else {
          tmpe.zeros(p2);
          exbeta = 1.0;
        }

        if(eq(time(l), t(q)) && delta(l)==1) {
          for (int j = 0; j < di; ++j) {
            S0 = SR0 - j*SD0/di;
            S1 = SR1 - j*SD1/di;
            w.row(l) += (tmpe-S1.t()/S0)/di;
            w.row(l) -= (tmpe-S1.t()/S0)*(di-j)*exbeta/n/S0/di;
          }
        } else {
          for (int j = 0; j < di; ++j) {
            S0 = SR0 - j*SD0/di;
            S1 = SR1 - j*SD1/di;
            w.row(l) -= (tmpe-S1.t()/S0)*exbeta/n/S0;
          }
        }
      }
    }
  }
  return w;
}

// [[Rcpp::export]]
List Cox_noX(arma::vec& gamma, const arma::vec& time, const arma::ivec& delta, 
             const arma::vec& W, const arma::vec& S, 
             const arma::vec& knots, 
             bool constantVE, double threshold=10^(-4), int maxit=500) {
  // Standard Cox model with piecewise log-linear hazard ratio r(t) with 5 pieces
  int p2 = gamma.n_elem, it = 0;
  double dif = 1.0, ll;
  arma::vec gamma0, U, temp, t = unique(time(find(delta)));
  arma::mat I, I_inv, w, covar(p2, p2);
  bool flag; // indicator whether I is invertible
  
  // begin: Newton-Raphson
  
  while (it<maxit && dif>threshold) {
    
    R_CheckUserInterrupt();
    
    it++;
    gamma0 = gamma;
    U.zeros(p2); 
    I.zeros(p2, p2);
    
    ll = getUI_noX(U, I, gamma, time, t, delta, W, S, knots, constantVE);
    
    flag = pinv(I_inv, I);
    
    if (!flag) {
      stop("Newton-Raphson terminated due to singular information matrix");
    }
    
    temp = I_inv*U;
    gamma += temp;
    dif = norm(gamma-gamma0, 2);
  }
  // end: Newton-Raphson
  
  // compute robust variance estimator
  
  ll = getUI_noX(U, I, gamma, time, t, delta, W, S, knots, constantVE);
  w = getw_noX(gamma, time, t, delta, W, S, knots, constantVE);
  covar = I_inv*w.t()*w*I_inv;
  
  List to_return(3); 
  // to_return: ll, covariance matrix
  to_return[0] = ll;
  to_return[1] = covar;
  to_return[2] = gamma;
  return to_return;
}
