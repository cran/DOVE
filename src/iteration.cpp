#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <ctime>
#include "BS.h"
#include "compare.h"

struct scoreData {
  arma::mat score;
  arma::mat dscore;
};



//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;



scoreData cov1Comp(arma::mat& exbt, 
                   int np, 
                   int m, 
                   int n, 
                   arma::mat& X, 
                   arma::mat& DZBBy) {

  arma::mat dscore1(np+m, np+m);
  arma::mat score1(n,np+m);

  arma::mat tmp, tmp2;
  arma::vec tt = sum(exbt,1);

  score1.cols(0,np-1) = DZBBy.cols(0,np-1) - tt % X.each_col();

  int k = 0;
  while (k < np) {
    // {n x m}
    tmp = X.col(k) % exbt.each_col();

    for (int j = k; j < np; ++j) {
      dscore1(k,j) = -accu(X.col(j) % tmp.each_col());
      dscore1(j,k) = dscore1(k,j);
    }
    for (int j = np; j < (np+m); ++j) {
      dscore1(k,j) = -sum(tmp.col(j-np));
      dscore1(j,k) = dscore1(k,j);
    }
    k++;
  }

  scoreData to_return;

  to_return.score = score1;
  to_return.dscore = dscore1;

  return to_return;

}

scoreData em1Comp(arma::mat& exbt, int np, int m, int n, arma::mat& DZBBy) {

  arma::mat dscore1(m,m);
  arma::mat score1(n,m);

  score1 = DZBBy.cols(np, np+m-1) - exbt;

  dscore1.diag() = -sum(exbt,0);

  scoreData to_return;

  to_return.score = score1;
  to_return.dscore = dscore1;

  return to_return;

}

scoreData score1(arma::mat& X, arma::vec& beta, arma::mat& dt, arma::mat& DZBBy) {

  int np = X.n_cols;
  int m = beta.n_elem - np;
  int n = dt.n_rows;

  // {n}
  arma::vec zbeta = X * beta.subvec(0,np-1);

  // {n x m}
  arma::mat EXbeta(n,m);
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < n; ++i) {
      EXbeta(i,j) = exp(zbeta(i) + beta(j+np));
    }
  }

  // {n x m}
  arma::mat exbt = EXbeta % dt;

  scoreData resc = cov1Comp(exbt, np, m, n, X, DZBBy);

  scoreData resm = em1Comp(exbt, np, m, n, DZBBy);

  arma::mat A(n,np), B(np+m,np), C(np,m);
  arma::mat S1 = resc.score, S2 = resm.score;
  arma::mat DS1 = resc.dscore, DS2 = resm.dscore;

  scoreData to_return;

  to_return.score = S1 + join_horiz(A, S2);
  to_return.dscore = DS1 + join_horiz(B, join_vert(C, DS2));

  return to_return;

}

scoreData score1_NoX(arma::vec& beta, arma::mat& dt, arma::mat& DZBBy) {

  int m = beta.n_elem;
  int n = dt.n_rows;

  // {n x m}
  arma::mat EXbeta(n,m);
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < n; ++i) {
      EXbeta(i,j) = exp(beta(j));
    }
  }

  // {n x m}
  arma::mat exbt = EXbeta % dt;

  scoreData resm = em1Comp(exbt, 0, m, n, DZBBy);

  return resm;

}

scoreData cov2Comp(arma::mat& X, 
                   arma::mat& iexb, 
                   arma::imat& loc, 
                   arma::mat& prepStep, 
                   int np, 
                   int m) {

  // iexb is {nVac x nEventsVac}

  int nEventsVac = iexb.n_cols;
  int nVac = iexb.n_rows;

  arma::mat dscore2(np+m, np+m);

  // {nEventsVac}
  arma::vec denom = sum(iexb,0).t();

  // {nEventsVac x np}
  arma::mat uuu = iexb.t() * X;
  uuu = uuu.each_col() / denom;

  // {nEventsVac x m}
  arma::mat ttt(nEventsVac,m);
  arma::vec tcol;
  for (arma::sword i = 0; i < m; ++i) {
    for (int j = 0; j < nEventsVac; ++j) {
      tcol = iexb.col(j);
      ttt(j,i) = accu(tcol.elem(find(loc.col(j) == (np+i+1)))) / denom(j);
    }
  }  

  arma::mat nomjk;
  arma::mat tempk(nEventsVac, nVac);

  arma::mat score2 = prepStep.cols(0,np-1) - uuu;
  arma::mat dscore2_npnp = uuu.t() * uuu;
  arma::mat dscore2_mnp = ttt.t() * uuu;
  double dnomjk;

  for (int k = 0; k < np; ++k) {

    // d exp(Z beta + gamma) / d beta_k
    // {nVac x nEventsVac}
    tempk = iexb.each_col() % X.col(k);

    nomjk = tempk.t() * X;

    for (int j = k; j < np; ++j) {
      dscore2(k,j) = -sum(nomjk.col(j) / denom) + dscore2_npnp(k,j);
      dscore2(j,k) = dscore2(k,j);
    }

    for (arma::sword j = 0; j < m; ++j) {

      dnomjk = 0.0;
      for (int i = 0; i < nEventsVac; ++i) {
        tcol = tempk.col(i);
        dnomjk = dnomjk + accu(tcol.elem(find(loc.col(i) == (np+j+1)))) / denom(i);
      }  

      dscore2(k,j+np) = -dnomjk + dscore2_mnp(j,k);
      dscore2(j+np,k) = dscore2(k,j+np);
    }
  }

  scoreData to_return;

  to_return.score = score2;
  to_return.dscore = dscore2;

  return to_return;
}

scoreData em2Comp(arma::mat& iexb, 
             arma::mat& prepStep,  
             arma::imat& loc,  
             int m,
             int np) {

  int j;

  int nEventsVac = iexb.n_cols;

  // {nEventsVac}
  arma::vec denom = sum(iexb, 0).t();

  // {nEventsVac x m}
  arma::mat ttt(nEventsVac,m);
  arma::vec tcol;
  for (arma::sword i = 0; i < m; ++i) {
    for (int j = 0; j < nEventsVac; ++j) {
      tcol = iexb.col(j);
      ttt(j,i) = accu(tcol.elem(find(loc.col(j) == (np+i+1)))) / denom(j);
    }
  }

  arma::mat score2 = prepStep.cols(np, np+m-1) - ttt;
  arma::mat dscore2 = ttt.t() * ttt;
  dscore2.diag() = dscore2.diag() - sum(ttt,0).t();

  scoreData to_return;

  to_return.score = score2;
  to_return.dscore = dscore2;

  return to_return;

}

scoreData score2(arma::mat& X,  
                 arma::vec& beta,  
                 arma::imat& indVac,  
                 arma::imat& loc, 
                 arma::mat& prepStep) {

  // number of covariates
  int np = X.n_cols;
  // number of intervals
  int m = beta.n_elem - np;

  // number of vaccinated
  int nVac = indVac.n_rows;
  // number of events for vaccinated
  int nEventsVac = indVac.n_cols;

  // {nVac} X beta
  arma::vec zbeta = X * beta.subvec(0,np-1);

  // assign appropriate beta to each bin
  // {nVac x nEventsVac}
  arma::mat tmp(nVac, nEventsVac);
  arma::uvec ids;
  for (int i = 0; i < m; ++i) {
    ids = find(loc == (np+i+1));
    tmp.elem(ids).fill(beta(np+i));
  }

  arma::mat EXbeta = exp(zbeta + tmp.each_col());

  // {nVac x nEventsVac}
  arma::mat iexb = indVac % EXbeta;

  scoreData resc = cov2Comp(X, iexb, loc, prepStep, np, m);

  scoreData resm = em2Comp(iexb, prepStep, loc, m, np);

  scoreData to_return;

  arma::mat A(np+m,np), B(np,m);
  arma::mat S1 = resc.score, S2 = resm.score;
  arma::mat DS1 = resc.dscore, DS2 = resm.dscore;

  to_return.score = join_horiz(S1, S2);
  arma::mat tt = join_vert(B,DS2);
  to_return.dscore = DS1 + join_horiz(A, tt);

  return to_return;
}

scoreData score2_NoX(arma::vec& beta,  
                     arma::imat& indVac,  
                     arma::imat& loc, 
                     arma::mat& prepStep) {

  int m = beta.n_elem;

  int nVac = indVac.n_rows;
  int nEventsVac = indVac.n_cols;

  // {nVac x nEventsVac}
  arma::mat tmp(nVac, nEventsVac);
  arma::uvec ids;
  for (int i = 0; i < m; ++i) {
    ids = find(loc == (i+1));
    tmp.elem(ids).fill(beta(i));
  }

  arma::mat EXbeta = exp(tmp);

  // {nVac x nEventsVac}
  arma::mat iexb = indVac % EXbeta;

  scoreData resm = em2Comp(iexb, prepStep, loc, m, 0);

  return resm;
}


// [[Rcpp::export]]
List iteration(arma::mat& dt, arma::mat& DZBBy, arma::mat& X, 
               arma::imat& indVac, arma::imat& loc, arma::mat& score2Prep,
               int m, arma::uvec& vac) {

  int np = X.n_cols, success;
  
  arma::vec oldbeta(np+m), newbeta(np+m), inv2(np+m);

  double err;

  oldbeta.zeros();

  arma::mat Xvac = X.rows(find(vac==1));
  arma::mat DS1, DS2, inv;
  arma::rowvec sumsc1, sumsc2;

  err = 1.0;
  success = 0;

  scoreData s1, s2;

  for (int i=1; i <= 500; ++i) {
    R_CheckUserInterrupt();
    s1 = score1(X, oldbeta, dt, DZBBy);

    R_CheckUserInterrupt();
    s2 = score2(Xvac, oldbeta, indVac, loc, score2Prep);

    sumsc1 = sum(s1.score,0);
    DS1 = s1.dscore;
    
    sumsc2 = sum(s2.score,0);
    DS2 = s2.dscore;

    inv = pinv(DS1 + DS2);
    inv2 = inv * (sumsc1 + sumsc2).t();

    newbeta = oldbeta - inv2;

    err = sum(abs(newbeta-oldbeta));
    if (err <= 1e-4) {
      success = i;
      break;
    }

    oldbeta = newbeta;

  }

  if (success == 0) throw std::range_error("method did not converge after 500 iterations");

  List to_return(6);

  to_return[0] = success;
  to_return[1] = oldbeta;
  to_return[2] = s1.score;
  to_return[3] = s1.dscore;
  to_return[4] = s2.score;
  to_return[5] = s2.dscore;


  return to_return;
}


// [[Rcpp::export]]
List iteration_NoX(arma::mat& dt, 
                   arma::mat& DZBBy,
                   arma::imat& indVac, 
                   arma::imat& loc, 
                   arma::mat& score2Prep,
                   int m) {

  int success;
  
  arma::vec oldbeta(m), newbeta(m), inv2(m);

  double err;

  oldbeta.zeros();

  arma::mat DS1, DS2, inv;
  arma::rowvec sumsc1, sumsc2;

  err = 1.0;
  success = 0;

  scoreData s1, s2;
  for (int i=1; i <= 500; ++i) {

    R_CheckUserInterrupt();
    s1 = score1_NoX(oldbeta, dt, DZBBy);

    R_CheckUserInterrupt();
    s2 = score2_NoX(oldbeta, indVac, loc, score2Prep);

    sumsc1 = sum(s1.score,0);
    DS1 = s1.dscore;
    
    sumsc2 = sum(s2.score,0);
    DS2 = s2.dscore;

    inv = pinv(DS1 + DS2);
    inv2 = inv * (sumsc1 + sumsc2).t();

    newbeta = oldbeta - inv2;

    err = sum(abs(newbeta-oldbeta));
    if (err <= 1e-4) {
      success = i;
      break;
    }
    oldbeta = newbeta;

  }

  if (success == 0) throw std::range_error("method did not converge after 500 iterations");

  List to_return(6);

  to_return[0] = success;
  to_return[1] = oldbeta;
  to_return[2] = s1.score;
  to_return[3] = s1.dscore;
  to_return[4] = s2.score;
  to_return[5] = s2.dscore;


  return to_return;
}
