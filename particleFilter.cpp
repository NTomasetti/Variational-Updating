// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace std;


// [[Rcpp::export]]
Rcpp::List particleFilter(vec y, vec draw, vec lambda, bool initial, int S = 50, bool prior = true){
  
  // Get variance matrix components
  vec mean = {lambda(0), lambda(1)};
  mat L(3, 3);
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j <= i; ++j){
      L(i, j) = lambda(3 + 3*i + j);
    }
  }
  mat Sigma = L * L.t();
  mat SigSubset = Sigma.submat(0, 0, 1, 1);
  mat SigInv = SigSubset.i();
  
  // Prior distribution
  double logprior = -0.5 * as_scalar((mean - draw).t() * SigInv * (mean - draw));
  
  // State prior
  
  double stateMean, stateSd;
  if(initial){
    stateMean = 0;
    stateSd = sqrt(draw(1) / (1 - pow(draw(0), 2)));
  } else {
    // Conditional distribution of x_T
    vec covar = {Sigma(0, 2), Sigma(1, 2)};
    stateMean = lambda(2) + as_scalar(covar.t() * SigInv * (draw - mean));
    double stateVar = Sigma(2, 2) - as_scalar(covar.t() * SigInv * covar);
     // Linear / Normal relationship for x_T+1
    stateMean = draw(0) * stateMean;
    stateVar = stateVar * pow(draw(0), 2) + draw(1);
    stateSd = sqrt(stateVar);
  }

  double logdens = 0;
  int T = y.n_elem;
  
  vec v(S), vResample(S), w(S), cdf(S);
  for(int t = 0; t < y.n_elem; ++t){
    
    if(t == 0){
      v = stateMean + stateSd * randn<vec>(S);
    } else {
      v = vResample * draw(0) + sqrt(draw(1)) * randn<vec>(S);
    }
 
    for(int i = 0; i < S; ++i){
      w(i) = 1.0 / sqrt(2 * 3.14159 * exp(v(i))) * exp(- pow(y(t), 2) / (2 * exp(v(i))));
    }
    logdens += log(sum(w) / S);
    
    w = w / sum(w);
    cdf = cumsum(w);
    
    for(int i = 0; i < S; ++i){
      double u = randu<double>();
      for(int k = 0; k < S; ++k){
        if(u < cdf(k)){
          vResample(i) = v(k);
          break;
        }
      }
    }
  }
  if(prior){
    return Rcpp::List::create(Rcpp::Named("dens") = logprior + logdens,
                              Rcpp::Named("state") = v(0));
  } else {
    return Rcpp::List::create(Rcpp::Named("dens") = logdens,
                              Rcpp::Named("state") = v(0));
  }
 
}
