// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace std;
using namespace arma;


// [[Rcpp::export]]
cube mixNormScore (mat lambda, mat draws, cube SigInv, mat mean){
  
  int dim = lambda.n_rows;
  int mix = lambda.n_cols;
  int n = draws.n_rows;
  int p = draws.n_cols;
  
  cube score(dim, mix, n, fill::zeros);
  
  for(int m = 0; m < mix; ++m){
    for(int i = 0; i < n; ++i){
      score.slice(i)(span(0, p-1), m) = SigInv.slice(m) * (draws.row(i).t() - mean.col(m));
      
      mat sigComponent = SigInv.slice(m) * (draws.row(i).t() - mean.col(m)) * (draws.row(i) - mean.col(m).t()) * SigInv.slice(m);  
      mat scoreSig = -SigInv.slice(m) + 0.5 * diagmat(SigInv.slice(m)) + sigComponent - 0.5 * diagmat(sigComponent);
      
      for(int j = 0; j < p; ++j){
        for(int k = 0; k <= j; ++k){
          vec sublambda;
          if(k == 0){
            sublambda = {4, 8, 12, 16};
          } else if(k == 1){
            sublambda = {9, 13, 17};
          } else if(k == 2){
            sublambda = {14, 18};
          } else {
            sublambda = 19;
          }
          for(int l = 0; l < sublambda.n_elem; ++l){
            score(p * (j+1) + k, m, i) += 2 * scoreSig(j, k + l) * lambda(sublambda(l), m);
          }
        }
      }
    }
  }
  return score;
}
