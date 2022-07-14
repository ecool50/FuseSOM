#include <Rcpp.h>
using namespace Rcpp;

/* 
this function was obtained from https://github.com/cran/cstab with some minor modifications
*/

// [[Rcpp::export]]
std::vector<bool> equal(std::vector<int> x) {
  int n_long = x.size(), n_short = x.size()-1;
  std::vector<bool> res((n_long*n_short)/2);
  int ind = 0;
  for(int i = 0; i < n_short; ++i){
    for(int j = (i+1); j < n_long; ++j){
      res[ind] = x[i] == x[j];
      //std::cout << x[i] << '_' << x[j] << '\n';
      ind++;
    }
  }
  return res;
}
