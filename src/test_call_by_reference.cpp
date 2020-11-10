#include <Rcpp.h>
using namespace Rcpp;

//' just testing some stuff
//'
//' blab
//' @export
// [[Rcpp::export]]
List test_call_by_reference(NumericMatrix x) {
  int ncol = x.cols();
  int nrow = x.rows();
  x(_, ncol + 1) = seq(1, nrow);

  List ret = List::create();
  ret["boing"] = 1;

  return(ret);
}


