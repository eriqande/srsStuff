#include <Rcpp.h>
#include <iostream>
#include <sstream>
#include <fstream>
using namespace Rcpp;

//' Read the Y dumpfile
//'
//' Ultimately I will want to compute pairwise pop-specific Fst's
//'
//' @param file path to the dumpfile
//' @param num_shorts the number of unsigned shorts that were written to the dumpfile
//' @param L the number of loci expected here
//' @param num_pops the number of populations expected.
//' @return This passes back a list of information. more specifics later
//'
//'
//' @export
// [[Rcpp::export]]
List pairwise_fst_from_dumpfile(
    CharacterVector file,
    int num_shorts,
    int L,
    int num_pops
  ) {

  IntegerMatrix Y0(L, num_pops);
  IntegerMatrix Y1(L, num_pops);

  int i,j,k;
  std::string fname = Rcpp::as<std::string>(file);
  std::ifstream dfile(fname.c_str(), std::ios::binary);
  if(!dfile) {
    Rcpp::stop("Failed to open the dumpfile!");
  }
  double denom;

  unsigned short tmpshort;

  // read everything into two matrices

  for(i=0;i<L;i++) {
    for(j=0;j<num_pops;j++) {
      dfile.read((char *) &tmpshort, sizeof(unsigned short));
      Y0(i, j) = tmpshort;
      dfile.read((char *) &tmpshort, sizeof(unsigned short));
      Y1(i, j) = tmpshort;
    }
    //Rcpp::Rcout << "i: "<< i << "   short_vec[i]: " << short_vec[i] << "\n";
  }

  // OK, now, I guess we go ahead and compute the pop-specific z-bars
  IntegerVector popIdx = seq(1, num_pops);
  NumericVector popL(num_pops);
  NumericVector popZsum(num_pops);
  for(i=0;i<L;i++) {
    for(j=0;j<num_pops;j++) {

      // Hey! Maybe I should weight popZsum differently.  i.e., add up the
      // numerator and the denominator separately, and then divide at the end.
      if(Y0(i,j) + Y1(i,j) > 1) {
        popL(j) += 1.0;
        popZsum(j) += (Y0(i,j) * (Y0(i,j) - 1.0) + Y1(i,j) * (Y1(i,j) - 1.0)) /
          ((Y0(i,j) + Y1(i,j)) * (Y0(i,j) + Y1(i,j) - 1.0));
      }
    }
  }

  // Now, we want to compute the z the "ancestor" of every possible pair
  // of populations.
  IntegerVector pop1(num_pops * (num_pops - 1) / 2);
  IntegerVector pop2(num_pops * (num_pops - 1) / 2);
  NumericVector pairL(num_pops * (num_pops - 1) / 2);
  NumericVector pairSum(num_pops * (num_pops - 1) / 2);
  int cnt = 0;
  for(j=0;j<num_pops;j++) {
    for(k=j+1;k<num_pops;k++) {
      pop1[cnt] = j + 1;
      pop2[cnt] = k + 1;
      cnt++;
    }
  }
  for(i=0;i<L;i++) {
    cnt = 0;
    for(j=0;j<num_pops-1;j++) {
      for(k=j+1;k<num_pops;k++) {
        denom =  ( (Y0(i,j) + Y1(i,j)) * (Y0(i,k) + Y1(i,k)) );
        if(denom > 0.00001) {
          pairL(cnt) += 1.0;
          pairSum(cnt) += ( (Y0(i,j) * Y0(i,k)) + (Y1(i,j) * Y1(i,k))) / denom;
        }
        cnt++;
      }
    }
  }




  DataFrame popz = DataFrame::create(
    _["idx"] = popIdx,
    _["L_used"] = popL,
    _["z_sum"] = popZsum
  );

  DataFrame pairz = DataFrame::create(
    _["pop1"] = pop1,
    _["pop2"] = pop2,
    _["pairL"] = pairL,
    _["pairSum"] = pairSum
  );


  List ret;

  ret["Y0"] = Y0;
  ret["Y1"] = Y1;
  ret["popZ"] = popz;
  ret["pairZ"] = pairz;

  return(ret);
}
