#include <Rcpp.h>
#include <iostream>
#include <sstream>
#include <fstream>
using namespace Rcpp;

//' Compute relative ibd measures from single read sampling
//'
//' This implements the calculations from the
//' "Unified Characterization of Population Structure
//' and Relatedness" framework of Weir and Goudet (2017), but it does it from
//' the allele depth information in a VCF file, by sampling a single
//' read from each individual.
//'
//' This method only works on biallelic markers (but it can work on biallelic indels
//' if you have them).  You must do all the filtering you want to do on your VCF file
//' (using bcftools, for example) and then you have to extract the read depths from that
//' file and put it into a file on your hard drive somewhere.  That file is the
//' input to this function.
//'
//' If you had all your data in a VCF file called raw.vcf, here is how you would
//' process it:
//'
//' 1. Filter it.  Something like this to get biallelic markers typed at over 50\%
//' of individuals (i.e., having at least one read at over 50\% of individuals), and
//' with a minor allele frequency of > 0.05:
//'
//' \code{bcftools view -m 2 -M 2 --min-af 0.05 --max-af 0.95 -i 'F_MISSING < 0.5' raw.vcf > filtered.vcf}
//'
//' 2. Then extract the allele depths out of that. The will look like 0,1 or 2,1 or,
//' if they are missing, they will just be a dot, ".".   So, each row has 2 + N whitespace
//' separated values
//' on it, (where N is the number of individuals). The first two words on each row
//' are Chrom and Pos.  This assumes that the allele depths
//' are comma-separated, and that sites missing reads are denoted by a period:
//'
//' \code{bcftools query -f '\%CHROM\t\%POS[\t\%AD]\n' filtered.vcf > allele_depths.txt}
//'
//' 3. You might want to check that allele_depths has the right number of "words" in it.  So, you
//' would do \code{wc allele_depths.txt}, and confirm that the second number in the output
//' is equal to \eqn{(N + 2) * L}, where \eqn{N} is the number of individuals (samples)
//' in the data set, filtered.vcf, and \eqn{L} is the number of markers in filtered.vcf.
//' In other words, the file should have L rows, each one with N + 2 columns: 2 for Chrom and
//' Pos and then 1 column for each of the N individuals.
//'
//' 4. Then, you call this function, giving it the path to the file allele_depths.txt,
//' and you also pass in a list of vectors that hold information about the
//' hierarchical levels you want to use to be able to compute relative measures of
//' identity by descent.
//'
//' @param file path to the file that holds the Chrom, Pos, [AD] file (i.e. allele_depths.txt
//' in the description above).
//' @param sample_names A character vector of the names of the sample, in the
//' order they appear in the VCF file.
//' @param freq_thresh loci with the frequency of either allele estimated (by the fraction
//' of sampled single reads of each type) less than freq_thresh will not be used.
//'
//' @return This passes back a list of information about the various allele
//' sharing statistics (and, eventually bootstrap information) that can be
//' used by a higher-level function to compute the F-statistics, etc.
//'
//'
//' @export
// [[Rcpp::export]]
List srs_identity(CharacterVector file,
               List groups,
               CharacterVector sample_names,
               int BootReps,
               double freq_thresh = 0.0) {
  int comma, i, j, r;
  std::string tempstr;
  std::string line;
  std::string Chrom;
  std::string Pos;
  std::string word;
  std::string fname = Rcpp::as<std::string>(file);
  int N = 0; // counts the individuals
  int L = -1; // counts the markers
  int rdr, rda; // read depth ref / read depth alt
  IntegerVector alleles(sample_names.size());  // one for each sample, -1 = no read; 0 = ref; 1 =
  IntegerMatrix allele_mat(BootReps, sample_names.size());
  int totAlt, totNotMissing; // for computing allele frequency from single reads
  int NumLociAllMissing = 0;  // keep track of how many loci have no reads at any individuals
  int NumLociMonomorphic = 0; // keep track of how many loci were monomorphic after read sampling. These can't be used
  int NumLociPolymorph = 0;  // how many markers were not missing at everyone and were polymorphic after sampling
  int NumLociSuitableMAF = 0; // how many polymorphic markers had high enough MAF?
  int NumLociTooLowFreq = 0; // How many polymorphic sites were skipped cuz of low MAF?

  // Stuff for keeping track of which individuals have reads of both alleles
  IntegerVector hetIdxs(sample_names.size());
  IntegerVector hetFlags(sample_names.size());
  int numHets;


  // here is a bunch of stuff for dealing with different levels of identity
  int num_id_levs = groups.length();
  if(num_id_levs < 2) Rcpp::stop("Hey, you have to have at least two levels of groups (indivs + pops)");
  IntegerVector id_lev_lengths(num_id_levs);
  IntegerVector lev_vec;

  // and this is for counting up the number of different
  // alleles in different bins in different hierarchical levels

  // Here, declare matrices to send some information back
  NumericMatrix retCov(sample_names.size(), sample_names.size());
  NumericMatrix retIBS(sample_names.size(), sample_names.size());
  NumericMatrix retM(sample_names.size(), sample_names.size());  // number non-missing + polymorphic loci
  // initialize those to 0.0s.  Rcpp11 has a nicer constructor, but not here, apparently.
  for(i=0;i<sample_names.size();i++) {
    for(j=0;j<sample_names.size();j++) {
      retCov(i,j) = 0.0;
      retIBS(i,j) = 0.0;
      retM(i,j) = 0.0;
    }
  }


  // compile information about how many elements are in each level of
  // the hierarchy of pops/regions/etc. for the identity stuff. And also
  // allocate memory for counting up 0 and 1 alleles in each of the
  // bins at each of the hierarhical levels
  std::vector< IntegerMatrix > bin_counts(num_id_levs) ;
  for(i=0;i<num_id_levs;i++) {
    lev_vec = as<IntegerVector>(groups[i]);
    id_lev_lengths(i) = lev_vec.length();
    bin_counts[i] = IntegerMatrix(id_lev_lengths(i), 2);
  }

  // Here are things for the fixed counts of 0's and 1's in populations, and also
  // the Randomized counts from Het Sampling.
  IntegerVector Fixed0s(id_lev_lengths(1));
  IntegerVector Fixed1s(id_lev_lengths(1));
  IntegerMatrix Rando0s(BootReps, id_lev_lengths(1));
  IntegerMatrix Rando1s(BootReps, id_lev_lengths(1));


  // open up a stream to read from the file
  std::ifstream infile (fname.c_str());

  ////////// this code goes line through line in the file and samples //////////
  ////////// a single read from each indiv                            //////////
  while (std::getline(infile, line)) {
    std::stringstream ss(line);
    L++; // increment the locus counter

    ss >> Chrom;  // eat the Chrom name
    ss >> Pos; // eat the position of the variant

    // now, cycle over remaning columns and read them
    N = -1;
    totAlt = 0;
    totNotMissing = 0;
    numHets = 0;
    Rcpp::Rcout << "Entering Row\n";
    while(ss >> word) {
      N++;  // for keeping track of which individual we are on.

      if(N + 1 > sample_names.size()) {
        Rcpp::Rcerr << "At Chrom:Pos " << Chrom << ":" << Pos << " have read " << N + 1 << " samples. Which is greater than  " <<  sample_names.size() << " from sample_names\n";
        Rcpp::stop("Insufficient number of sample names.  Bailing out!");
      }

      if(word == ".") {  // no reads from the individual at this site
        rdr = 0L;
        rda = 0L;
      } else {
        comma = word.find(",");
        rdr = std::stoi(word.substr(0, comma)) ;
        rda = std::stoi(word.substr(comma + 1));
      }
      // at this point, rda and rdr are the read depths for individual N at marker L
      // (both of which are base-0 indexed.)

      // now we sample one of the reads
      hetFlags(N) = 0L;
      if(rdr == 0 && rda == 0) {
        alleles[N] = -1;
        allele_mat(0, N) = -1;
      } else if(rdr > 0 && rda == 0) {
        alleles[N] = 0;
        allele_mat(0, N) = 0;
      } else if(rdr == 0 && rda > 0) {
        alleles[N] = 1;
        allele_mat(0, N) = 1;
      } else {
        alleles[N] = runif(1, 0, 1)[0] < ((double) rda / (rdr + rda)); // sample alternate allele with prob rda/(rda+rdr)
        allele_mat(_, N) = runif(BootReps, 0, 1) < ((double) rda / (rdr + rda)); // sample alternate allele with prob rda/(rda+rdr)
        hetIdxs(numHets) = N;
        hetFlags(numHets++) = 1L;
      }

      if(alleles[N] != -1) {
        totAlt += alleles[N];
        totNotMissing++;
      }
      //Rcpp::Rcout << L << "  " << N << "  " << word << "  " << rdr << "  " << rda << "  " << comma << "      " << wdr << "  " << wda << "     " <<  alleles[N] << "\n";

    } // close loop that reads columns from each line
    N++;  // add one at the end.

    // throw an error if we didn't get the number expected.
    if(N != sample_names.size()) {
      Rcpp::Rcerr << "At Chrom:Pos " << Chrom << ":" << Pos << " just read " << N << " samples. Expected " <<  sample_names.size() << " from sample_names\n";
      Rcpp::stop("Didn't find the right number of samples.");
    }

    // Now, only the do the next parts for markers that are not all missing and are polymorphic
    if(totNotMissing == 0) {
      NumLociAllMissing++;
    } else if(totAlt == 0 || totAlt == totNotMissing) {
      NumLociMonomorphic++;
    } else {
      NumLociPolymorph++;

      ///////////////////  HERE IS WHERE THE ACTUAL CALCULATIONS FOR EACH LOCUS OCCUR /////////

      // First we compute the Fixed and Rando counts for the populations (first level above individuals)
      // zero out to accumulate sums
      //Rcpp::Rcout << "Zeroing stuff\n";
      for(i=0;i<id_lev_lengths(1);i++) {
        Fixed0s(i) = 0L;
        Fixed1s(i) = 0L;
        for(r=0;r<BootReps;r++) {
          Rando0s(r, i) = 0L;
          Rando1s(r, i) = 0L;
        }
      }
      //Rcpp::Rcout << "Done Zeroing stuff\n";
      // Now, cycle over the individuals, and add things up according to which population each is in
      lev_vec = as<IntegerVector>(groups[0]); // this gives us the idx of the population of each indvidual
      for(i=0;i<N;i++) {
        Fixed0s(lev_vec(i)) += (hetFlags(i) == 0) && allele_mat(0, i) == 0;
        Fixed1s(lev_vec(i)) += (hetFlags(i) == 0) && allele_mat(0, i) == 1;
      }
      for(i=0;i<numHets;i++) {
        for(r=0;r<BootReps;r++) {
        Rando0s(r, lev_vec(hetIdxs(i))) += allele_mat(r, hetIdxs(i)) == 0;
        Rando1s(r, lev_vec(hetIdxs(i))) += allele_mat(r, hetIdxs(i)) == 1;
        }
      }

      // COOL!  For all of the identity calculations above the level of individuals
      // we can now use Fixed0s, Fixed1s, and Rando0s and Rando1s to get the quantities
      // we need at each level and for each resample/bootstrap replicate.



      for(j=0;j<num_id_levs;j++) {
        for(i=0;i<id_lev_lengths(j);i++) {
            bin_counts[j](i, 0) = 0L;
            bin_counts[j](i, 1) = 0L;
        }
      }

      Rcpp::Rcout << "alleles: " << alleles << "\n";

      Rcpp::Rcout << "Filling stuff\n";
      // Then put the values in there
      for(j=0;j<num_id_levs - 1;j++) { // j cycles over the levels in the hierarchy
        for(i=0;i<id_lev_lengths(j);i++) {  // i cycles over the elements in that hierarchy
          if(j==0) {
            if(alleles[i] != -1) {
              bin_counts[j](i, alleles[i]) = 1L;
              lev_vec = as<IntegerVector>(groups[j]);
              bin_counts[j+1](lev_vec(i), alleles[i])++;
            }
          } else {
            lev_vec = as<IntegerVector>(groups[j]);
            bin_counts[j+1](lev_vec(i), 0) += bin_counts[j](i, 0);
            bin_counts[j+1](lev_vec(i), 1) += bin_counts[j](i, 1);
          }
        }
      }
      //return(List::create(bin_counts[0], bin_counts[1], bin_counts[2], bin_counts[3]));

     /////////////////// END OF ACTUAL CALCULATIONS FOR EACH LOCUS ///////////////////////////
    } // close if marker is polymorphic and not all missing

    return(List::create(allele_mat, numHets, hetIdxs, Fixed0s, Fixed1s, Rando0s, Rando1s));
  } // close loop over markers (rows in file)


  // store to return

  NumericVector retNums(6);
  retNums[0] = L + 1;
  retNums[1] = NumLociPolymorph;
  retNums[2] = NumLociMonomorphic;
  retNums[3] = NumLociAllMissing;
  retNums[4] = NumLociTooLowFreq;
  retNums[5] = NumLociSuitableMAF;
  retNums.names() = CharacterVector::create("all_sites",
                "polyorphic_sites",
                "monomorphic_sites",
                "sites_missing_reads_from_everyone",
                "polymorphic_sites_skipped_for_low_MAF",
                "num_polymorph_sites_passing_MAF");

  List ret;
  ret["IBS"] = retIBS;
  ret["Cov"] = retCov;
  ret["M"] = retM;
  ret["sample_names"] = sample_names;
  ret["freq_thresh"] = freq_thresh;
  ret["site_counts"] = retNums;

  return(ret);
}

