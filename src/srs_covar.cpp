#include <Rcpp.h>
#include <iostream>
#include <sstream>
#include <fstream>
using namespace Rcpp;

//' Compute covariance between individuals based on sampling a single read at each site
//'
//' This implements the covariance calculation method in ANGSD described at
//' [http://www.popgen.dk/angsd/index.php/PCA_MDS](http://www.popgen.dk/angsd/index.php/PCA_MDS),
//' but it does it from
//' the allele depth information in a VCF file.  No need to go back to your BAMS, dude!
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
//'     ```
//'     bcftools view \
//'       -m 2 -M 2 \
//'       --min-af 0.05 --max-af 0.95 \
//'       -i 'F_MISSING < 0.5' \
//'       raw.vcf > filtered.vcf
//'     ```
//'
//' 2. Then extract the allele depths out of that. The will look like 0,1 or 2,1 or,
//' if they are missing, they will just be a dot, ".".   So, each row has 2 + N whitespace
//' separated values
//' on it, (where N is the number of individuals). The first two words on each row
//' are Chrom and Pos.  This assumes that the allele depths
//' are comma-separated, and that sites missing reads are denoted by a period:
//'     ```
//'     bcftools query \
//'       -f '%CHROM\t%POS[\t%AD]\n' \
//'       filtered.vcf > allele_depths.txt
//'     ```
//'
//' 3. You might want to check that allele_depths has the right number of "words" in it.  So, you
//' would do `wc allele_depths.txt`, and confirm that the second number in the output
//' is equal to \eqn{(N + 2) * L}, where \eqn{N} is the number of individuals (samples)
//' in the data set, `filtered.vcf`, and \eqn{L} is the number of markers in `filtered.vcf`.
//' In other words, the file should have \eqn{L} rows, each one with \eqn{N + 2} columns: 2 for Chrom and
//' Pos and then 1 column for each of the \eqn{N} individuals.
//'
//' 4. Then, you call this function, giving it the path to the file allele_depths.txt,
//' and you also pass in a vector of names for the individuals.
//' @param file path to the file that holds the Chrom, Pos, and allele depths file (i.e. `allele_depths.txt`
//' in the description above).
//' @param sample_names A character vector of the names of the sample, in the
//' order they appear in the VCF file.
//' @param freq_thresh loci with the frequency of either allele estimated (by the fraction
//' of sampled single reads of each type) less than freq_thresh will not be used.
//'
//' @return This passes back a list that includes an \eqn{N x N} covariance matrix ($Cov); a martrix of
//' proportion of sampled reads identical by state between individuals ($IBS); a matrix of
//' number of sites having at least one read in both individuals of the pair ($M); the sample names,
//' the frequency threshold used in this function; and a report about the total number of loci
//' investigated and used.
//'
//'
//' @export
// [[Rcpp::export]]
List srs_covar(CharacterVector file,
                                             CharacterVector sample_names,
                                             double freq_thresh = 0.0) {
  int comma, i, j;
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
  int totAlt, totNotMissing; // for computing allele frequency from single reads
  int NumLociAllMissing = 0;  // keep track of how many loci have no reads at any individuals
  int NumLociMonomorphic = 0; // keep track of how many loci were monomorphic after read sampling. These can't be used
  int NumLociPolymorph = 0;  // how many markers were not missing at everyone and were polymorphic after sampling
  int NumLociSuitableMAF = 0; // how many polymorphic markers had high enough MAF?
  int NumLociTooLowFreq = 0; // How many polymorphic sites were skipped cuz of low MAF?
  double fa, r, a;  // for alt allele freq, and the number of ref and alt alleles

  // Here, declare matrices to send some information back
  NumericMatrix retCov(sample_names.size(), sample_names.size());
  NumericMatrix retIBS(sample_names.size(), sample_names.size());
  NumericMatrix retM(sample_names.size(), sample_names.size());  // number non-missing + polymorphic loci
  NumericMatrix retMt_S(sample_names.size(), sample_names.size()); // to store the mean proportion of pairs
  // drawn from separate individuals that are IBS.  For the Weir and Goudet calculation

  // initialize those to 0.0s.  Rcpp11 has a nicer constructor, but not here, apparently.
  for(i=0;i<sample_names.size();i++) {
    for(j=0;j<sample_names.size();j++) {
      retCov(i,j) = 0.0;
      retIBS(i,j) = 0.0;
      retM(i,j) = 0.0;
      retMt_S(i,j) = 0.0;
    }
  }

  // open up a stream to read from the file
  std::ifstream infile (fname.c_str());

  while (std::getline(infile, line)) {
    std::stringstream ss(line);
    L++; // increment the locus counter

    ss >> Chrom;  // eat the Chrom name
    ss >> Pos; // eat the position of the variant

    // now, cycle over remaning columns and read them
    N = -1;
    totAlt = 0;
    totNotMissing = 0;
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
      if(rdr == 0 && rda == 0) {
        alleles[N] = -1;
      } else if(rdr > 0 && rda == 0) {
        alleles[N] = 0;
      } else if(rdr == 0 && rda > 0) {
        alleles[N] = 1;
      } else {
        alleles[N] = runif(1, 0, 1)[0] < ((double) rda / (rdr + rda)); // sample alternate allele with prob rda/(rda+rdr)
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

      // compute the alt allele frequency from the sampled reads
      fa = (double)totAlt / totNotMissing;
      a = (double)totAlt;
      r = (double)(totNotMissing - totAlt);

      if(fa > freq_thresh && fa < 1 - freq_thresh) {  // here, impose a MAF cutoff if desired
        NumLociSuitableMAF++;
        for(i=0;i<N;i++) {
          for(j=0;j<=i;j++) {  // just fill out the lower triangle here
            if(alleles[i] >= 0 && alleles[j] >= 0) {
              retCov(i, j) += (alleles[i] - fa) * (alleles[j] - fa) / (fa * (1.0 - fa));
              retIBS(i, j) += (alleles[i] == alleles[j]);
              retM(i,j) += 1.0;
              retMt_S(i,j) += (r * (r - 1) + a * (a - 1)) / (totNotMissing * (totNotMissing - 1));
            }
          }
        }
      } else {
        NumLociTooLowFreq++;
      }
    } // close if marker is polymorphic and not all missing
  } // close loop over markers (rows in file)

  // Now, divide all those entries with the appropriate M, and fill out the
  // upper triangle too, making it symmetric.
  for(i=0;i<N;i++) {
    for(j=0;j<=i;j++) {
      retCov(i, j) = retCov(i, j) / retM(i,j);
      retIBS(i, j) = retIBS(i, j) / retM(i,j);
      retMt_S(i,j) = retMt_S(i,j) / retM(i,j);

      if(i != j) {
        retCov(j, i) = retCov(i, j);
        retIBS(j, i) = retIBS(i, j);
        retM(j, i) = retM(i, j);
        retMt_S(j, i) = retMt_S(i,j);
      }
    }
  }
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
  ret["Mt_S"] = retMt_S;
  ret["sample_names"] = sample_names;
  ret["freq_thresh"] = freq_thresh;
  ret["site_counts"] = retNums;

  return(ret);
}

