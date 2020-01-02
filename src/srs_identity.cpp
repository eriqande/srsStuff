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
//' Note that it would not be too hard to modify this to deal with multiple alleles.
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
//' @param pops_of_indivs a vector of 0-based indices of the populations each of the individuals
//' belongs to. Must be in the order that individuals appear in the VCF file.
//' @param num_pops The number of internal nodes in the tree that are directly above the individuals.
//' @param num_internal_nodes the number of internal nodes in the tree (this includes the populations)
//' @param daughters a list of vectors.  Each one is the 0-based indexes of the daughters
//' of the non-population internal nodes. Here 0 corresponds to the first internal node (the first population).
//' @param sample_names A character vector of the names of the sample, in the
//' order they appear in the VCF file.
//' @return This passes back a list of information about the various allele
//' sharing statistics (and, eventually bootstrap information) that can be
//' used by a higher-level function to compute the F-statistics, etc.
//'
//'
//' @export
// [[Rcpp::export]]
List srs_identity(CharacterVector file,
                  IntegerVector pops_of_indivs,
                  int num_pops,
                  int num_internal_nodes,
                  List daughters,
                  CharacterVector sample_names,
                  int BootReps) {
  int comma, i, j, r, k, d, dp, k1, k2;
  std::string tempstr;
  std::string line;
  std::string Chrom;
  std::string Pos;
  std::string word;
  std::string fname = Rcpp::as<std::string>(file);
  int N = 0; // counts the individuals
  int L = -1; // counts the markers
  int rdr, rda; // read depth ref / read depth alt
  IntegerVector alleles(pops_of_indivs.size());  // one for each sample, -1 = no read; 0 = ref; 1 =
  IntegerMatrix allele_mat(BootReps, pops_of_indivs.size());
  int totAlt, totNotMissing; // for computing allele frequency from single reads
  double denom;

  // Stuff for keeping track of which individuals have reads of both alleles
  IntegerVector hetIdxs(pops_of_indivs.size());   // vector of length num_indivs, but we only use the first numHets elements, where
                                                  // numHets is the number of individuals with reads of both alleles at this locus.
                                                  // It stores the idxs of the individuals that have reads of both alleles at the locus.
  IntegerVector hetFlags(pops_of_indivs.size());  // vector of length num_indivs, each element 0, or, if indiv has reads from each allele it's a 1
  int numHets;

  // Here are things for the fixed counts of 0's and 1's in populations, and also
  // the Randomized counts from Het Sampling.
  IntegerVector Fixed0s(num_pops);
  IntegerVector Fixed1s(num_pops);
  IntegerMatrix Rando0s(BootReps, num_pops);
  IntegerMatrix Rando1s(BootReps, num_pops);

  // Here are the YO and Y1 matrices which will be used for each locus
  // to count up the number of alleles beneath each internal node
  NumericMatrix Y0(BootReps, num_internal_nodes);
  NumericMatrix Y1(BootReps, num_internal_nodes);
  IntegerVector Ytot0(BootReps);  // these next two are for counting or 0s and 1s over all population in each boot rep
  IntegerVector Ytot1(BootReps);  // so as to ignore reps in which there is no variation.

  // This is what we will end up returning.  Note, the first row will be the actual estimate
  // and we will effect that by setting bootn[0] = 1L for that always
  NumericMatrix zsum(BootReps, num_internal_nodes);

  // And here we have a matrix to store the number of loci (L) to divide
  // each zsum by.  (It can vary across nodes/pops and reps, so we keep
  // a whole matrix of these)
  IntegerMatrix zL(BootReps, num_internal_nodes);

  // initialize both of those to 0.0s.  Rcpp11 has a nicer constructor, but not here, apparently.
  for(i=0;i<BootReps;i++) {
    for(j=0;j<num_internal_nodes;j++) {
      zsum(i,j) = 0.0;
      zL(i, j) = 0L;
    }
  }



  // This is to hold the number of times the current locus is selected
  // for use in each of the BootReps bootstrap replicates.
  IntegerVector bootn(BootReps);

  double zbar;  // a temp variable that we will want to use to store things.

  IntegerVector daught;  // to hold the results that come out of the daughters list for each level
  int ndaught; // to store the number of daughters
  int npairs; // to store the number of pairs of daughters

  // open up a stream to read from the file
  std::ifstream infile (fname.c_str());

  ////////// this code goes line by line in the file and samples //////////
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
    //Rcpp::Rcout << "Entering Row\n";
    while(ss >> word) {
      N++;  // for keeping track of which individual we are on.

      if(N + 1 > pops_of_indivs.size()) {
        Rcpp::Rcerr << "At Chrom:Pos " << Chrom << ":" << Pos << " have read " << N + 1 << " samples. Which is greater than  " <<  pops_of_indivs.size() << " from pops_of_indivs\n";
        Rcpp::stop("Insufficient number of pops_of_indivs.  Bailing out!");
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

      // now we sample one of the reads.  This is a little weird looking.  For individuals that
      // have only one type of read (or none) we just put their values into the first row of allele_mat.
      hetFlags(N) = 0L;
      if(rdr == 0 && rda == 0) {
        allele_mat(0, N) = -1L;
      } else if(rdr > 0 && rda == 0) {
        allele_mat(0, N) = 0L;
      } else if(rdr == 0 && rda > 0) {
        allele_mat(0, N) = 1L;
      } else {
        allele_mat(_, N) = runif(BootReps, 0, 1) < ((double) rda / (rdr + rda)); // sample alternate allele with prob rda/(rda+rdr)
        hetIdxs(numHets++) = N;
        hetFlags(N) = 1L;
      }

      if(alleles[N] != -1) {
        totAlt += alleles[N];
        totNotMissing++;
      }
      //Rcpp::Rcout << L << "  " << N << "  " << word << "  " << rdr << "  " << rda << "  " << comma << "      " << wdr << "  " << wda << "     " <<  alleles[N] << "\n";

    } // close loop that reads columns from each line
    N++;  // add one at the end.

    // throw an error if we didn't get the number expected.
    if(N != pops_of_indivs.size()) {
      Rcpp::Rcerr << "At Chrom:Pos " << Chrom << ":" << Pos << " just read " << N << " samples. Expected " <<  pops_of_indivs.size() << " from pops_of_indivs\n";
      Rcpp::stop("Didn't find the right number of samples.");
    }


      ///////////////////  HERE IS WHERE THE ACTUAL CALCULATIONS FOR EACH LOCUS OCCUR /////////

      // First we compute the Fixed and Rando counts for the populations (first level above individuals)
      // zero out to accumulate sums
      //Rcpp::Rcout << "Zeroing stuff\n";
      for(i=0;i<num_pops;i++) {
        Fixed0s(i) = 0L;
        Fixed1s(i) = 0L;
        for(r=0;r<BootReps;r++) {
          Rando0s(r, i) = 0L;
          Rando1s(r, i) = 0L;
        }
      }
      //Rcpp::Rcout << "Done Zeroing stuff\n";
      // Now, cycle over the individuals, and add things up according to which population each is in
      for(i=0;i<N;i++) {
        Fixed0s(pops_of_indivs(i)) += (hetFlags(i) == 0) && allele_mat(0, i) == 0;
        Fixed1s(pops_of_indivs(i)) += (hetFlags(i) == 0) && allele_mat(0, i) == 1;
      }
      for(i=0;i<numHets;i++) {
        for(r=0;r<BootReps;r++) {
          Rando0s(r, pops_of_indivs(hetIdxs(i))) += allele_mat(r, hetIdxs(i)) == 0;
          Rando1s(r, pops_of_indivs(hetIdxs(i))) += allele_mat(r, hetIdxs(i)) == 1;
        }
      }

      // COOL!  For all of the identity calculations above the level of individuals
      // we can now use Fixed0s, Fixed1s, and Rando0s and Rando1s to get the quantities
      // we need at each level and for each resample/bootstrap replicate.

      // determine how many times this locus is used in the poisson bootstrap samples
      bootn = rpois(BootReps, 1);
      bootn[0] = 1L;  // this is important.  The first row of zsums will be the actual estimate, and the remaining ones will be the bootstrap samples.

      // Here we cycle over all the pops and get the Y0 and Y1 matrices filled
      // for 0 up to num_pops - 1. And, while we are it, we will fill in the zsum values
      // for those populations
      for(r=0;r<BootReps;r++) {
        Ytot0(r) = 0L;
        Ytot1(r) = 0L;
        for(i=0;i<num_pops;i++) {  // cycle over pops once to fill the Y0 and Y1 matrices
          Y0(r, i) = Fixed0s(i) + Rando0s(r, i);
          Y1(r, i) = Fixed1s(i) + Rando1s(r, i);
          Ytot0(r) += Y0(r, i);
          Ytot1(r) += Y1(r, i);
        }

         // now, cycle over the pops again to compute zsums, etc.
        for(i=0;i<num_pops;i++) {
          // Now, only do something with polymorphic loci
          if( !(Ytot0(r) == 0 || Ytot1(r) == 0) ) {
            // we only actually add stuff to the sum when there are at least 2 gene copies that have been sampled
            // from each population.  We should come back and keep track of how many loci actually got used in there,
            // and also keep track of the averge sample size (number of individuals with reads) across all the internal
            // nodes.  Note that Y0(r, i) + Y1(r, i) should be constant over all r, so I could probably speed up this if() (though not really worth it, I expect)
            if((Y0(r, i) + Y1(r, i)) >= 2) {
              zbar =  (  Y0(r, i) * (Y0(r, i) - 1.0) +  Y1(r, i) * (Y1(r, i) - 1.0) ) /
                ( (Y0(r, i) + Y1(r, i) ) * (Y0(r, i) + Y1(r, i) - 1.0) );  // this is what gets added to the sum for this rep

              zsum(r, i) += bootn(r) * zbar;
              zL(r, i) += bootn(r);
            }
          }
        }
      }

      // NOW DEAL WITH THE REMAINING INTERNAL NODES
      // cycle over the bootstrap samples in the outer loop so that we can drop whole
      // loci at this level.
      for(r=0;r<BootReps;r++) {

        // And now we can cycle over the remaining internal nodes and compute Y0s and Y1s
        // for them, and add them to the zbars.
        if( !(Ytot0(r) == 0 || Ytot1(r) == 0) ) {
          for(j=0;j<num_internal_nodes - num_pops;j++) {
            i = j + num_pops;  // this is the index of the internal node
            daught = as<IntegerVector>(daughters(j));
            ndaught = daught.size();
            //Rcpp::Rcout << ndaught << ":" << daught << "\n";




            // first, go ahead and compute Y0 and Y1 for these internal nodes, as these
            // will get used for the parent of this node
            Y0(r, i) = 0;  // initialize to accumulate a sum
            Y1(r, i) = 0;
            for(k=0;k<ndaught;k++) {
              Y0(r, i) += Y0(r, daught[k]);
              Y1(r, i) += Y1(r, daught[k]);
            }

            // now compute the zbar and add it to the zsum
            zbar = 0.0;
            npairs = 0; // we will explicitly count them up, since some pairs might not have any data
            for(k1=0;k1<ndaught-1;k1++) {
              for(k2=k1+1;k2<ndaught;k2++) {
                d = daught[k1];
                dp = daught[k2];
                denom = ( (Y0(r, d) + Y1(r, d)) * (Y0(r, dp) + Y1(r, dp)) );
                if(denom > 0.0001) { // greater than 0, but these are ints coded as doubles, so give it a little buffer
                  zbar += (Y0(r, d) * Y0(r, dp) + Y1(r, d) * Y1(r, dp)) / denom;
                  npairs++;
                }
              }
            }
            // and then add that to zsum after dividing by the number of pairs and multiplying it by bootn
            if(npairs > 0.0001) {
              zsum(r, i) += bootn(r) * zbar / npairs;
              zL(r, i) += bootn(r);
            }
          }
        }
      }


      //return(List::create(bin_counts[0], bin_counts[1], bin_counts[2], bin_counts[3]));

     /////////////////// END OF ACTUAL CALCULATIONS FOR EACH LOCUS ///////////////////////////

    //return(List::create(allele_mat, numHets, hetIdxs, Fixed0s, Fixed1s, Rando0s, Rando1s, Y0, Y1, zsum, bootn, zL));
  } // close loop over markers (rows in file)


  // Now normalize the zsums and return all the stuff that we want to:
  for(i=0;i<BootReps;i++) {
    for(j=0;j<num_internal_nodes;j++) {
      zsum(i, j) /= zL(i, j);
    }
  }

  //return(List::create(zsum, bootn, zL));
  // store to return
  List ret;
  ret["zbar"] = zsum;
  ret["zL"] = zL;

  return(ret);
}

