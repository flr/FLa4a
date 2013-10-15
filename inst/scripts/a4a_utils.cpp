#include <RcppArmadillo.h>
#include <cmath>


using namespace std;
using namespace arma;
using namespace Rcpp;


RcppExport SEXP matmul(SEXP Sfpar, SEXP SdesignF) 
{
BEGIN_RCPP

    arma::mat designF = as<arma::mat>(SdesignF);
    arma::vec fpar = as<arma::vec>(Sfpar);
    
    arma::mat expandedF = designF * fpar;

    return List::create( Named( "f" ) = expandedF);

END_RCPP
}


RcppExport SEXP calcLogLik(SEXP Sobs, SEXP SstockInfo, SEXP Sdesign, SEXP Spar) 
{
BEGIN_RCPP

    // read in observations and covariates (fleet, year, age)
    DataFrame DFobs(Sobs);
    arma::vec obs = log(as<arma::vec>(DFobs["obs"]));
    arma::uvec Fleet = as<arma::uvec>(DFobs["fleet"]);
    arma::uvec Year = as<arma::uvec>(DFobs["year"]);
    arma::uvec Age = as<arma::uvec>(DFobs["age"]);

    // dimensions of stock
    int noobs = obs.n_rows;
    int noFleets = max(Fleet) - min(Fleet) + 1;
    int noYears  = max(Year) - min(Year) + 1;
    int noAges   = max(Age) - min(Age) + 1;
    
    // convert for ease of use
    Age = Age - min(Age);
    Year = Year - min(Year);
    Fleet = Fleet - min(Fleet);


    // convert stock information: m, mat, stock weight etc.
    List LstockInfo(SstockInfo);
    arma::mat M  = as<arma::mat>(LstockInfo["m"]);
    arma::mat pm = as<arma::mat>(LstockInfo["mat"]);
    arma::mat sw = as<arma::mat>(LstockInfo["stock.wt"]);
    arma::mat cw = as<arma::mat>(LstockInfo["catch.wt"]);
    arma::vec fbar = as<arma::vec>(LstockInfo["fbar"]);
    arma::vec isPlusGrp = as<arma::vec>(LstockInfo["plusgroup"]);
    arma::vec surveyTimes = as<arma::vec>(LstockInfo["surveytime"]);


    // convert design matrices
    List Ldesign(Sdesign);
    arma::mat designF = as<arma::mat>(Ldesign["Xf"]);
    arma::mat designQ = as<arma::mat>(Ldesign["Xq"]);
    arma::mat designR = as<arma::mat>(Ldesign["Xr"]);
    arma::mat designCV = as<arma::mat>(Ldesign["Xcv"]);

    
    // convert parameters NOTE all are on the log scale
    List Lpar(Spar);
    arma::mat fpar = as<arma::mat>(Lpar["f"]);
    arma::mat qpar = as<arma::mat>(Lpar["q"]);
    arma::mat rpar = as<arma::mat>(Lpar["r"]);
    arma::mat ry1 = as<arma::mat>(Lpar["ry1"]);
    arma::mat cvpar = as<arma::mat>(Lpar["cv"]);
    int nofits = fpar.n_rows;


    //
    // Now for the model
    //
    arma::vec pred = zeros<vec>(noobs,1);
    arma::mat n(noAges, noYears);
    arma::vec ll = zeros<vec>(nofits,1);

    for (int j = 0; j < nofits; j++) {

    // do F prediction
    arma::mat f = designF * fpar.row(j).t();
    f.reshape(noAges, noYears);

    // do q prediction: result is a n x 1 matrix, assigned to a slice of a cube. we then reshape the cube.
    arma::cube q(noYears * noAges * (noFleets - 1), 1, 1);
    q.slice(0) = designQ * qpar.row(j).t();
    q.reshape(noAges, noYears, noFleets - 1);

    // do R prediction
    arma::mat r = designR * rpar.row(j).t();

    // do CV prediction: result is a n x 1 matrix, assigned to a slice of a cube. we then reshape the cube.
    arma::cube cv(noYears * noAges * noFleets, 1, 1);
    cv.slice(0) = designCV * cvpar.row(j).t();
    cv.reshape(noAges, noYears, noFleets);

    // estimate n at age
    n.each_row() = r.t();
    //n.cols(0) = join_cols(r(0),ry1);
    for (int a = 1; a < noAges; a++) n(a,0) = ry1(a-1);

    for(int a = 1; a < noAges; a++) {
      for(int y = 1; y < noYears; y++) {
        n(a,y) = n(a-1,y-1) - exp(f(a-1,y-1)) - M(a-1,y-1);
        if( (a+1 == noAges) && (isPlusGrp(0) > 0.5) ) {
          n(a,y) = log( exp(n(a,y)) + exp( n(a,y-1) - exp(f(a,y-1)) - M(a,y-1) ));
        }
      }
    }

    // predict catches and indices and calculate likelihood
    arma::vec locll = zeros<vec>(1,1);
    arma::vec locZ(1), locCV(1);
    for(int i = 0; i < noobs; i++)
    {
      locZ = exp(f(Age(i), Year(i))) + M(Age(i),Year(i));
      if (Fleet(i) == 0) { //    catches predicted 
        pred(i) = f(Age(i), Year(i)) - log(locZ(0)) + log(1.0 - exp(-locZ(0))) + n(Age(i), Year(i));
      } else {           //    survey predicted 
        pred(i) = q(Age(i), Year(i), Fleet(i) - 1) - locZ(0) * surveyTimes(Fleet(i) - 1) + n(Age(i),Year(i)); 
      }
      locCV =  cv(Age(i), Year(i), Fleet(i));
      locll += -0.5 * log(2.0 * datum::pi) - locCV - 0.5 * exp(-2.0 * locCV) * pow(obs(i) - pred(i), 2.0);    
    }

    ll(j) = locll(0);

    }
    
    return List::create(Named("ll") = ll);

END_RCPP
}



//    NumericVector surveyTimes(Ssurvtime);

//    NumericVector SRcv(SSRcv);

//    NumericVector fpriorFlag(SfpriorFlag);
//    NumericVector qpriorFlag(SqpriorFlag);
//    NumericMatrix covFP(ScovFP);


//    double rec_loga, rec_logb, logSdLogR; //(-1)
//    double fdev[maxYear-minYear+1][maxAge-minAge+1];


   //  Parameters are defined!  Now some initialising!!

   //  Now the model!!


// if(fpriorFlag > .5) {
//    dvar_vector fzeroes(ageRange(1),ageRange(2));
//    fzeroes=0;
//    nll+=sum(nLogNormal(trans(fdev),fzeroes,(dvar_matrix)covFP));
//  }

//  if(qpriorFlag>.5){
//    for(int ff=2; ff<=maxFleet; ++ff){
//      dvar_vector qzeroes(ageRangeS(ff,1),ageRangeS(ff,2));
//      qzeroes=0;
//      nll+=sum(nLogNormal(trans(qdev(ff)),qzeroes,(dvar_matrix)covQP(ff)));
//    }
//  }


//  if(SRphase>0){ 
//    dvariable predLogR; 
//    dvariable varLogR; 
//    for(int y=minYear+1; y<=maxYear; ++y){
//      predLogR=rec_loga+log(ssb(y-1))-log(exp(rec_logb)+ssb(y-1));
//      varLogR=exp(2.0*logSdLogR);
//      nll+=nldnorm(r(y),predLogR,varLogR);    
//    }
//  }


//  imatrix ageRangeQP(2,maxFleet,1,2)
//  ivector minAgeQP=column(ageRangeS,1);
//  ivector maxAgeQP=column(ageRangeS,2);
//  3darray covQP(2,maxFleet,minAgeQP,maxAgeQP,minAgeQP,maxAgeQP)



