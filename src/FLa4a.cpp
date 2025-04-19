#include <TMB.hpp>

template <class Type>
vector<Type> ssbFUN(matrix<Type> logN, matrix<Type> logF, matrix<Type> M, matrix<Type> SW, matrix<Type> MO, matrix<Type> PF, matrix<Type> PM)
{
  int nrow = logN.rows();
  int ncol = logN.cols();
  vector<Type> ret(nrow);
  ret.setZero();
  for (int y = 0; y < nrow; ++y)
  {
    for (int a = 0; a < ncol; ++a)
    {
      ret(y) += SW(y, a) * MO(y, a) * exp(logN(y, a)) * exp(-PF(y, a) * exp(logF(y, a)) - PM(y, a) * M(y, a));
    }
  }
  return ret;
}

template <class Type>
Type objective_function<Type>::operator()()
{

  DATA_VECTOR(logObs)
  DATA_IMATRIX(aux)

  DATA_INTEGER(minYear)
  DATA_INTEGER(minAge)
  DATA_IVECTOR(surveyMinAges)
  DATA_IVECTOR(surveyMaxAges)
  DATA_IVECTOR(fleetTypes)
  DATA_VECTOR(sampleTimes)

  DATA_MATRIX(M)
  DATA_MATRIX(SW)
  DATA_MATRIX(MO)
  DATA_MATRIX(PF)
  DATA_MATRIX(PM)

  DATA_MATRIX(designF)
  DATA_MATRIX(designQ)
  DATA_MATRIX(designN1)
  DATA_MATRIX(designR)
  DATA_MATRIX(designV)
  DATA_MATRIX(designRa)
  DATA_MATRIX(designRb)

  // What model are we working with:
  DATA_INTEGER(RmodelId)
  DATA_SCALAR(spr0)
  DATA_INTEGER(estimateCV)

  // data summaries
  int nobs = logObs.size();
  int nrow = M.rows();
  int ncol = M.cols();
  int nsurvey = surveyMinAges.size();
  int hasSRR = RmodelId > 0;

  // PARAMETER_MATRIX(logN)
  PARAMETER_VECTOR(Fpar)
  PARAMETER_VECTOR(Qpar)
  PARAMETER_VECTOR(N1par)
  PARAMETER_VECTOR(Rpar)
  PARAMETER_VECTOR(Vpar)
  PARAMETER_VECTOR(Rapar)
  PARAMETER_VECTOR(Rbpar)
  PARAMETER(logSrCV)

  /// expand F
  vector<Type> expandedF(nobs);
  expandedF = designF * Fpar;

  matrix<Type> logF(nrow, ncol);
  for (int y = 0; y < nrow; ++y)
  {
    for (int a = 0; a < ncol; ++a)
    {
      logF(y, a) = expandedF(y * ncol + a);
    }
  }

  /// expand Q - block diagonal design matrix
  /// one block for each survey
  vector<Type> expandedQ(nobs);
  expandedQ = designQ * Qpar;
  array<Type> logQ(nsurvey, nrow, ncol);
  for (int s = 0; s < nsurvey; ++s)
  {
    for (int y = 0; y < nrow; ++y)
    {
      for (int a = 0; a < ncol; ++a)
      {
        logQ(s, y, a) = expandedQ(s * nrow * ncol + y * ncol + a);
      }
    }
  }

  /// expand V - block diagonal design matrix
  /// one block for each fleet + survey
  vector<Type> expandedV(nobs);
  expandedV = designV * Vpar;
  array<Type> logV(nsurvey + 1, nrow, ncol);
  for (int s = 0; s < nsurvey + 1; ++s)
  {
    for (int y = 0; y < nrow; ++y)
    {
      for (int a = 0; a < ncol; ++a)
      {
        logV(s, y, a) = expandedV(s * nrow * ncol + y * ncol + a);
      }
    }
  }

  /// expand N1
  vector<Type> logN1(nobs);
  logN1 = designN1 * N1par;

  /// expand R
  vector<Type> logR(nobs);
  logR = designR * Rpar;

  /// population
  matrix<Type> logN(nrow, ncol);

  // first year
  logN(0, 0) = logR(0);
  for (int a = 1; a < ncol; ++a)
  {
    logN(0, a) = logN1(a - 1);
  }

  // other years
  for (int y = 1; y < nrow; ++y)
  {
    logN(y, 0) = logR(y);
    for (int a = 1; a < ncol; ++a)
    {
      logN(y, a) = logN(y - 1, a - 1) - exp(logF(y - 1, a - 1)) - M(y - 1, a - 1);
      if (a == ncol - 1)
      {
        logN(y, a) = log(
            exp(logN(y, a)) +
            exp(logN(y - 1, a) - exp(logF(y - 1, a)) - M(y - 1, a)));
      }
    }
  }

  /// ssb
  vector<Type> ssb = ssbFUN(logN, logF, M, SW, MO, PF, PM);

  /// obs part

  vector<Type> logPred(nobs);
  vector<Type> logObsSd(nobs);

  Type Z;
  int y, a, f;
  for (int i = 0; i < nobs; ++i)
  {
    f = aux(i, 0) - 1;
    y = aux(i, 1) - minYear;
    a = aux(i, 2) - minAge;
    Z = exp(logF(y, a)) + M(y, a);
    logObsSd(i) = exp(logV(f, y, a));
    switch (fleetTypes(f))
    {
    case 1:
      logPred(i) = logN(y, a) - log(Z) + log(1 - exp(-Z)) + logF(y, a);
      break;

    case 2:
      logPred(i) = logQ(f - 1, y, a) + logN(y, a) - Z * sampleTimes(f - 1);
      break; // f - 1, because there is one commercial fleet

    default:
      std::cout << "Error: This fleet type not implemented yet." << std::endl;
      exit(EXIT_FAILURE);
      break;
    }
  }

  /// likelihood
  vector<Type> nllpart(nsurvey + 2);
  for (int i = 0; i < nllpart.size(); ++i)
  {
    nllpart(i) = 0;
  }
  for (int i = 0; i < nobs; ++i)
  {
    f = aux(i, 0) - 1;
    nllpart(f) += -dnorm(logObs(i), logPred(i), logObsSd(i), true); // add weights to variance
  }
  Type jnll = -sum(dnorm(logObs, logPred, logObsSd, true)); // add weights to variance

  /// sr model penalty
  if (hasSRR == 1)
  {
    // expand Ra and Rb
    vector<Type> logRa(nobs);
    logRa = designRa * Rapar;

    vector<Type> logRb(nobs);
    logRb = designRb * Rbpar;

    // careful with geomean! it uses all years
    vector<Type> predLogR(nrow - minAge);
    vector<Type> obsLogR(nrow - minAge);

    Type varSDPredR = pow(log(pow(exp(logSrCV), 2) + 1), 0.5);
    Type h;
    Type v;

    // calculate predictions of recruitment for year = minYear + minAge to maxYear
    // this misses the first minAge years.
    for (int y = 0; y < nrow - minAge; ++y)
    {
      // un penalised estimates of recriutment (our "observations")
      obsLogR(y) = logR(y + minAge);
      switch (RmodelId)
      {
      case 1:
        // beverton holt
        predLogR(y) = logRa(y + minAge) + log(ssb(y)) - log(exp(logRb(y + minAge)) + ssb(y));
        break;

      case 2:
        // ricker
        predLogR(y) = logRa(y + minAge) + log(ssb(y)) - exp(logRb(y + minAge)) * ssb(y);
        break;

      case 3:
        // smooth hockey stick (Mesnil and Rochet, gamma = 0.1)
        predLogR(y) = logRa(y + minAge) + log(ssb(y) + sqrt(exp(2.0 * logRb(y + minAge)) + 0.0025) - sqrt(pow(ssb(y) - exp(logRb(y + minAge)), 2.0) + 0.0025));
        break;

      case 4:
        // geometric mean
        predLogR(y) = logRa(y + minAge);
        break;

      case 5:
        // bevholt with steepness: ra is a transform of h; rb is a transform of v
        // spr0 is provided by user
        h = exp(logRa(y + minAge)) / (1 + exp(logRa(y + minAge))) * 0.8 + 0.2;
        v = exp(logRb(y + minAge));
        predLogR(y) = log(6 * h * v * ssb(y)) - log(spr0 * ((h + 1) * v + (5 * h - 1) * ssb(y)));
        break;

      default:
        std::cout
            << "Stock recruitment code not implemented yet." << std::endl;
        exit(EXIT_FAILURE);
        break;
      }
    }
    // calculate the likelihood of the recruitment estimates
    nllpart(nsurvey + 1) = -sum(dnorm(obsLogR, predLogR, varSDPredR, true));
    jnll += -sum(dnorm(obsLogR, predLogR, varSDPredR, true));

    if (RmodelId == 4)
    {
      // include in first years
      for (int y = 1; y < minAge; ++y)
      {
        nllpart(nsurvey + 1) += -dnorm(logR(y), logRa(y), varSDPredR, true);
        jnll += -dnorm(logR(y), logRa(y), varSDPredR, true);
      }
    }

    // if (estimateCV > 0)
    // {
    //   if (RmodelId == 4)
    //   {
    //     // include in first years
    //     for (int y = 0; y < minAge-1; ++y)
    //     {
    //       nllpart(nsurvey + 1) += -dnorm(logR(y), logRa(y), varSDPredR, true);
    //       jnll += -dnorm(logR(y), logRa(y), varSDPredR, true);
    //     }
    //   }
    //   else
    //   {
    //     first years are geomean if CV is being estimated
    //     for (int y = 0; y < minAge; ++y)
    //     {
    //       Type meanLogR = predLogR.sum() / predLogR.size();
    //       nllpart(nsurvey + 1) += -dnorm(logR(y), meanLogR, varSDPredR, true);
    //       jnll += -dnorm(logR(y), logRa(y), varSDPredR, true);
    //     }
    //   }
    // }
  }

  vector<Type> F = exp(expandedF);
  matrix<Type> N = exp(logN.array());
  vector<Type> Q = exp(expandedQ);

  REPORT(jnll);
  REPORT(nllpart)
  REPORT(logPred);
  REPORT(logObsSd);
  REPORT(logF);
  REPORT(logN1);
  REPORT(logR);
  REPORT(logN);
  REPORT(logQ);
  REPORT(logV);
  REPORT(logSrCV);

  ADREPORT(F);
  ADREPORT(N);
  ADREPORT(Q);

  return jnll;

  /* Quick Reference
     ===============

     ** Macros to read data and declare parameters:

     _Template_Syntax_              _C++_type_                     _R_type_
     DATA_VECTOR(name)              vector<Type>                   vector
     DATA_MATRIX(name)              matrix<Type>                   matrix
     DATA_SCALAR(name)              Type                           numeric(1)
     DATA_INTEGER(name)             int                            integer(1)
     DATA_FACTOR(name)              vector<int>                    factor
     DATA_SPARSE_MATRIX(name)       Eigen::SparseMatrix<Type>      dgTMatrix
     DATA_ARRAY(name)               array<Type>                    array
     PARAMETER_MATRIX(name)         matrix<Type>                   matrix
     PARAMETER_VECTOR(name)         vector<Type>                   vector
     PARAMETER_ARRAY(name)          array<Type>                    array
     PARAMETER(name)                Type                           numeric(1)

     ** Macro to report intermediate expressions back to R:

     REPORT(x)
     ADREPORT(x)

     ** Basic constructors:

     vector<Type> v(n1);
     matrix<Type> m(n1,n2);
     array<Type> a(n1,n2,n3)

     ** Basic operations:

     v+v,v-v,v*v,v/v                Pointwise binary operations
     m*v                            Matrix-vector multiply
     a.col(i)                       R equivalent of a[,,i]
     a.col(i).col(j)                R equivalent of a[,j,i]
     a(i,j,k)                       R equivalent of a[i,j,k]
     exp(v)                         Pointwise math
     m(i,j)                         R equivalent of m[i,j]
     v.sum()                        R equivalent of sum(v)
     m.transpose()                  R equivalent of t(m)

     ** Distributions:

     Type dnbinom2(const Type &x, const Type &mu, const Type &var, int give_log=0)
     Type dpois(const Type &x, const Type &lambda, int give_log=0)
     Type dlgamma(Type y, Type shape, Type scale, int give_log=0)
     Type dnorm(Type x, Type mean, Type sd, int give_log=0)

     ** Parallel accumulator declaration (only methods "+=" and "-="):

     parallel_accumulator<Type> res(this);

  */
}
