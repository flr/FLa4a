//  --------------------------------------------------------------------------
// Copyright (c) 2008,2009,2010,2011,2012, 2013, 2014, 2015, 2016, 2017, 2018,
// Anders Nielsen <an@aqua.dtu.dk> and Colin Millar <colinpmillar@gmail.com>.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//   * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//   * Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//   * Neither the name of the assessment tool a4a nor the
//     names of its contributors may be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANDERS NIELSEN OR COLIN MILLAR BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//  --------------------------------------------------------------------------

GLOBALS_SECTION
  #include <time.h>
  time_t StartTime;   // Time variables, for use in timing loop
  #include <df1b2fun.h>
  #include <string>
  using namespace std;
  //#include "nLogNormal.h"
  //ofstream clogf("program.log");
  //#define TRACE(object) clogf<<"line "<<__LINE__<<", file "<<__FILE__<<", "\
  //                           <<#object" =\n"<<object<<endl<<endl;

  /* differentiable function for the likelihood */

  dvariable nldnorm(const double& obs, const dvariable& mu,
                    const dvariable& var) {
    return 0.5 * (log(2.0 * M_PI * var) + square(obs - mu) / var);
  }

  /* differentiable function for the likelihood */

  dvariable nldnorm(const dvariable& obs, const dvariable& mu,
                    const dvariable& var) {
    return 0.5 * (log(2.0 * M_PI * var) + square(obs - mu) / var);
  }


// *********************************
//
DATA_SECTION
//
// *********************************

  !! time(&StartTime);

  //
  // This first block is for the data
  //

  // age range and year range for analysis (inlcuding forecasts and hindcasts)
  init_vector ageRange(1,2) // full age range
  //!!TRACE(ageRange)
  init_vector yearRange(1,2) // full year range
  //!!TRACE(yearRange)
  // The number of surveys and when they take place */
  init_int nsurveys // number of surveys
  //!!TRACE (nsurveys)
  init_vector surveyMinAge(1,nsurveys) // min age of surveys
  //!!TRACE (surveyMinAge)
  init_vector surveyMaxAge(1,nsurveys) // max age of surveys
  //!!TRACE (surveyMaxAge)
  init_vector surveyTimes(1,nsurveys) // when does survey take place
  //!!TRACE (surveyTimes)
  // The fbar range and plus group information
  init_vector fbarRange(1,2) // fbar age range
  //!!TRACE(fbarRange)
  init_int isPlusGrp // is oldest age a plus group; 0 = NO; 1 = YES
  //!!TRACE(isPlusGrp)
  // The number of observations and the observation data
  init_int noobs // number of observations
  //!!TRACE(noobs)
  init_matrix obs(1,noobs,1,5) // fleet year age obs weight
  //!!TRACE(obs)
  // The number of auxilliary data points (we need this info for missing values
  // and predictions) and the auxilliary data. Explicit covariates will enter in
  // the design matrices... but this could pose a problem for prediction...
  init_int noaux // number of auxilliary data points
  //                ( should be diff(yearRange) * diff(ageRange) )
  //!!TRACE(noaux)
  init_matrix aux(1,noaux,1,7)
  //  year  age     m       m.spwn  harvest.spwn    mat * stkwt  stkwt
  //!!TRACE(aux)

  int idx

  //
  // useful summaries
  //
  int minYear
  !! minYear = yearRange(1);
  //!!TRACE(minYear)
  int maxYear
  !! maxYear = yearRange(2);
  //!!TRACE(maxYear)
  int minAge
  !! minAge = ageRange(1);
  //!!TRACE(minAge)
  int maxAge
  !! maxAge = ageRange(2);
  //!!TRACE(maxAge)

  //
  // The following blocks read the configuration files for the different
  // sub-models
  //

  // the Fishing mortality model

  !! ad_comm::change_datafile_name("fmodel.cfg");
  // First the fixed effects
  init_int noFpar
  //!!TRACE(noFpar)
  init_int noExpandedF
  //!!TRACE(noExpandedF)
  init_matrix designF(1,noExpandedF,1,noFpar)
  //!!TRACE(designF)
  // Then the (psuedo-)random effects vcov mat
  //init_int fpriorFlag
  //!! TRACE(fpriorFlag)
  //int fpriorPhase
  //!! if(fpriorFlag<0.5){fpriorPhase=-1;}else{fpriorPhase = 2;}
  //!! TRACE(fpriorPhase)
  //init_matrix covFP(1,noExpandedF,1,noExpandedF)
  //!! TRACE(covFP)


  // the Survey catchability model

  !! ad_comm::change_datafile_name("qmodel.cfg");
  // First the fixed effects
  init_int noQpar
  //!!TRACE(noQpar)
  init_int noExpandedQ
  //!!TRACE(noExpandedQ)
  init_matrix designQ(1,noExpandedQ,1,noQpar)
  //!!TRACE(designQ)
  // Then the (psuedo-)random effects vcov mat
  //init_int qpriorFlag
  //!! TRACE(qpriorFlag)
  //int qpriorPhase
  //!! if(qpriorFlag<0.5){qpriorPhase=-1;}else{qpriorPhase = 2;}
  //!! TRACE(qpriorPhase)
  //init_matrix covQP(1,noExpandedQ,1,noExpandedQ)
  //!! TRACE(covQP)


  // the variance model

  !! ad_comm::change_datafile_name("vmodel.cfg");
  // First the fixed effects
  init_int noVpar
  //!!TRACE(noVpar)
  init_int noExpandedV
  //!!TRACE(noExpandedV)
  init_matrix designV(1,noExpandedV,1,noVpar)
  //!!TRACE(designV)


  // the initial age structure model

  !! ad_comm::change_datafile_name("ny1model.cfg");
  // First the fixed effects
  init_int noNy1par
  //!!TRACE(noNy1par)
  init_int noExpandedNy1
  //!!TRACE(noExpandedNy1)
  init_matrix designNy1(1,noExpandedNy1,1,noNy1par)
  //!!TRACE(designNy1)

  // the internal recruitment model

  !! ad_comm::change_datafile_name("rmodel.cfg");
  // First the fixed effects
  init_int noRpar
  //!!TRACE(noRpar)
  init_int noExpandedR
  //!!TRACE(noExpandedR)
  init_matrix designR(1,noExpandedR,1,noRpar)
  //!!TRACE(designR)

  // the recruitment model

  !! ad_comm::change_datafile_name("srrmodel.cfg");
  // What model are we working with:
  init_int Rmodel
  //!!TRACE(Rmodel)
  // what CV is specified
  init_number srCV
  //!!TRACE(srCV)
  int SRaphase
  !! if (srCV > 0) { SRaphase = 2; } else { SRaphase = -1; }
  //!! TRACE(SRaphase)
  int SRbphase // swith of b if using geomean model
  !! if (srCV < 0 | Rmodel == 3) { SRbphase = -1; } else { SRbphase = 2; }
  !! if (srCV < 0 | Rmodel == 4) { SRbphase = -1; } else { SRbphase = 2; }
  //!! TRACE(SRbphase)
  init_number spr0 // only used with SV models
  //!! TRACE(spr0)
  // First the fixed effects for the a param (the level)
  init_int noRapar
  //!! TRACE(noRapar)
  init_int noExpandedRa
  //!! TRACE(noExpandedRa)
  init_matrix designRa(1,noExpandedRa,1,noRapar)
  //!! TRACE(designRa)
  // Than the fixed effects for the b param (the shape)
  init_int noRbpar
  //!! TRACE(noRbpar)
  init_int noExpandedRb
  //!! TRACE(noExpandedRb)
  init_matrix designRb(1,noExpandedRb,1,noRbpar)
  //!! TRACE(designRb)
  // Then the (psuedo-)random effects vcov mat on the a param
  //init_int rapriorFlag
  //!! TRACE(rapriorFlag)
  //int rapriorPhase
  //!! if(rapriorFlag<0.5){rapriorPhase=-1;}else{rapriorPhase = 2;}
  //!! TRACE(rapriorPhase)
  //init_matrix covRaP(1,noExpandedRa,1,noExpandedRa)
  //!! TRACE(covRaP)


  // MCMC report
  !!CLASS ofstream NMCMCreport("NMCMCreport.csv",ios::trunc);
  !!CLASS ofstream FMCMCreport("FMCMCreport.csv",ios::trunc);
  !!CLASS ofstream QMCMCreport("QMCMCreport.csv",ios::trunc);
  !!CLASS ofstream VMCMCreport("VMCMCreport.csv",ios::trunc);

// *********************************
//
PARAMETER_SECTION
//
// *********************************

  // the paramters of the fixed effects
  init_vector fpar(1,noFpar,3)
  init_vector qpar(1,noQpar,2)
  init_vector vpar(1,noVpar,1)
  init_vector ny1par(1,noNy1par,4)
  init_vector rpar(1,noRpar,1)
  init_vector rapar(1,noRapar,SRaphase)
  init_vector rbpar(1,noRbpar,SRbphase)

  // the parameters of the random effects
  //init_number logSdLogR(-1)
  //init_vector fdev(1,noExpandedF,fpriorPhase)
  //init_vector qdev(1,noExpandedQ,qpriorPhase)
  //init_vector radev(1,noExpandedRa,rapriorPhase)

  objective_function_value nll;

  // derived model quantities
  vector expandedF(1,noExpandedF)
  vector expandedQ(1,noExpandedQ)
  vector expandedV(1,noExpandedV)
  vector expandedNy1(1,noExpandedNy1)
  vector expandedR(1,noExpandedR)
  vector expandedRa(1,noExpandedRa)
  vector expandedRb(1,noExpandedRb)

  matrix  f(minYear,maxYear,minAge,maxAge)
  matrix  m(minYear,maxYear,minAge,maxAge)
  matrix  fspwn(minYear,maxYear,minAge,maxAge)
  matrix  mspwn(minYear,maxYear,minAge,maxAge)
  matrix  matWt(minYear,maxYear,minAge,maxAge)
  matrix  stkWt(minYear,maxYear,minAge,maxAge)
  matrix  logStkWt(minYear,maxYear,minAge,maxAge)
  3darray q(1,nsurveys,minYear,maxYear,minAge,maxAge)
  vector  r(minYear,maxYear) // predicted recruitment
  vector  ra(minYear,maxYear) // SRR params
  vector  rb(minYear,maxYear) // SRR params
  3darray v(1,nsurveys+1,minYear,maxYear,minAge,maxAge)
  matrix  n(minYear,maxYear,minAge,maxAge)
  vector  pred(1,noobs)

  vector ssb(minYear,maxYear)
  sdreport_number ssbmaxYear

//  vector fbar(minYear,maxYear)

//  sdreport_vector stateofstock(1,2)
//  sdreport_number r_last
//  sdreport_vector fpar_vec(1,noFpar)


// *********************************
//
INITIALIZATION_SECTION
//
// *********************************
fpar 0;
qpar 0;
vpar 0;
ny1par 0;
rpar 0;
rapar 0;
rbpar 0;

// *********************************
//
PRELIMINARY_CALCS_SECTION
//
// *********************************

// we might want to set up initial values for the stock recruitment model...
  //if(srCV>0){
  //  logSdLogR=log(srCV);
  //}

  //fpar=log(0.1)*max(designF);
  //qpar=log(1E-5)*max(designF);

  //
  // auxilliary data
  //
  idx = 0;
  for (int y = minYear; y <= maxYear; ++y) {
    for (int a = minAge; a <= maxAge; ++a) {
      idx = idx + 1;
      m(y,a) = aux(idx,3);
      mspwn(y,a) = aux(idx,4);
      fspwn(y,a) = aux(idx,5);
      matWt(y,a) = aux(idx,6);
      stkWt(y,a) = aux(idx,7);
      logStkWt(y,a) = log(aux(idx,7));
    }
  }



// *********************************
//
PROCEDURE_SECTION
//
// *********************************
  time_t currentTime;
  time(& currentTime);
  if (difftime(currentTime,StartTime) > 3600) { // Terminate after 60 minutes
    cerr<<endl;
    cerr<<"############################################################"<<endl;
    cerr<<"############################################################"<<endl;
    cerr<<"############################################################"<<endl;
    cerr<<"     MAX TIME ALLOWED EXCEEDED - MODEL DID NOT FINISH"       <<endl;
    cerr<<"############################################################"<<endl;
    cerr<<"############################################################"<<endl;
    cerr<<"############################################################"<<endl;
    cerr<<endl;
    ad_exit(1);
  }
  nll = 0.0;

  //
  // Fishing mortality model
  //
  expandedF = designF * fpar; //+ fdev;
  idx = 0;
  for (int y = minYear; y <= maxYear; ++y) {
    for (int a = minAge; a <= maxAge; ++a) {
      f(y,a) = expandedF(++idx);
    }
  }


  //
  // survey catchability model
  //
  expandedQ = designQ * qpar; //+ qdev;
  idx = 0;
  for (int ff = 1; ff <= nsurveys; ++ff) {
    for (int y = minYear; y <= maxYear; ++y) {
      for (int a = minAge; a <= maxAge; ++a) {
        q(ff,y,a) = expandedQ(++idx);
      }
    }
  }


  //
  // variance model
  //
  expandedV = designV * vpar;
  idx = 0;
  for (int ff = 1; ff <= nsurveys + 1; ++ff) {
    for (int y = minYear; y <= maxYear; ++y) {
      for (int a = minAge; a <= maxAge; ++a) {
        v(ff,y,a) = expandedV(++idx);
      }
    }
  }


  //
  // initial age structure model
  //
  expandedNy1 = designNy1 * ny1par;
  idx = 0;
  for (int a = minAge + 1; a <= maxAge; ++a) {
    n(minYear,a) = expandedNy1(++idx);
  }

  //
  // internal r model
  //
  expandedR = designR * rpar;
  idx = 0;
  for (int y = minYear; y <= maxYear; ++y) {
    r(y) = expandedR(++idx);
  }

  //
  // the full population structure
  //
  n.colfill(minAge,r);
  for (int a = minAge + 1; a <= maxAge; ++a) {
    for (int y = minYear + 1; y <= maxYear; ++y) {
      n(y,a) = n(y - 1,a - 1) -
                mfexp(f(y - 1,a - 1)) -
                mfexp(m(y - 1,a - 1));
      if ((a == maxAge) && (isPlusGrp > 0.5)) {
        n(y,a) = log(mfexp(n(y,a)) +
                      mfexp(n(y - 1,a) -
                            mfexp(f(y - 1,a)) -
                            mfexp(m(y - 1,a))
                            )
                      );
      }
    }
  }

  //
  // fbar and ssb
  //
  for (int y = minYear; y <= maxYear; ++y) {
    ssb(y) = sum(elem_prod(mfexp(n(y) -
                                 mfexp(f(y)) * fspwn(y) -
                                 mfexp(m(y)) * mspwn(y)),
                           matWt(y)));
//    fbar(y) = 0.0;
//    for(int a = fbarRange(1); a <= fbarRange(2); ++a) {
//      fbar(y) += mfexp(f(y,a));
//    }
//    fbar(y) /= fbarRange(2) - fbarRange(1) + 1;
  }

  ssbmaxYear = ssb(maxYear);
//  stateofstock(1) = ssb(maxYear);
//  stateofstock(2) = fbar(maxYear);

  //
  // The main likelihood
  //
  int locFleet, locSurvey, locYear, locAge;
  dvector obsVec(1,5);
  int minSurveyAge, maxSurveyAge;
  double locObs, locWgt;
  dvariable locZ;
  dvariable locVar;
  for (int i = 1; i <= noobs; ++i) {
    obsVec   = obs(i);
    locFleet = obsVec(1);
    locYear  = obsVec(2);
    locAge   = obsVec(3);
    locObs   = obsVec(4);
    locWgt   = obsVec(5);
    locSurvey = locFleet - 2;

    // here we split - if locAge == -1 then we have a biomass index
    // or total catch weight obs

    if (locAge >= 0)
    { // standard observation
      locZ = mfexp(f(locYear,locAge)) + mfexp(m(locYear,locAge));
      if (locFleet <= 2)
      { // catches
        pred(i) = f(locYear,locAge) -
                  log(locZ) +
                  log(1.0 - mfexp(-locZ)) +
                  n(locYear,locAge);
        locVar = mfexp(2.0 * v(locFleet,locYear,locAge));
      }
      else
      { // survey
        pred(i) = q(locSurvey,locYear,locAge) -
                  locZ * surveyTimes(locSurvey) +
                  n(locYear,locAge);
        locVar = mfexp(2.0 * v(locSurvey + 1,locYear,locAge));
      }
    }
    else
    { // if age is < 0, an observation of biomass / total catch weight has been specified
      if (locFleet <= 2)
      { // catches
        pred(i) = 0; // not sure i need to but best to be safe
        for (int a = minAge; a <= maxAge; ++a)
        {
          locZ = mfexp(f(locYear,a)) + mfexp(m(locYear,a));
          pred(i) +=
            mfexp(
              f(locYear,a) -
              log(locZ) +
              log(1.0 - mfexp(-locZ)) +
              n(locYear,a) +
              logStkWt(locYear,a)
            );
        }
        pred(i) = log(pred(i));
        // note variance are stored in the minimum age column for surveys
        // but what do we do for catch weights??
        locVar = mfexp(2.0 * v(locFleet-1,locYear,minAge));
      }
      else
      { // survey
        pred(i) = 0; // not sure i need to but best to be safe
        for (int a = surveyMinAge(locSurvey);
             a <= surveyMaxAge(locSurvey); ++a)
        {
          locZ = mfexp(f(locYear,a)) + mfexp(m(locYear,a));
          pred(i) +=
            mfexp(
              q(locSurvey,locYear,a) +
              logStkWt(locYear,a) +
              n(locYear,a) -
              surveyTimes(locSurvey) * locZ
            );
        }
        pred(i) = log(pred(i));
        // note variance are stored in the minimum age column for
        // biomass surveys -
        // but what do we do for catch weights??
        locVar = mfexp(2.0 * v(locSurvey + 1,locYear,minAge));
      }
    }
    nll += locWgt * nldnorm(locObs, pred(i), locVar);
  }


  //
  // stock recruit model
  //
  if (SRaphase > 0) { // then include a SRR model

    //
    // recruitment model
    //
    expandedRa = designRa * rapar; //+ radev;
    expandedRb = designRb * rbpar;
    idx = 0;
    for (int y = minYear; y <= maxYear; ++y) {
      idx = idx + 1;
      ra(y) = expandedRa(idx); // a and b are correlated though...
      rb(y) = expandedRb(idx);
    }

    //
    // weighted likelihood
    //
    dvariable predLogR;
    dvariable varLogR;
    dvariable h;
    dvariable v;

    if (Rmodel == 1) { // beverton holt
      for (int y = minYear + minAge; y <= maxYear; ++y) {
        predLogR = ra(y) +
                   log(ssb(y - minAge)) -
                   log(mfexp(rb(y)) + ssb(y - minAge));
        varLogR = log(pow(srCV, 2) + 1);
        nll += nldnorm(r(y), predLogR, varLogR);
      }
    }
    if (Rmodel == 2) { // ricker
      for (int y = minYear + minAge; y <= maxYear; ++y) {
        predLogR = ra(y) +
                   log(ssb(y - minAge)) -
                   mfexp(rb(y)) * ssb(y - minAge);
        varLogR = log(pow(srCV, 2) + 1);
        nll += nldnorm(r(y), predLogR, varLogR);
      }
    }
    if (Rmodel == 3) { // smooth hockey stick (Mesnil and Rochet, gamma = 0.1)
      for (int y = minYear + minAge; y <= maxYear; ++y) {
        predLogR = ra(y) +
                   log(ssb(y - minAge) +
                       sqrt(mfexp(2.0 * rb(y)) + 0.0025) -
                       sqrt(pow(ssb(y - minAge) - mfexp(rb(y)), 2.0) + 0.0025));
        varLogR = log(pow(srCV, 2) + 1);
        nll += nldnorm(r(y), predLogR, varLogR);
      }
    }
    if (Rmodel == 4) { // geomean
      for(int y = minYear + 1; y <= maxYear; ++y){
        predLogR = ra(y);
        varLogR = log(pow(srCV, 2) + 1);
        nll += nldnorm(r(y), predLogR, varLogR);
      }
    }
    if (Rmodel == 5) { // bevholt with steepness: ra is a transform of h; rb is a transform of v
      for(int y = minYear + minAge; y <= maxYear; ++y){
        h = mfexp(ra(y)) / (1 + mfexp(ra(y))) * 0.8 + 0.2;
        v = mfexp(rb(y));
        predLogR = log(6 * h * v * ssb(y - minAge)) -
                   log(spr0 * ((h + 1) * v + (5 * h - 1) * ssb(y - minAge))); // spr0 is provided by user
        varLogR = log(pow(srCV, 2) + 1);
        nll += nldnorm(r(y), predLogR, varLogR);
      }
    }
  }


  //
  // random bits as weighted likelihood components
  //
  //if (fpriorFlag > 0.5) {
  //  dvar_vector fzeroes(1, noExpandedF);
  //  fzeroes = 0;
  //  nll += nLogNormal(fdev, fzeroes, (dvar_matrix) covFP);
  //}

  //if (qpriorFlag > 0.5) {
  //  dvar_vector qzeroes(1, noExpandedQ);
  //  qzeroes = 0;
  //  nll += nLogNormal(qdev, qzeroes, (dvar_matrix) covQP);
  //}

  //if (rapriorFlag > 0.5) {
  //  dvar_vector razeroes(1, noExpandedRa);
  //  razeroes = 0;
  //  nll += nLogNormal(radev, razeroes, (dvar_matrix) covRaP);
  //}


  //
  // output the residuals
  //
  //if (sd_phase()) {
  //  ofstream res("a4a.res");
  //  res<<"fleet\tyear\tage\tobs\tpred\tsd\tres"<<endl<<residuals<<endl;
  //  res.close();
  // }

  //
  // MCMC report
  //
  if (mceval_phase()) {
    NMCMCreport << mfexp(n) << endl;
    FMCMCreport << mfexp(f) << endl;
    QMCMCreport << q << endl;
    VMCMCreport << v << endl;
  }


// *********************************
//
RUNTIME_SECTION
//
// *********************************

convergence_criteria 1E-1, 1E-2, 1E-3, 1E-12
maximum_function_evaluations 10, 20, 30, 10000

// *********************************
//
REPORT_SECTION
//
// *********************************

  ofstream nout("n.out");
  nout<<mfexp(n)<<endl;
  nout.close();
  ofstream fout("f.out");
  fout<<mfexp(f)<<endl;
  fout.close();
  ofstream qout("q.out");
  qout<<q<<endl;
  qout.close();
  ofstream vout("v.out");
  vout<<v<<endl;
  vout.close();
//  ofstream ssbout("ssb.out");
//  ssbout<<ssb<<endl;
//  ssbout.close();

// *********************************
//
TOP_OF_MAIN_SECTION
//
// *********************************

  arrmblsize=2000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(150000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(800000);
  gradient_structure::set_MAX_NVAR_OFFSET(100000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
