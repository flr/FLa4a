/*
* Copyright (c) 2012-2014 European Union
* European Commission Joint Research Centre G.04.
*
* Authors: Anders Nielsen <an@aqua.dtu.dk>
* 	       Colin Millar Colin Millar <millarc@marlab.ac.uk>
*
* Licensed under the EUPL, Version 1.1 or later (the "Licence");
* 
* You may not use this work except in compliance with the Licence.
* You may obtain a copy of the Licence at:
*
* http://ec.europa.eu/idabc/eupl
*
* Unless required by applicable law or agreed to in writing,
* software distributed under the Licence is distributed on an "AS IS" basis,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the Licence for the specific language governing permissions and limitations
* under the Licence.
*/


////snip001////
GLOBALS_SECTION 
  #include <time.h>
  time_t StartTime;   //Time variables, for use in timing loop
  #include <df1b2fun.h>
  #include <string>
  using namespace std;
  #include "nLogNormal.h"
  ofstream clogf("program.log");
  #define TRACE(object) clogf<<"line "<<__LINE__<<", file "<<__FILE__<<", "\
                             <<#object" =\n"<<object<<endl<<endl; 

  dvariable nldnorm(const double& obs, const dvariable& mu, 
                    const dvariable& var){
    return 0.5*(log(2.0*M_PI*var)+square(obs-mu)/var);
  }

  dvariable nldnorm(const dvariable& obs, const dvariable& mu, 
                    const dvariable& var){
    return 0.5*(log(2.0*M_PI*var)+square(obs-mu)/var);
  }
////end-snip001////

////doc001////
  /* 
  The \verb|GLOBAL\_SECTION| describes variables and functions which can be called anywhere in the AD Model Builder program. 
  Standard $C^{++}$ functions for computing running times (\verb|time.h|) and for manipulating text strings (\verb|string|) 
  are included in the program. 
  A custom built library (\verb|nLogNormal.h|) to efficiently evaluate the multivariate normal distribution is included. A global macro 
  is defined called \verb|TRACE|, which allows the developers to monitor data objects, parameters, and other derived quantities anywhere 
  in the program. Finally, the negative log likelihood density function (\verb|nldnorm|) of the uni-variate normal distribution is 
  defined (for two calling situations), which makes it easily accessible in the rest of the program. 
  */ 
////end-doc001////

////snip002////
DATA_SECTION
  !! time(&StartTime);
  init_int noobs
  !! TRACE(noobs)
  init_matrix obs(1,noobs,1,4) // f y a o 
  !! TRACE(obs)
  !!////end-snip002////
  !!////doc002////
  !!/* 
  !!  The \verb|DATA\_SECTION| is where all the data objects are read in from the data files. 
  !!  The main data object here is a big matrix (\verb|obs|) containing as many rows as there are 
  !!  catch and survey observations available (\verb|noobs|). The columns of the matrix are 1) fleet number with the convention 
  !!  that fleet number 1 is catches,  2) observation year,  3) observation age, and 4) the observation 
  !!  itself. The first line of the \verb|DATA\_SECTION| initializes the running time used internally by the program. 
  !!  Wherever the \verb|TRACE| macro is called it writes the object to a file called \verb|program.log|, which can 
  !!  then later be used to validate that all inputs are perceived correctly by the program. The double \verb|!| is used 
  !!  here to indicate that the code is to be interpreted as straight $C^{++}$ code.        
  !!*/ 
  !!////end-doc002////
  !!////snip003////
  int minYear
  !! minYear=min(column(obs,2));
  !! TRACE(minYear)
  int maxYear
  !! maxYear=max(column(obs,2));
  !! TRACE(maxYear)
  int minAge
  !! minAge=min(column(obs,3));
  !! TRACE(minAge)
  int maxAge
  !! maxAge=max(column(obs,3));
  !! TRACE(maxAge)
  int noAges
  !! noAges=maxAge-minAge+1;
  !! TRACE(noAges)
  int minFleet
  !! minFleet=min(column(obs,1));
  !! TRACE(minFleet)
  int maxFleet
  !! maxFleet=max(column(obs,1));
  !! TRACE(maxFleet)
  !!////end-snip003////
  !!////doc003////
  !!/* 
  !!  From the main input matrix \verb|obs|, a number of useful summary variables are now extracted. 
  !!  As an example observe that the minimum observation year is computed as the minimum of the second column. 
  !!  Nothing new is read in from the files here, but just declared and computed for convenience.   
  !!*/ 
  !!////end-doc003////

  !!////snip004////
  init_matrix natMor(minYear,maxYear,minAge,maxAge)
  !! TRACE(natMor)
  matrix m(minYear,maxYear,minAge,maxAge)
  !! m=log(natMor);
  !! TRACE(m)
  init_matrix mat(minYear,maxYear,minAge,maxAge)
  !! TRACE(mat)
  init_matrix cw(minYear,maxYear,minAge,maxAge)
  !! TRACE(cw)
  init_matrix sw(minYear,maxYear,minAge,maxAge)
  !! TRACE(sw)
  !!////end-snip004////
  !!////doc004////
  !!/* 
  !!  The natural mortality, the maturity, the catch weights, and the stock weights are now read in as simple matrices.       
  !!  The natural mortality on logarithmic scale is calculated and stored in a matrix called \verb|m|.    
  !!*/ 
  !!////end-doc004////

  !!////snip005////
  !! ad_comm::change_datafile_name("model.cfg");
  init_matrix varKey(1,maxFleet,minAge,maxAge);
  !! TRACE(varKey)
  int noVarPar
  !! noVarPar=max(varKey);
  !! TRACE(noVarPar)
  init_vector surveyTimes(2,maxFleet)
  !! TRACE(surveyTimes)
  init_vector fbarRange(1,2)
  !! TRACE(fbarRange)
  init_int isPlusGrp
  !! TRACE(isPlusGrp)
  init_number SRcv;
  !! TRACE(SRcv)
  int SRphase
  !! if(SRcv<0){SRphase=-1;}else{SRphase = 2;}
  !! TRACE(SRphase)
  init_int fpriorFlag
  !! TRACE(fpriorFlag)
  int fpriorPhase
  !! if(fpriorFlag<0.5){fpriorPhase=-1;}else{fpriorPhase = 2;}
  !! TRACE(fpriorPhase)
  init_int qpriorFlag
  !! TRACE(qpriorFlag)
  int qpriorPhase
  !! if(qpriorFlag<0.5){qpriorPhase=-1;}else{qpriorPhase = 2;}
  !! TRACE(qpriorPhase)
  !!////end-snip005////
  !!////doc005////
  !!/* 
  !! The input file is now changed to \verb|model.cfg|, which is a file containing a few model configuration 
  !! options, or other inputs. The first thing read in is the observation variance parameter coupling matrix, 
  !! which maps each fleet and age combination to its corresponding variance parameter. The survey times, 
  !! as fractions into the year. The range of the average fishing mortalities reported is read in, followed by
  !! whether the oldest age is a plus group. The a number called the CV of the stock 
  !! recruitment relationship is set. This number (if positive) sets the coefficient of variation in a penalized 
  !! likelihood approach, and if negative turns the assumption of a stock recruitment relationship off, which 
  !! in AD Model Builder is obtained by setting the estimation phase to a negative number. Similarly a flag is 
  !! set for turning on and off the prior on the fishing mortalities and catchabilities.     
  !!*/ 
  !!////end-doc005////

  !!////snip006////
  !! ad_comm::change_datafile_name("fmodel.cfg");
  init_ivector ageRange(1,2)
  !! TRACE(ageRange)
  init_ivector yearRange(1,2)
  !! TRACE(yearRange)
  init_int noFpar
  !! TRACE(noFpar)
  int noExpandedF
  !! noExpandedF=(ageRange(2)-ageRange(1)+1)*(yearRange(2)-yearRange(1)+1);
  !! TRACE(noExpandedF)
  init_matrix designF(1,noExpandedF,1,noFpar)
  !! TRACE(designF)
  !!////end-snip006////
  !!////doc006////
  !!/* 
  !! The input file is now changed to \verb|fmodel.cfg|, which is a file containing the configuration      
  !! for the fishing mortality model assumed. The file is created outside of the AD Model builder program 
  !! and contains description of the dimensions, and a design matrix mapping the parameters into the fishing  
  !! mortality of each year and age combination. 
  !!*/ 
  !!////end-doc006//// 

  !!////snip006a////
  !! ad_comm::change_datafile_name("fprior.cfg");
  init_ivector ageRangeFP(1,2)
  !! TRACE(ageRangeFP)
  int dimFP
  !! dimFP=ageRangeFP(2)-ageRangeFP(1)+1;
  !! TRACE(dimFP)
  init_matrix covFP(1,dimFP,1,dimFP)
  !! TRACE(covFP)
  !!////end-snip006a////
  !!////doc006a////
  !!/* 
  !! This code chunk reads in the \verb|fprior.cfg| configuration file, which sets the prior covariance 
  !! of the fishing mortality deviations, if activated. Activation is via the \verb|model.cfg| file.  The current 
  !! implementation requires a positive definite covariance matrix which provides the model to be applied
  !! to each year separately.  
  !!*/ 
  !!////end-doc006a//// 


  !!////snip007////
  imatrix ageRangeS(2,maxFleet,1,2)
  imatrix yearRangeS(2,maxFleet,1,2)
  ivector noQpar(2,maxFleet)
  ivector noExpandedQ(2,maxFleet)
  !! ifstream qin("qmodel.cfg");
  !! string tmpLine;
  !! getline(qin,tmpLine); 
  !! for(int f=2; f<=maxFleet; ++f){
  !!   for(int i=1; i<=2; ++i){getline(qin,tmpLine);}
  !!   qin>>ageRangeS(f);
  !!   for(int i=1; i<=2; ++i){getline(qin,tmpLine);}
  !!   qin>>yearRangeS(f);
  !!   for(int i=1; i<=2; ++i){getline(qin,tmpLine);}
  !!   qin>>noQpar(f);
  !!   for(int i=1; i<=2; ++i){getline(qin,tmpLine);}
  !!   noExpandedQ(f)=((ageRangeS(f,2)-ageRangeS(f,1)+1)*(yearRangeS(f,2)-yearRangeS(f,1)+1));
  !!   for(int i=1; i<=noExpandedQ(f); ++i){getline(qin,tmpLine);}
  !! }
  !! TRACE(ageRangeS)
  !! TRACE(yearRangeS)
  !! TRACE(noQpar)
  !! qin.close();
  !! ivector ones(2,maxFleet); 
  !! ones=1;
  3darray designQ(2,maxFleet,ones,noExpandedQ,ones,noQpar) 
  !! qin.open("qmodel.cfg");
  !! getline(qin,tmpLine); 
  !! for(int f=2; f<=maxFleet; ++f){
  !!   for(int i=1; i<=8; ++i){getline(qin,tmpLine);}   
  !!   qin>>designQ(f);     
  !!   getline(qin,tmpLine);
  !! }
  !! TRACE(designQ)
  !! qin.close();
  !!////end-snip007////
  !!////doc007////
  !!/* 
  !! This code chunk achieves the same for the catchability model setup for each survey fleet, as the previous    
  !! code chunk did for the fishing mortality. The input is read from the file \verb|qmodel.cfg|. The code 
  !! is a little bit more low-level $C^{++}$, because the logical format of the input (listing one fleet 
  !! configuration after another did not allow all dimensions to be known prior to the time where the \verb|3darray| 
  !! had to be defined. This forced the implementation above where the file is passed twice. First to collect 
  !! the dimension information, next to read in all the separate design matrices.      
  !!*/ 
  !!////end-doc007////

  !!////snip007a////
  imatrix ageRangeQP(2,maxFleet,1,2)
  !! ifstream qpin("qprior.cfg");
  !! getline(qpin,tmpLine); 
  !! for(int f=2; f<=maxFleet; ++f){
  !!   for(int i=1; i<=2; ++i){getline(qpin,tmpLine);}
  !!   qpin>>ageRangeQP(f);
  !!   for(int i=1; i<=2; ++i){getline(qpin,tmpLine);}
  !!   for(int i=ageRangeQP(f,1); i<=ageRangeQP(f,2); ++i){getline(qpin,tmpLine);}
  !! }
  !! TRACE(ageRangeQP)
  !! qpin.close();
  !! ivector minAgeQP=column(ageRangeS,1);
  !! ivector maxAgeQP=column(ageRangeS,2);
  3darray covQP(2,maxFleet,minAgeQP,maxAgeQP,minAgeQP,maxAgeQP)
  !! qpin.open("qprior.cfg");
  !! getline(qpin,tmpLine); 
  !! for(int f=2; f<=maxFleet; ++f){
  !!   for(int i=1; i<=4; ++i){getline(qpin,tmpLine);}   
  !!   qpin>>covQP(f);     
  !!   getline(qpin,tmpLine);
  !! }
  !! TRACE(covQP)
  !! qpin.close();
  !!////end-snip007a////
  !!////doc007a////
  !!/* 
  !! This code chunk reads in the \verb|qprior.cfg| configuration file, which sets the prior covariance 
  !! of the catchability deviations, if activated. Activation is via the \verb|model.cfg| file.  The current 
  !! implementation requires a positive definite covariance matrix which provides the model be applied
  !! to each year separately.  
  !!*/ 
  !!////end-doc007a////

  !!////snip007b////
  !! ad_comm::change_datafile_name("rmodel.cfg");
  init_ivector yearRangeR(1,2)
  !! TRACE(yearRangeR)
  init_int noRpar
  !! TRACE(noRpar)
  int noExpandedR
  !! noExpandedR=yearRangeR(2)-yearRangeR(1)+1;
  !! TRACE(noExpandedR)
  init_matrix designR(1,noExpandedR,1,noRpar)
  !! TRACE(designR)
  !!////end-snip007b////
  !!////doc007b////
  !!/* 
  !! This code chunk achieves the same for the recruitment model setup, as the previous    
  !! code chunks did for the fishing mortality and survey catchability. The input is read from the file \verb|rmodel.cfg|. 
  !! The code is almost exactly the same as for the fishing mortality model      
  !!*/ 
  !!////end-doc007b////


  !!////snip008////
  matrix residuals(1,noobs,1,7)
  !!////end-snip008////
  !!////doc008////
  !!/* 
  !! A matrix is defined to store the residual information, such that it can be written to a file 
  !! after the model has finished. The columns of this are: fleet, year, age, observation, prediction, 
  !! standard deviation, and normalized residual. Observation, prediction, standard deviation, and residual 
  !! are on logarithmic scale. 
  !!*/ 
  !!////end-doc008////
 
  !!////snip009////
PARAMETER_SECTION
  init_vector rpar(1,noRpar)
  init_vector ry1(minAge+1,maxAge)
  init_matrix qpar(2,maxFleet,1,noQpar)
  init_vector fpar(1,noFpar)
  init_vector logSdlogObs(1,noVarPar)
  init_number rec_loga(SRphase)
  init_number rec_logb(SRphase)
  init_number logSdLogR(-1)
  init_matrix fdev(minYear,maxYear,minAge,maxAge,fpriorPhase)

  !!////end-snip009////
  !!////doc009////
  !!/* 
  !! The \verb|PARAMETER\_SECTION| is where all model parameters are defined. The parameters of this model are the numbers at are for the  
  !! first year a and for the first age, catchability parameters, fishing mortality parameters, observation standard deviations, and 
  !! parameters of the stock recruitment (if active). 
  !!*/ 
  !!////end-doc009////

  !!////snip010////
  objective_function_value nll;
  !!////end-snip010////
  !!////doc010////
  !!/* 
  !! This line, which is part of all AD Model Builder programs, is where the name function of the function to be minimized is defined. 
  !!*/ 
  !!////end-doc010////

  !!////snip011////
  vector expandedF(1,noExpandedF)
  matrix expandedQ(2,maxFleet,1,noExpandedQ)
  vector expandedR(1,noExpandedR)
  matrix f(minYear,maxYear,minAge,maxAge)
  !! ivector minYearS=column(yearRangeS,1);
  !! ivector maxYearS=column(yearRangeS,2);
  !! ivector minAgeS=column(ageRangeS,1);
  !! ivector maxAgeS=column(ageRangeS,2);
  3darray q(2,maxFleet,minYearS,maxYearS,minAgeS,maxAgeS)
  init_3darray qdev(2,maxFleet,minYearS,maxYearS,minAgeS,maxAgeS,qpriorPhase)
  vector r(minYear,maxYear)
  matrix n(minYear,maxYear,minAge,maxAge)
  vector pred(1,noobs)
  !!////end-snip011////
  !!////doc011////
  !!/* 
  !! In addition to defining parameters to be estimated the \verb|PARAMETER\_SECTION| can also be used to define variables 
  !! to store intermediate calculations. Here are container objects set up to hold fishing mortalities, catchabilities, 
  !! stock numbers, and predictions corresponding to the observations.  
  !!*/ 
  !!////end-doc011////

  !!////snip012////
  sdreport_vector ssb(minYear,maxYear)
  sdreport_vector fbar(minYear,maxYear)
  !!////end-snip012////
  !!////doc012////
  !!/* 
  !! Finally two vectors are set up to hold the spawning stock biomass and average fishing mortalities. These are set up 
  !! as \verb|sdreport_vector|, because that causes the standard deviations to be calculated and reported for these  
  !! quantities.
  !!*/ 
  !!////end-doc012////

  !! ////snip013////
PRELIMINARY_CALCS_SECTION
  
  rec_loga=12; 
  rec_logb=0;
  if(SRcv>0){
    logSdLogR=log(SRcv);
  }
  fdev.initialize();
  qdev.initialize();
  ////end-snip013////
  ////doc013////
  /* 
  The \verb|PRELIMINARY\_CALCS\_SECTION| is used to set initial values and possibly make a few calculations before 
  the parameter fitting loop is started. Here only a few initial values for the stock recruitment function 
  is set. For the main part of the parameters it was not necessary to set initial values for this model.
  \verb|fdev| and \verb|qdev| are initialised to zero here in case they are not used.
  */ 
  ////end-doc013////


  ////snip014////
PROCEDURE_SECTION
  time_t currentTime;
  time(& currentTime);
  if(difftime(currentTime,StartTime)>600){ // Terminate after 10 minutes 
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
  nll=0.0;
  ////end-snip014////
  ////doc014////
  /* 
  The \verb|PROCEDURE\_SECTION| is where the negative log-likelihood function $\ell(\theta,X)$ is calculated.
  Here a check is done for runtime, and if that exceeds 10 minutes the program is terminated with an 
  error message. The negative log-likelihood is initialized to zero before the different terms are added.      
  */ 
  ////end-doc014////

  ////snip015////
  if(fpriorFlag>.5){
    dvar_vector fzeroes(ageRange(1),ageRange(2));
    fzeroes=0;
    nll+=sum(nLogNormal(trans(fdev),fzeroes,(dvar_matrix)covFP));
  }

  if(qpriorFlag>.5){
    for(int ff=2; ff<=maxFleet; ++ff){
      dvar_vector qzeroes(ageRangeS(ff,1),ageRangeS(ff,2));
      qzeroes=0;
      nll+=sum(nLogNormal(trans(qdev(ff)),qzeroes,(dvar_matrix)covQP(ff)));
    }
  }
  ////end-snip015////
  ////doc015////
  /* 
  Here the fishing mortality and or the catchability deviation prior is added --- if activated.
  The efficiency of this part relies heavily on code in the included header file to efficiently
  compute the multivariate negative log-density for a normal distribution.   
  */ 
  ////end-doc015////

  ////snip016////
  expandedF=designF*fpar; 
  int idx=0;
  for(int a=minAge; a<=maxAge; ++a){
    for(int y=minYear; y<=maxYear; ++y){
      f(y,a)=expandedF(++idx) + fdev(y,a);
    }
  }
  ////end-snip016////
  ////doc016////
  /* 
  The fishing mortality for all years and ages are computed here, which is based on the design matrix, which is read in and the 
  actual F parameters estimated in the program. 
  */ 
  ////end-doc016////

  ////snip017////  
  for(int ff=2; ff<=maxFleet; ++ff){
    expandedQ(ff)=designQ(ff)*qpar(ff); 
    idx=0;
    for(int a=ageRangeS(ff,1); a<=ageRangeS(ff,2); ++a){
      for(int y=yearRangeS(ff,1); y<=yearRangeS(ff,2); ++y){
        q(ff,y,a)=expandedQ(ff,++idx) + qdev(ff,y,a);
      }
    }
  }  
  ////end-snip017////
  ////doc017////
  /* 
  The catchabilities for all survey fleets, years and ages are computed here. This is based on the design matrices, which are read in and the 
  actual catchability parameters estimated in the program. 
  */ 
  ////end-doc017////

  ////snip018////
  expandedR=designR*rpar; 
  idx=0;
  for(int y=minYear; y<=maxYear; ++y){
    r(y)=expandedR(++idx);
  }
  ////end-snip018////
  ////doc018////
  /* 
  The recruitment for all years are computed here. This is based on the read in design matrix and the recruitment parameters estimated in the program. 
  */ 
  ////end-doc018////

  ////snip019////
  n.colfill(minAge,r);
  for(int a=minAge+1; a<=maxAge; ++a){
    n(minYear,a)=ry1(a);
    for(int y=minYear+1; y<=maxYear; ++y){
      n(y,a)=n(y-1,a-1)-exp(f(y-1,a-1))-exp(m(y-1,a-1));
      if((a==maxAge) && (isPlusGrp>0.5)){
        n(y,a)=log(exp(n(y,a))+exp(n(y-1,a)-exp(f(y-1,a))-exp(m(y-1,a))));
      }
    }
  }
  ////end-snip019////
  ////doc019////
  /* 
  All N-at-age are computed via the stock equation from the recruits and the first years N's which are the model parameters.
  This is where the plus group calculations take place if the plus group flag \verb|isPlusGrp| is set to 1 
  */ 
  ////end-doc019////

  ////snip020//// 
  for(int y=minYear; y<=maxYear; ++y){
    ssb(y)=sum(elem_prod(elem_prod(exp(n(y)),sw(y)),mat(y)));
    fbar(y)=0.0;
    for(int a=fbarRange(1); a<=fbarRange(2); ++a){
      fbar(y)+=exp(f(y,a));
    }
    fbar(y)/=fbarRange(2)-fbarRange(1)+1;
  }  
  ////end-snip020////
  ////doc020////
  /* 
  The spawning stock biomass and average fishing mortalities are computed from the calculated stock sizes and fishing mortalities. 
  */ 
  ////end-doc020////

  ////snip021////
  if(SRphase>0){ 
    dvariable predLogR; 
    dvariable varLogR; 
    for(int y=minYear+1; y<=maxYear; ++y){
      predLogR=rec_loga+log(ssb(y-1))-log(exp(rec_logb)+ssb(y-1));
      varLogR=exp(2.0*logSdLogR);
      nll+=nldnorm(r(y),predLogR,varLogR);    
    }
  }
  ////end-snip021////
  ////doc021////
  /* 
  If activated, the stock recruitment relationship is penalized towards a Beverton-Holt relationship. 
  */ 
  ////end-doc021////

  ////snip022////
  int locFleet,locYear,locAge;
  dvector obsVec(1,4);
  double locObs;
  dvariable locZ;
  dvariable locVar;
  for(int i=1; i<=noobs; ++i){
    obsVec=obs(i);
    locFleet=obsVec(1);
    locYear=obsVec(2);
    locAge=obsVec(3);
    locObs=log(obsVec(4));
    locZ=exp(f(locYear,locAge))+exp(m(locYear,locAge));
    if(locFleet==1){ //    catches predicted 
      pred(i)=f(locYear,locAge)-log(locZ)+log(1.0-exp(-locZ))+n(locYear,locAge);
    }else{           //    survey predicted 
      pred(i)=q(locFleet,locYear,locAge)-locZ*surveyTimes(locFleet)+n(locYear,locAge); 
    }
    locVar=exp(2.0*logSdlogObs(varKey(locFleet,locAge)));
    nll+=nldnorm(locObs,pred(i),locVar);    
    if(sd_phase()){
      residuals(i,1)=locFleet;
      residuals(i,2)=locYear;
      residuals(i,3)=locAge;
      residuals(i,4)=locObs;
      residuals(i,5)=value(pred(i));
      residuals(i,6)=value(sqrt(locVar));
      residuals(i,7)=value((locObs-pred(i))/sqrt(locVar));
    }
  }
  ////end-snip022////
  ////doc022////
  /* 
  For each observation of catch and survey index in the data set, the expected catch and survey index is calculated based on 
  stock sizes, fishing mortalities, and survey catchabilities. From the variance parameters a variance 
  for each observation is assigned, and the normal likelihood contribution (for the log-observation) is added 
  to the negative log likelihood. Finally the residuals are stored for later reporting.     
  */ 
  ////end-doc022////

  ////snip023////
  if(sd_phase()){
    ofstream res("a4a.res");
    res<<"fleet\tyear\tage\tobs\tpred\tsd\tres"<<endl<<residuals<<endl;
    res.close();
  }
  ////end-snip023////
  ////doc023////
  /* 
  The collected residual information is written to its output file. 
  */ 
  ////end-doc023////

  ////snip024////
REPORT_SECTION
  ofstream nout("n.out");
  nout<<exp(n)<<endl;
  nout.close();
  ofstream fout("f.out");
  fout<<exp(f)<<endl;
  fout.close();
  ofstream qout("q.out");
  qout<<q+qdev<<endl;
  qout.close();
  ////end-snip024////
  ////doc024////
  /* 
  The \verb|REPORT\_SECTION| is where custom output is to be written after the likelihood has been optimized w.r.t. the model 
  parameters. For easy later access and plotting the stock sizes, fishing mortalities, and catchabilities are written to 
  separate files.  
  */ 
  ////end-doc024////

  ////snip025////
TOP_OF_MAIN_SECTION
  arrmblsize=2000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(150000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(800000);
  gradient_structure::set_MAX_NVAR_OFFSET(100000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
  ////end-snip025////
  ////doc025////
  /* 
  The \verb|TOP\_OF\_MAIN\_SECTION| is where various buffer sizes are set.  Here we have increased the number of parameter allowed in the model from the default \verb|500| to \verb|5000|.  
  */ 
  ////end-doc025////

