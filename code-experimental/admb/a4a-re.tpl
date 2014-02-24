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

GLOBALS_SECTION 
  #include <time.h>
  time_t StartTime;   //Time variables, for use in timing loop
  #include <df1b2fun.h>
  #include "nLogNormal.h"
  ofstream clogf("program.log");
  #define TRACE(object) clogf<<"line "<<__LINE__<<", file "<<__FILE__<<", "<<#object" =\n"<<object<<endl<<endl; 

  dvariable nldnorm(const double& obs, const dvariable& mu, const dvariable& var){
    return 0.5*(log(2.0*M_PI*var)+square(obs-mu)/var);
  }

  df1b2variable nldnorm(const double& obs, const df1b2variable& mu, const df1b2variable& var){
    return 0.5*(log(2.0*M_PI*var)+square(obs-mu)/var);
  }

  dvariable nldnorm(const dvariable& obs, const dvariable& mu, const dvariable& var){
    return 0.5*(log(2.0*M_PI*var)+square(obs-mu)/var);
  }

  df1b2variable nldnorm(const df1b2variable& obs, const df1b2variable& mu, const df1b2variable& var){
    return 0.5*(log(2.0*M_PI*var)+square(obs-mu)/var);
  }

  ////////////////////////////////////////////////////////////////////////////////////

DATA_SECTION
  !! time(&StartTime);
  init_int noobs
  !! TRACE(noobs)
  init_matrix obs(1,noobs,1,4) // f y a o 
  !! TRACE(obs)
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

  !! ad_comm::change_datafile_name("model.cfg");
  init_imatrix keyQ(2,maxFleet,minAge,maxAge)
  !! TRACE(keyQ)
  int noQpar
  !! noQpar=max(keyQ);
  !! TRACE(noQpar)
  init_vector surveyTimes(2,maxFleet)
  !! TRACE(surveyTimes)
  init_vector fbarRange(1,2)
  !! TRACE(fbarRange)

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

  matrix residuals(1,noobs,1,7)
 
PARAMETER_SECTION

  init_vector ry1(minAge+1,maxAge)
  init_vector q(1,noQpar)
  init_vector fpar(1,noFpar)
  init_vector logSdlogObs(1,maxFleet)

  init_number rec_loga
  init_number rec_logb
  init_number logSdLogR

  random_effects_vector r(minYear,maxYear)
  objective_function_value jnll;
 
  vector expandedF(1,noExpandedF)
  matrix f(minYear,maxYear,minAge,maxAge)
  matrix n(minYear,maxYear,minAge,maxAge)
  vector pred(1,noobs)  
  sdreport_vector ssb(minYear,maxYear)
  sdreport_vector fbar(minYear,maxYear)

PRELIMINARY_CALCS_SECTION
  rec_loga=1; 
  rec_logb=-12;

PROCEDURE_SECTION
  time_t currentTime;
  time(& currentTime);
  if(difftime(currentTime,StartTime)>1800){ // Terminate after 30 minutes 
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
  jnll=0.0;

  expandedF=designF*fpar;
 
  int idx=0;
  for(int a=minAge; a<=maxAge; ++a){
    for(int y=minYear; y<=maxYear; ++y){
      f(y,a)=expandedF(++idx);
    }
  }
 
  //n.colfill(minAge,r);
  for(int y=minYear; y<=maxYear; ++y){n(y,minAge)=r(y);}

  for(int a=minAge+1; a<=maxAge; ++a){
    n(minYear,a)=ry1(a);
    for(int y=minYear+1; y<=maxYear; ++y){
      n(y,a)=n(y-1,a-1)-exp(f(y-1,a-1))-exp(m(y-1,a-1));
      if(a==maxAge){
        n(y,a)=log(exp(n(y,a))+exp(n(y-1,a)-exp(f(y-1,a))-exp(m(y-1,a))));
      }
    }
  }
  
  for(int y=minYear; y<=maxYear; ++y){
    ssb(y)=sum(elem_prod(elem_prod(exp(n(y)),sw(y)),mat(y)));
    fbar(y)=0.0;
    for(int a=fbarRange(1); a<=fbarRange(2); ++a){
      fbar(y)+=exp(f(y,a));
    }
    fbar(y)/=fbarRange(2)-fbarRange(1)+1;
  }  

  //S-R stuff 
  dvariable predLogR; 
  dvariable varLogR; 
  for(int y=minYear+1; y<=maxYear; ++y){
    predLogR=rec_loga+log(ssb(y-1))-log(1.0+exp(rec_logb)*ssb(y-1));
    varLogR=exp(2.0*logSdLogR);
    jnll+=nldnorm(r(y),predLogR,varLogR);    
  }

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
    if(locFleet==1){ //catches
      pred(i)=f(locYear,locAge)-log(locZ)+log(1.0-exp(-locZ))+n(locYear,locAge);
    }else{ //survey
      pred(i)=q(keyQ(locFleet,locAge))-locZ*surveyTimes(locFleet)+n(locYear,locAge); 
    }
    locVar=exp(2.0*logSdlogObs(locFleet));
    jnll+=nldnorm(locObs,pred(i),locVar);    
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
  if(sd_phase()){
    ofstream res("a4a.res");
    res<<"fleet\tyear\tage\tobs\tpred\tsd\tres"<<endl<<residuals<<endl;
    res.close();
  }

REPORT_SECTION
  ofstream nout("n.out");
  nout<<exp(n)<<endl;
  nout.close();
  ofstream fout("f.out");
  fout<<exp(f)<<endl;
  fout.close();

TOP_OF_MAIN_SECTION
  arrmblsize=2000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(150000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(800000);
  gradient_structure::set_MAX_NVAR_OFFSET(100000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);

