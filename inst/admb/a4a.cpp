  #include <time.h>
  time_t StartTime;   //Time variables, for use in timing loop
  #include <df1b2fun.h>
  #include <string>
  using namespace std;
  //#include "nLogNormal.h"
  //ofstream clogf("program.log");
  //#define TRACE(object) clogf<<"line "<<__LINE__<<", file "<<__FILE__<<", "\
  //                           <<#object" =\n"<<object<<endl<<endl;
  /* differentiable function for the likelihood */ 
  dvariable nldnorm(const double& obs, const dvariable& mu, 
                    const dvariable& var){
    return 0.5*(log(2.0*M_PI*var)+square(obs-mu)/var);
  }
  /* differentiable function for the likelihood */ 
  dvariable nldnorm(const dvariable& obs, const dvariable& mu, 
                    const dvariable& var){
    return 0.5*(log(2.0*M_PI*var)+square(obs-mu)/var);
  }
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <a4a.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
 time(&StartTime);
  ageRange.allocate(1,2,"ageRange");
  yearRange.allocate(1,2,"yearRange");
  nsurveys.allocate("nsurveys");
  surveyMinAge.allocate(1,nsurveys,"surveyMinAge");
  surveyMaxAge.allocate(1,nsurveys,"surveyMaxAge");
  surveyTimes.allocate(1,nsurveys,"surveyTimes");
  fbarRange.allocate(1,2,"fbarRange");
  isPlusGrp.allocate("isPlusGrp");
  noobs.allocate("noobs");
  obs.allocate(1,noobs,1,5,"obs");
  noaux.allocate("noaux");
  aux.allocate(1,noaux,1,7,"aux");
 minYear = yearRange(1);
 maxYear = yearRange(2);
 minAge = ageRange(1);
 maxAge = ageRange(2);
 ad_comm::change_datafile_name("fmodel.cfg");
  noFpar.allocate("noFpar");
  noExpandedF.allocate("noExpandedF");
  designF.allocate(1,noExpandedF,1,noFpar,"designF");
 ad_comm::change_datafile_name("qmodel.cfg");
  noQpar.allocate("noQpar");
  noExpandedQ.allocate("noExpandedQ");
  designQ.allocate(1,noExpandedQ,1,noQpar,"designQ");
 ad_comm::change_datafile_name("vmodel.cfg");
  noVpar.allocate("noVpar");
  noExpandedV.allocate("noExpandedV");
  designV.allocate(1,noExpandedV,1,noVpar,"designV");
 ad_comm::change_datafile_name("ny1model.cfg");
  noNy1par.allocate("noNy1par");
  noExpandedNy1.allocate("noExpandedNy1");
  designNy1.allocate(1,noExpandedNy1,1,noNy1par,"designNy1");
 ad_comm::change_datafile_name("rmodel.cfg");
  noRpar.allocate("noRpar");
  noExpandedR.allocate("noExpandedR");
  designR.allocate(1,noExpandedR,1,noRpar,"designR");
 ad_comm::change_datafile_name("srrmodel.cfg");
  Rmodel.allocate("Rmodel");
  srCV.allocate("srCV");
 if(srCV>0){SRaphase = 2;}else{SRaphase = -1;}
 if(srCV < 0 | Rmodel == 3){SRbphase = -1;}else{SRbphase = 2;}
 if(srCV < 0 | Rmodel == 4){SRbphase = -1;}else{SRbphase = 2;}
  spr0.allocate("spr0");
  noRapar.allocate("noRapar");
  noExpandedRa.allocate("noExpandedRa");
  designRa.allocate(1,noExpandedRa,1,noRapar,"designRa");
  noRbpar.allocate("noRbpar");
  noExpandedRb.allocate("noExpandedRb");
  designRb.allocate(1,noExpandedRb,1,noRbpar,"designRb");
  pad_NMCMCreport = new ofstream("NMCMCreport.csv",ios::trunc);;
  pad_FMCMCreport = new ofstream("FMCMCreport.csv",ios::trunc);;
  pad_QMCMCreport = new ofstream("QMCMCreport.csv",ios::trunc);;
  pad_VMCMCreport = new ofstream("VMCMCreport.csv",ios::trunc);;
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  fpar.allocate(1,noFpar,3,"fpar");
  qpar.allocate(1,noQpar,2,"qpar");
  vpar.allocate(1,noVpar,1,"vpar");
  ny1par.allocate(1,noNy1par,4,"ny1par");
  rpar.allocate(1,noRpar,1,"rpar");
  rapar.allocate(1,noRapar,SRaphase,"rapar");
  rbpar.allocate(1,noRbpar,SRbphase,"rbpar");
  nll.allocate("nll");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  expandedF.allocate(1,noExpandedF,"expandedF");
  #ifndef NO_AD_INITIALIZE
    expandedF.initialize();
  #endif
  expandedQ.allocate(1,noExpandedQ,"expandedQ");
  #ifndef NO_AD_INITIALIZE
    expandedQ.initialize();
  #endif
  expandedV.allocate(1,noExpandedV,"expandedV");
  #ifndef NO_AD_INITIALIZE
    expandedV.initialize();
  #endif
  expandedNy1.allocate(1,noExpandedNy1,"expandedNy1");
  #ifndef NO_AD_INITIALIZE
    expandedNy1.initialize();
  #endif
  expandedR.allocate(1,noExpandedR,"expandedR");
  #ifndef NO_AD_INITIALIZE
    expandedR.initialize();
  #endif
  expandedRa.allocate(1,noExpandedRa,"expandedRa");
  #ifndef NO_AD_INITIALIZE
    expandedRa.initialize();
  #endif
  expandedRb.allocate(1,noExpandedRb,"expandedRb");
  #ifndef NO_AD_INITIALIZE
    expandedRb.initialize();
  #endif
  f.allocate(minYear,maxYear,minAge,maxAge,"f");
  #ifndef NO_AD_INITIALIZE
    f.initialize();
  #endif
  m.allocate(minYear,maxYear,minAge,maxAge,"m");
  #ifndef NO_AD_INITIALIZE
    m.initialize();
  #endif
  fspwn.allocate(minYear,maxYear,minAge,maxAge,"fspwn");
  #ifndef NO_AD_INITIALIZE
    fspwn.initialize();
  #endif
  mspwn.allocate(minYear,maxYear,minAge,maxAge,"mspwn");
  #ifndef NO_AD_INITIALIZE
    mspwn.initialize();
  #endif
  matWt.allocate(minYear,maxYear,minAge,maxAge,"matWt");
  #ifndef NO_AD_INITIALIZE
    matWt.initialize();
  #endif
  stkWt.allocate(minYear,maxYear,minAge,maxAge,"stkWt");
  #ifndef NO_AD_INITIALIZE
    stkWt.initialize();
  #endif
  q.allocate(1,nsurveys,minYear,maxYear,minAge,maxAge,"q");
  #ifndef NO_AD_INITIALIZE
    q.initialize();
  #endif
  r.allocate(minYear,maxYear,"r");
  #ifndef NO_AD_INITIALIZE
    r.initialize();
  #endif
  ra.allocate(minYear,maxYear,"ra");
  #ifndef NO_AD_INITIALIZE
    ra.initialize();
  #endif
  rb.allocate(minYear,maxYear,"rb");
  #ifndef NO_AD_INITIALIZE
    rb.initialize();
  #endif
  v.allocate(1,nsurveys+1,minYear,maxYear,minAge,maxAge,"v");
  #ifndef NO_AD_INITIALIZE
    v.initialize();
  #endif
  n.allocate(minYear,maxYear,minAge,maxAge,"n");
  #ifndef NO_AD_INITIALIZE
    n.initialize();
  #endif
  pred.allocate(1,noobs,"pred");
  #ifndef NO_AD_INITIALIZE
    pred.initialize();
  #endif
  ssb.allocate(minYear,maxYear,"ssb");
  #ifndef NO_AD_INITIALIZE
    ssb.initialize();
  #endif
  ssbmaxYear.allocate("ssbmaxYear");
}

void model_parameters::initializationfunction(void)
{
  fpar.set_initial_value(0);
  qpar.set_initial_value(0);
  vpar.set_initial_value(0);
  ny1par.set_initial_value(0);
  rpar.set_initial_value(0);
  rapar.set_initial_value(0);
  rbpar.set_initial_value(0);
}

void model_parameters::preliminary_calculations(void)
{

  admaster_slave_variable_interface(*this);
  
  //if(srCV>0){
  //  logSdLogR=log(srCV);
  //}
  //fpar=log(0.1)*max(designF);
  //qpar=log(1E-5)*max(designF);
  //
  // auxilliary data
  //
  idx = 0;
  for (int y=minYear; y<=maxYear; ++y) {
    for (int a=minAge; a<=maxAge; ++a) {
      idx = idx + 1;
      m(y,a) = aux(idx,3);
      mspwn(y,a) = aux(idx,4);
      fspwn(y,a) = aux(idx,5);
      matWt(y,a) = aux(idx,6);
      stkWt(y,a) = aux(idx,7);
    }
  }
}

void model_parameters::userfunction(void)
{
  nll =0.0;
  ofstream& NMCMCreport= *pad_NMCMCreport;
  ofstream& FMCMCreport= *pad_FMCMCreport;
  ofstream& QMCMCreport= *pad_QMCMCreport;
  ofstream& VMCMCreport= *pad_VMCMCreport;
  time_t currentTime;
  time(& currentTime);
  if(difftime(currentTime,StartTime)>3600){ // Terminate after 60 minutes 
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
  //
  // Fishing mortality model
  //
  expandedF = designF*fpar; //+ fdev; 
  idx = 0;
  for (int y=minYear; y<=maxYear; ++y) {
    for (int a=minAge; a<=maxAge; ++a) {
      f(y,a) = expandedF(++idx);
    }
  }
  //
  // survey catchability model
  //
  expandedQ = designQ*qpar; //+ qdev;
  idx = 0;
  for(int ff = 1; ff <= nsurveys; ++ff){ 
    for(int y = minYear; y <= maxYear; ++y){
      for(int a = minAge; a <= maxAge; ++a){
        q(ff,y,a) = expandedQ(++idx);
      }
    }
  }    
  //
  // variance model
  //
  expandedV = designV*vpar;
  idx = 0;
  for(int ff = 1; ff <= nsurveys + 1; ++ff){ 
    for(int y = minYear; y <= maxYear; ++y){
      for(int a = minAge; a <= maxAge; ++a){
        v(ff,y,a) = expandedV(++idx);
      }
    }
  }    
  //
  // initial age structure model
  //
  expandedNy1 = designNy1*ny1par;
  idx = 0;
  for(int a = minAge+1; a <= maxAge; ++a){
    n(minYear,a) = expandedNy1(++idx);
  }
  //
  // internal r model
  //
  expandedR = designR*rpar;
  idx = 0;
  for(int y = minYear; y <= maxYear; ++y){
    r(y) = expandedR(++idx);
  }
  //
  // the full population structure
  //
  n.colfill(minAge,r);
  for(int a=minAge+1; a<=maxAge; ++a){
    for(int y=minYear+1; y<=maxYear; ++y){
      n(y,a)=n(y-1,a-1)-mfexp(f(y-1,a-1))-mfexp(m(y-1,a-1));
      if((a==maxAge) && (isPlusGrp > 0.5)){
        n(y,a)=log(mfexp(n(y,a))+mfexp(n(y-1,a)-mfexp(f(y-1,a))-mfexp(m(y-1,a))));
      }
    }
  }
  //
  // fbar and ssb
  //
  for(int y=minYear; y<=maxYear; ++y){
	ssb(y) = sum(elem_prod(mfexp(n(y)-mfexp(f(y))*fspwn(y)-mfexp(m(y))*mspwn(y)), matWt(y))); 
  }  
  ssbmaxYear = ssb(maxYear);
  //
  // The main likelihood
  //
  int locFleet,locYear,locAge;
  dvector obsVec(1,5);
	int minSurveyAge, maxSurveyAge;
  double locObs;
  dvariable locZ;
  dvariable locVar;
  for (int i=1; i<=noobs; ++i) {
    obsVec   = obs(i);
    locFleet = obsVec(1);
    locYear  = obsVec(2);
    locAge   = obsVec(3);
    locObs   = obsVec(4); 
    // here we split - if locAge == -1 then we have a biomass index and use add to a different likelihood component.
    if (locAge >= 0) { // standard observation
      locZ     = mfexp(f(locYear,locAge)) + mfexp(m(locYear,locAge));
      if (locFleet==1) { //    catches predicted 
        pred(i) = f(locYear,locAge)-log(locZ)+log(1.0-mfexp(-locZ))+n(locYear,locAge);
      } else {          //    survey predicted 
        pred(i) = q(locFleet-1,locYear,locAge) - locZ * surveyTimes(locFleet-1) + n(locYear,locAge); 
      }
      locVar = mfexp(2.0 * v(locFleet, locYear, locAge));
      nll += obsVec(5) * nldnorm(locObs, pred(i), locVar); // or do we multiply the variance directly...    
    } else { // an observation of biomass
      pred(i) = 0; // not sure i need to but best to be safe
      for(int a=surveyMinAge(locFleet-1); a<=surveyMaxAge(locFleet-1); ++a) {
        locZ     = mfexp(f(locYear,a)) + mfexp(m(locYear,a));
        pred(i) += mfexp(q(locFleet-1, locYear, a)) * stkWt(locYear, a) * mfexp(n(locYear,a) - surveyTimes(locFleet-1) * locZ);
      }
      locVar = mfexp(2.0 * v(locFleet, locYear, minAge)); // note variance are stored in the minimum age column
      nll += obsVec(5) * nldnorm(locObs, log(pred(i)), locVar); // or do we multiply the variance directly...    
    }
  }
  //
  // stock recruit model
  //
  if(SRaphase > 0) { // then include a SRR model 
    //
    // recruitment model
    //
    expandedRa = designRa*rapar; //+ radev;
    expandedRb = designRb*rbpar;
    idx = 0;
    for (int y = minYear; y <= maxYear; ++y) {
      idx = idx+1;
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
      for(int y=minYear+minAge; y<=maxYear; ++y){
        predLogR = ra(y) + log(ssb(y-minAge)) - log(mfexp(rb(y)) + ssb(y-minAge));
        varLogR = log(pow(srCV,2)+1);
        nll += nldnorm(r(y), predLogR, varLogR);    
      }
    }
    if (Rmodel == 2) { // ricker
      for(int y=minYear+minAge; y<=maxYear; ++y){
        predLogR = ra(y) + log(ssb(y-minAge)) - mfexp(rb(y)) * ssb(y-minAge);
        varLogR = log(pow(srCV,2)+1);
        nll += nldnorm(r(y), predLogR, varLogR);    
      }
    }
    if (Rmodel == 3) { // smooth hockey stick (Mesnil and Rochet, gamma = 0.1)
      for(int y=minYear+minAge; y<=maxYear; ++y){
        predLogR = ra(y) + log(ssb(y-minAge) + sqrt(mfexp(2.0*rb(y)) + 0.0025) - sqrt(pow(ssb(y-minAge) - mfexp(2.0*rb(y)), 2.0) + 0.0025));
        varLogR = log(pow(srCV,2)+1);
        nll += nldnorm(r(y), predLogR, varLogR);    
      }
    }
    if (Rmodel == 4) { // geomean
      for(int y=minYear+1; y<=maxYear; ++y){
        predLogR = ra(y);
        varLogR = log(pow(srCV,2)+1);
        nll += nldnorm(r(y), predLogR, varLogR);    
      }
    }
    if (Rmodel == 5) { // bevholt with steepness: ra is a transform of h; rb is a transform of v
      for(int y=minYear+minAge; y<=maxYear; ++y){
        h = mfexp(ra(y)) / (1 + mfexp(ra(y))) * 0.8 + 0.2;
        v = mfexp(rb(y));
        predLogR =  log(6 * h * v * ssb(y-minAge)) - log( spr0 * ((h + 1)*v + (5*h - 1)*ssb(y-minAge)) ); // spr0 is provided by user
        varLogR = log(pow(srCV,2)+1);
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
}

void model_parameters::set_runtime(void)
{
  dvector temp("{1E-1,1E-2,1E-3,1E-12}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
  dvector temp1("{10,20,30,10000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
}

void model_parameters::report()
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
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
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{
  delete pad_NMCMCreport;
  pad_NMCMCreport = NULL;
  delete pad_FMCMCreport;
  pad_FMCMCreport = NULL;
  delete pad_QMCMCreport;
  pad_QMCMCreport = NULL;
  delete pad_VMCMCreport;
  pad_VMCMCreport = NULL;
}

void model_parameters::final_calcs(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  arrmblsize=2000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(150000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(800000);
  gradient_structure::set_MAX_NVAR_OFFSET(100000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
