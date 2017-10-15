#include <TMB.hpp>                                

template<class Type>
Type objective_function<Type>::operator() ()
{
	/************DATA_SECTION**************/
	DATA_FACTOR(ageRange); // full age range
	DATA_FACTOR(yearRange); // full year range
	DATA_FACTOR(surveyMinAge); // min age of surveys
	DATA_FACTOR(surveyMaxAge); // max age of surveys
	DATA_VECTOR(surveyTimes); // when survey took place
	DATA_VECTOR(fbarRange); // fbar age range
	DATA_MATRIX(obs); // fleet year age obs weight (TODO this might cause problems integers and doubles together)
	DATA_FACTOR(locFleetVec); //first column of obs must be read in as integer, unless conversion operator exists
	DATA_FACTOR(locYearVec); //second column of obs must be read in as integer, unless conversion operator exists
	DATA_FACTOR(locAgeVec); //third column of obs must be read in as integer, unless conversion operator exists
	DATA_MATRIX(aux); // year, age, m, m.spwn, harvest.spwn, mat*stkwt, stkwt
	DATA_MATRIX(designF); // Fishing mortality model
	DATA_MATRIX(designQ); // Survey catchability model
	DATA_MATRIX(designV); // Variance model
	DATA_MATRIX(designNy1); // Initial age structure model
	DATA_MATRIX(designR); // Internal Recruitment model
	DATA_MATRIX(designRa); // [noExpandedRa,noRapar]
	DATA_MATRIX(designRb); // [noExpandedRb,noRbpar]
	
	DATA_SCALAR(srCV); // what CV is specified
	DATA_SCALAR(spr0);// only used with SV models
	
	//control variables
	DATA_INTEGER(isPlusGrp); // is oldest age a plus group; 0 = NO; 1 = YES
	DATA_INTEGER(Rmodel); //recruitment model; 1=BevHolt, 2=Richer, 3=smooth hockey stick, 4=Geomean, 5=BevHolt with steepness
	
	int noSurveys = surveyMinAge.size();
	int noObs = obs.rows();
	int noAux = aux.rows();
	int minYear = yearRange(0);
	int maxYear = yearRange(1);
	int minAge = ageRange(0);
	int maxAge = ageRange(1);
	
	int noAges = maxAge-minAge+1;
	int noYears= maxYear-minYear+1;
		
	//storage objects
	int idx;
	matrix<Type>  f(noYears, noAges);
	matrix<Type>  m(noYears, noAges);
	matrix<Type>  fspwn(noYears, noAges);
	matrix<Type>  mspwn(noYears, noAges);
	matrix<Type>  matWt(noYears, noAges);
	matrix<Type>  stkWt(noYears, noAges);
	array<Type> q(noSurveys,noYears, noAges);
	vector<Type>  r(noYears); // predicted recruitment
	vector<Type>  ra(noYears); // SRR params
	vector<Type>  rb(noYears); // SRR params
	array<Type> v(noSurveys+1,noYears,noAges);
	matrix<Type>  n(noYears,noAges);
	vector<Type>  pred(noObs);
	vector<Type> ssb(noYears);
	Type ssbmaxYear;
	Type predLogR;
	Type varLogR;
	Type h;
	Type tmpv;
	
	//reorganize data
	idx = 0;
	for (int y=0; y<=noYears-1; y++) {
		for (int a=0; a<=noAges-1; a++) {
			m(y,a) = aux(idx,2);
			mspwn(y,a) = aux(idx,3);
			fspwn(y,a) = aux(idx,4);
			matWt(y,a) = aux(idx,5);
			stkWt(y,a) = aux(idx,6);
			idx++;
		}
	}
	/************PARAMETER_SECTION**************/
	
	PARAMETER_VECTOR(fpar);
	PARAMETER_VECTOR(qpar);
	PARAMETER_VECTOR(vpar);
	PARAMETER_VECTOR(ny1par);
	PARAMETER_VECTOR(rpar);
	PARAMETER_VECTOR(rapar);
	PARAMETER_VECTOR(rbpar);

	/************PROCEDURE_SECTION**************/
	Type nll;
	nll=0;
	
	// fishing mortality model
	vector<Type> expandedF = designF*fpar; //+fdev
	idx = 0;
	for (int y=0; y<=noYears-1; y++) {
		for (int a=0; a<=noAges-1; a++) {
			f(y,a) = expandedF(idx++);
		}
	}
		
	// survey catchability model
	vector<Type> expandedQ = designQ*qpar; //+qdev
	idx = 0;
	for(int ff = 0; ff <= noSurveys-1; ff++){ 
		for(int y = 0; y <= noYears-1; y++){
			for(int a = 0; a <= noAges-1; a++){
				q(ff,y,a) = expandedQ(idx++);
			}
		}
	}

	// variance model
	vector<Type> expandedV = designV*vpar;
	idx = 0;
	for(int ff = 0; ff <= noSurveys; ff++){ 
		for(int y = 0; y <= noYears-1; y++){
			for(int a = 0; a <= noAges-1; a++){
				v(ff,y,a) = expandedV(idx++);
			}
		}
	}
	// initial age structure model
	vector<Type> expandedNy1 = designNy1*ny1par;
	idx = 0;
	for(int a = 1; a <= noAges-1; a++){// Start at minAge + 1
		n(0,a) = expandedNy1(idx++);
	}

	// internal r model. 
	vector<Type> expandedR = designR*rpar;
	idx = 0;
	for(int y = 0; y <= noYears-1; y++){//This loop is unnecessary.
		r(y) = expandedR(idx++);
	}
	// full population structure
	n.col(0)=r;
	for(int a=1; a<=noAges-1; a++){
		for(int y=1; y<=noYears-1; y++){
			n(y,a)=n(y-1,a-1)-exp(f(y-1,a-1))-exp(m(y-1,a-1));
			if((a==noAges-1) && (isPlusGrp > 0.5)){
				n(y,a)=log(exp(n(y,a))+exp(n(y-1,a)-exp(f(y-1,a))-exp(m(y-1,a))));
			}
		}
	}
	// fbar and ssb
	for(int y=0; y<=noYears-1; y++){
        ssb(y) = ((n.row(y).array() - f.row(y).array().exp() * fspwn.row(y).array() - m.row(y).array().exp() * mspwn.row(y).array() ).exp() * matWt.row(y).array() ).sum();
	}	
	ssbmaxYear = ssb(noYears-1);

	// main likelihood
	int locFleet;
	int locYear;
	int locAge;
	int minSurveyAge;
	int maxSurveyAge;
	Type locObs;
	Type locZ;
	Type locVar;
	for (int i=0; i<=noObs-1; i++) {
		locFleet = locFleetVec(i);
		locYear  = locYearVec(i)-minYear;
		locAge   = locAgeVec(i)-minAge;
		locObs   = obs(i,3);

		//here we split - if locAge == -1 then we have a biomass index and use add to a different likelihood component.
		if (locAge >= 0) { // standard observation
			locZ = exp(f(locYear,locAge)) + exp(m(locYear,locAge));
			if (locFleet==1) { //    catches predicted 
				pred(i) = f(locYear,locAge)-log(locZ)+log(Type(1.0)-exp(-locZ))+n(locYear,locAge);
			} 
			else {          //    survey predicted. Fleet >=2
				pred(i) = q(locFleet-2,locYear,locAge) - locZ * surveyTimes(locFleet-2) + n(locYear,locAge); 
			}
			locVar = exp(Type(2.0) * v(locFleet-1, locYear, locAge));
			nll -= obs(i,4) * dnorm(locObs, pred(i), sqrt(locVar), true); // or do we multiply the variance directly...    
		}
		 else { // an observation of biomass
		//TODO, be careful in this case. what is Fleet? Is it always >=2?
			pred(i) = Type(0); // not sure if need to but best to be safe
			for(int a = surveyMinAge(locFleet-2)-minAge; a<=surveyMaxAge(locFleet-2)-minAge; a++) {
				locZ     = exp(f(locYear,a)) + exp(m(locYear,a));
				pred(i) += exp(q(locFleet-2, locYear, a)) * stkWt(locYear, a) * exp(n(locYear,a) - surveyTimes(locFleet-2) * locZ);
			}
			locVar = exp(Type(2.0) * v(locFleet, locYear, minAge)); // note variance are stored in the minimum age column
			nll -= obs(i,4) * dnorm(locObs, log(pred(i)), sqrt(locVar), true); // or do we multiply the variance directly...    
		}
	}

    if (srCV >= 0) {
	// stock recruit model
	vector<Type> expandedRa = designRa*rapar; //+ radev;
	vector<Type> expandedRb = designRb*rbpar;
	idx = 0;
	for (int y = minAge; y <= noYears-1; y++) {
		ra(y) = expandedRa(idx); // a and b are correlated though...
		rb(y) = expandedRb(idx);
		idx++;
	}
	varLogR = log(pow(srCV,Type(2.0))+Type(1.0));//const for any Rmodel. Only do once.
	//there could be a probelm because r(y) is predicted and not observed
	if (Rmodel == 1) { // beverton holt
		for(int y=minAge; y<=noYears-1; y++){
			predLogR = ra(y) + log(ssb(y-minAge)) - log(exp(rb(y)) + ssb(y-minAge));
			nll += -dnorm(r(y), predLogR, sqrt(varLogR), true);    
		}
	}
	if (Rmodel == 2) { // ricker
		for(int y=minAge; y<=noYears-1; y++){
			predLogR = ra(y) + log(ssb(y-minAge)) - exp(rb(y)) * ssb(y-minAge);
			nll += -dnorm(r(y), predLogR, sqrt(varLogR), true);    
		}
	}
	if (Rmodel == 3) { // smooth hockey stick (Mesnil and Rochet, gamma = 0.1)
		for(int y=minAge; y<=noYears-1; y++){
			predLogR = ra(y) + log(ssb(y-minAge) + sqrt(exp(Type(2.0)*rb(y)) + Type(0.0025)) - sqrt(pow(ssb(y-minAge) - exp(2.0*rb(y)), Type(2.0)) + Type(0.0025)));
			nll += -dnorm(r(y), predLogR, sqrt(varLogR), true);   
		}
	}
	if (Rmodel == 4) { // geomean
		for(int y=1; y<=noYears-1; y++){
			predLogR = ra(y);
			nll += -dnorm(r(y), predLogR, sqrt(varLogR), true);    
		}
	}
    }


    
//	if (Rmodel == 5) { // bevholt with steepness: ra is a transform of h; rb is a transform of v
//		for(int y=minAge; y<=noYears-1; y++){
//			h = exp(ra(y)) / (Type(1.0) + exp(ra(y))) * Type(0.8) + Type(0.2);
//			tmpv = exp(rb(y));
//			predLogR =  log(Type(6.0) * h * tmpv * ssb(y-minAge)) - log( spr0 * ((h + Type(1.0))*v + (Type(5.0)*h - Type(1.0))*ssb(y-minAge)) ); // spr0 is provided by user
//			nll += -dnorm(r(y), predLogR,  sqrt(varLogR), true);   
//		}
//	}
	
	/************REPORT_SECTION**************/
	
	array<Type> expn(noYears, noAges);//
	expn=n.array().exp();
	array<Type> expf(noYears, noAges);
	expf=f.array().exp();
	ADREPORT(expn);
	ADREPORT(expf);
	ADREPORT(q);
	ADREPORT(v);
	ADREPORT(ssbmaxYear);
	
	return nll;
//	return res;
  
}
