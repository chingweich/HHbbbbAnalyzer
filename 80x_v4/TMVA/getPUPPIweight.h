#ifndef getPUPPIweight_h
#define getPUPPIweight_h
#incldue "TF1.h"
float getPUPPIweight(float puppipt, float puppieta,TF1 *tf1[]){
	float genCorr  = 1.;
	float recoCorr = 1.;
	float totalWeight = 1.;
	
  
  
	genCorr =  tf1[0]->Eval( puppipt );
	if( fabs(puppieta)  <= 1.3 ) recoCorr = tf1[1]->Eval( puppipt );
	else if( fabs(puppieta) > 1.3 ) recoCorr = tf1[2]->Eval( puppipt );
  
	totalWeight = genCorr * recoCorr;

	return totalWeight;
	
}

#endif