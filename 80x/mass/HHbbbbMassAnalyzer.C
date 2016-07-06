//Histogram 
#include <TH1.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TProfile.h>
//vector ,string ,stream
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

//root feature
#include <TLegend.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include "TSystem.h"
#include "TStyle.h"
#include <TClonesArray.h>
#include <TSystem.h>

//math 
#include <cmath>
#include <algorithm>

//other including
//#include "setNCUStyle.C"
#include "../untuplizer.h"
//#include "jetEnergyScale.h"


#include "HHbbbbMassAnalyzerBase.C"

using namespace std;

void HHbbbbMassAnalyzer(){

	HHbbbbMassAnalyzerBase(1,2,"Bulk1200.root","Bulk1200",1);
	HHbbbbMassAnalyzerBase(1,2,"Bulk1400.root","Bulk1400");
	HHbbbbMassAnalyzerBase(1,2,"Bulk1600.root","Bulk1600");
	HHbbbbMassAnalyzerBase(1,2,"Bulk1800.root","Bulk1800");
	HHbbbbMassAnalyzerBase(1,2,"Bulk2500.root","Bulk2500");
	HHbbbbMassAnalyzerBase(1,2,"Bulk3000.root","Bulk3000");
	HHbbbbMassAnalyzerBase(1,2,"Bulk4000.root","Bulk4000");
	HHbbbbMassAnalyzerBase(1,2,"Bulk4500.root","Bulk4500",2);
	
	
	HHbbbbMassAnalyzerBase(1,2,"BulkWW600.root","BulkWW600",1);
	HHbbbbMassAnalyzerBase(1,2,"BulkWW800.root","BulkWW800");
	HHbbbbMassAnalyzerBase(1,2,"BulkWW1000.root","BulkWW1000");
	HHbbbbMassAnalyzerBase(1,2,"BulkWW1400.root","BulkWW1400");
	HHbbbbMassAnalyzerBase(1,2,"BulkWW1600.root","BulkWW1600");
	HHbbbbMassAnalyzerBase(1,2,"BulkWW1800.root","BulkWW1800");
	HHbbbbMassAnalyzerBase(1,2,"BulkWW2000.root","BulkWW2000");
	HHbbbbMassAnalyzerBase(1,2,"BulkWW2500.root","BulkWW2500");
	HHbbbbMassAnalyzerBase(1,2,"BulkWW3000.root","BulkWW3000");
	HHbbbbMassAnalyzerBase(1,2,"BulkWW3500.root","BulkWW3500");
	HHbbbbMassAnalyzerBase(1,2,"BulkWW4000.root","BulkWW4000");
	HHbbbbMassAnalyzerBase(1,2,"BulkWW4500.root","BulkWW4500",2);
	
}