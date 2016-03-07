#include <TLegend.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <TH1D.h>
#include <TH2D.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TROOT.h>
#include "TImage.h"
#include "TSystem.h"
#include "TStyle.h"
#include "untuplizer.h"
#include <TClonesArray.h>
#include <fstream>
#include <cmath>
#include <TSystem.h>
#include <string>
#include <sstream>
#include "setNCUStyle.C"
#include "jetEnergyScale.h"
#include <algorithm>
#include "BTagCalibrationStandalone.h"
#include "BTagCalibrationStandalone.cc"

void test(){
	
	
	BTagCalibration calib("csvv1", "CSVv2.csv");
	BTagCalibrationReader reader(&calib,               // calibration instance
                             BTagEntry::OP_LOOSE,  // operating point
                             "comb",               // measurement type
                             "central");           // systematics type

	double jet_scalefactor = reader.eval(BTagEntry::FLAV_B,0.6,400); 
	cout<<jet_scalefactor<<endl;
}