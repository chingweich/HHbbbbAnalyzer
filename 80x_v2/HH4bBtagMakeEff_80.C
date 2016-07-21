//Histogram 
#include <TH1.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TGraphErrors.h>

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
#include "untuplizer.h"
//#include "jetEnergyScale.h"

#include "standalone_LumiReWeighting.cc"
#include "standalone_LumiReWeighting.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"
#include "HH4bBtagMakeEffBase_80.C"
//#include "HHbbbbBtagEffBase_76.C"
//#include "HHbbbbAnalyzerBaseC_76.C"


using namespace std;

void HH4bBtagMakeEff_80(int a){

	string st1[40]={
		/*0-11*/
		"/data7/chchen/80x/BulkGravTohhTohbbhbb_narrow_M/BulkGravTohhTohbbhbb_narrow_M-1000_13TeV-madgraph.root",
		"/data7/chchen/80x/BulkGravTohhTohbbhbb_narrow_M/BulkGravTohhTohbbhbb_narrow_M-1200_13TeV-madgraph.root",
		"/data7/chchen/80x/BulkGravTohhTohbbhbb_narrow_M/BulkGravTohhTohbbhbb_narrow_M-1400_13TeV-madgraph.root",
		"/data7/chchen/80x/BulkGravTohhTohbbhbb_narrow_M/BulkGravTohhTohbbhbb_narrow_M-1600_13TeV-madgraph.root",
		"/data7/chchen/80x/BulkGravTohhTohbbhbb_narrow_M/BulkGravTohhTohbbhbb_narrow_M-1800_13TeV-madgraph.root",
		"/data7/chchen/80x/BulkGravTohhTohbbhbb_narrow_M/BulkGravTohhTohbbhbb_narrow_M-2000_13TeV-madgraph.root",
		"/data7/chchen/80x/BulkGravTohhTohbbhbb_narrow_M/BulkGravTohhTohbbhbb_narrow_M-2500_13TeV-madgraph.root",
		"/data7/chchen/80x/BulkGravTohhTohbbhbb_narrow_M/BulkGravTohhTohbbhbb_narrow_M-3000_13TeV-madgraph.root",
		"/data7/chchen/80x/BulkGravTohhTohbbhbb_narrow_M/BulkGravTohhTohbbhbb_narrow_M-3500_13TeV-madgraph.root",
		"/data7/chchen/80x/BulkGravTohhTohbbhbb_narrow_M/BulkGravTohhTohbbhbb_narrow_M-4000_13TeV-madgraph.root",
		"/data7/chchen/80x/BulkGravTohhTohbbhbb_narrow_M/BulkGravTohhTohbbhbb_narrow_M-4500_13TeV-madgraph.root",
		/*12-21*/
		"/data7/chchen/0606/RadionTohhTohbbhbb_narrow_M/RadionTohhTohbbhbb_narrow_M-1000_13TeV-madgraph.root",
		"/data7/chchen/0606/RadionTohhTohbbhbb_narrow_M/RadionTohhTohbbhbb_narrow_M-1200_13TeV-madgraph.root",
		"/data7/chchen/0606/RadionTohhTohbbhbb_narrow_M/RadionTohhTohbbhbb_narrow_M-1400_13TeV-madgraph.root",
		"/data7/chchen/0606/RadionTohhTohbbhbb_narrow_M/RadionTohhTohbbhbb_narrow_M-1600_13TeV-madgraph.root",
		"/data7/chchen/0606/RadionTohhTohbbhbb_narrow_M/RadionTohhTohbbhbb_narrow_M-1800_13TeV-madgraph.root",
		"/data7/chchen/0606/RadionTohhTohbbhbb_narrow_M/RadionTohhTohbbhbb_narrow_M-2000_13TeV-madgraph.root",
		"/data7/chchen/0606/RadionTohhTohbbhbb_narrow_M/RadionTohhTohbbhbb_narrow_M-2500_13TeV-madgraph.root",
		"/data7/chchen/0606/RadionTohhTohbbhbb_narrow_M/RadionTohhTohbbhbb_narrow_M-3000_13TeV-madgraph.root",
		"/data7/chchen/0606/RadionTohhTohbbhbb_narrow_M/RadionTohhTohbbhbb_narrow_M-3500_13TeV-madgraph.root",
		"/data7/chchen/0606/RadionTohhTohbbhbb_narrow_M/RadionTohhTohbbhbb_narrow_M-4000_13TeV-madgraph.root",
		"/data7/chchen/0606/RadionTohhTohbbhbb_narrow_M/RadionTohhTohbbhbb_narrow_M-4500_13TeV-madgraph.root",
		/*22-32*/
		"/data7/chchen/80x/QCD/700/NCUGlobalTuples_",
		"/data7/chchen/80x/QCD/700_1/NCUGlobalTuples_",
		"/data7/chchen/80x/QCD/1000/NCUGlobalTuples_",
		"/data7/chchen/80x/QCD/1000_1/NCUGlobalTuples_",
		"/data7/chchen/80x/QCD/1500/NCUGlobalTuples_",
		"/data7/chchen/80x/QCD/1500_1/NCUGlobalTuples_",
		"/data7/chchen/80x/QCD/2000/NCUGlobalTuples_",
		"/data7/chchen/80x/QCD/2000_1/NCUGlobalTuples_",
		
		
	};
	string  fileName[40]={
	"B1000","B1200","B1400","B1600","B1800","B2000","B2500","B3000","B3500","B4000","B4500",
	"R1000","R1200","R1400","R1600","R1800","R2000","R2500","R3000","R3500","R4000","R4500",
	"QCD700_1","QCD700_2","QCD1000_1","QCD1000_2","QCD1500_1","QCD1500_2","QCD2000_1","QCD2000_2"
	};
	
	int aa[40]={
		355,673,115,235,90,181,47,94
	};
	string dataPath="/data7/chchen/80x/JetHT/NCUGlobalTuples_";
	if(a==21){
		for(int j=0;j<11;j++){
			HH4bBtagMakeEffBase_80(1,2,st1[j],fileName[j],"");
			//HHbbbbBtagEffBase_76(1,2,st1[j],fileName[j],"");
		}
	}
	else if(a==30){
		HH4bBtagMakeEffBase_80(1,500,dataPath,"data1");	
		//HHbbbbBtagEffBase_76(1,500,dataPath,"data1");	
	}
	else if (a>21){
		HH4bBtagMakeEffBase_80(1,aa[a-22],st1[a],fileName[a],"");
	}
	/*
	
	string  fileName[40]={
	"B1000","B1200","B1400","B1600","B1800","B2000","B2500","B3000","B3500","B4000","B4500",
	"R1000","R1200","R1400","R1600","R1800","R2000","R2500","R3000","R3500","R4000","R4500",
	"QCD100_1","QCD100_2","QCD200","QCD300","QCD500","QCD700","QCD1000_1","QCD1000_2","QCD1500","QCD2000_1","QCD2000_2"
	};
	int aa0[40]={
		1,1,1,1,1,1,1,1,1,1,1,
		1,1,1,1,1,1,1,1,1,1,1,
		1,1000,1,1,1,
		1,1,1,1,1,1,
	};
	int aa[40]={
		2,2,2,2,2,2,2,2,2,2,2,
		2,2,2,2,2,2,2,2,2,2,2,
		1000,1925,459,498,473,
		357,126,244,100,56,103,
	};
	string data1="/data7/syu/NCUGlobalTuples/76X_data/1ff2146/JetHT-Run2015C-16Dec2015-v1/160422_131222/0000/NCUGlobalTuples_";
	string data2="/data7/syu/NCUGlobalTuples/76X_data/1ff2146/JetHT-Run2015D-16Dec2015-v1/160422_131233/0000/NCUGlobalTuples_";
	string data3="/data7/syu/NCUGlobalTuples/76X_data/1ff2146/JetHT-Run2015D-16Dec2015-v1/160422_131233/0001/NCUGlobalTuples_";
	if(a==33){	
	HHbbbbBtagEffBase_76(1,501,data2,"data2");	
	HHbbbbBtagEffBase_76(501,1000,data2,"data3");	
	}
	if(a==34){
		HHbbbbBtagEffBase_76(1,56,data1,"data1");	
	    HHbbbbBtagEffBase_76(1000,1500,data3,"data4");	
	}
	if(a==35)HHbbbbBtagEffBase_76(1501,1968,data3,"data5");	
	if (a==36){
		
		
	}
	if(a>=22 && a<=26)cout<<""<<endl;
	//else if(a<33)HHbbbbBtagMakeEff_76(aa0[a],aa[a],st1[a],fileName[a]);
	else if(a<33){
		HHbbbbBtagMakeEff_76(aa0[a],aa[a],st1[a],fileName[a],"");
		//HHbbbbBtagMakeEff_76(aa0[a],aa[a],st1[a],fileName[a],"JESUp");
		//HHbbbbBtagMakeEff_76(aa0[a],aa[a],st1[a],fileName[a],"JESDown");
		HHbbbbBtagEffBase_76(aa0[a],aa[a],st1[a],fileName[a],"");
		//HHbbbbBtagEffBase_76(aa0[a],aa[a],st1[a],fileName[a],"tau21Up");
		//HHbbbbBtagEffBase_76(aa0[a],aa[a],st1[a],fileName[a],"tau21Down");
	}
	*/
}