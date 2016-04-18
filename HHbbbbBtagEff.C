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
//#include "BTagCalibrationStandalone.h"
#include "HHbbbbBtagMakeEff.C"
//#include "dataPrinter.C"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"
#include "HHbbbbBtagEffBase.C"
#include "HHbbbbBtagEffBaseO.C"
//#include "HHbbbbBtagAdditionThinJet.C"
#include "HHbbbbAnalyzerBaseC.C"
TCanvas* c1;

void HHbbbbBtagEff(int a){
	//setNCUStyle();
	
	//setNCUStyle();
	
	string st1[40]={
		/*0-11*/
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-1000_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-1200_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-1400_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-1600_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-1800_13TeV-madgraph.root",
	    "/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-2000_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-2500_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-3000_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-3500_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-4000_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-4500_13TeV-madgraph.root",
		/*12-21*/
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-1000_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-1200_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-1400_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-1600_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-1800_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-2000_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-2500_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-3000_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-3500_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-4000_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-4500_13TeV-madgraph.root",
		/*22-30*/
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/QCD_HTBinned/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/160306_060112/0000/NCUGlobalTuples_",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/QCD_HTBinned/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/160306_060112/0001/NCUGlobalTuples_",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/QCD_HTBinned/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/160306_060145/0000/NCUGlobalTuples_",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/QCD_HTBinned/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/160306_060215/0000/NCUGlobalTuples_",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/QCD_HTBinned/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/160306_060251/0000/NCUGlobalTuples_",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/QCD_HTBinned/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/160308_130001/0000/NCUGlobalTuples_",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/QCD_HTBinned/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/160308_130035/0000/NCUGlobalTuples_",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/QCD_HTBinned/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/160308_130118/0000/NCUGlobalTuples_",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/QCD_HTBinned/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/160308_130153/0000/NCUGlobalTuples_",
		
	};
	string  masspoint[11]={"1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500"};
	string  fileName[40]={"B1000","B1200","B1400","B1600","B1800","B2000","B2500","B3000","B3500","B4000","B4500",
	"R1000","R1200","R1400","R1600","R1800","R2000","R2500","R3000","R3500","R4000","R4500",
	"QCD100_1","QCD100_2","QCD200","QCD300","QCD500","QCD700","QCD1000","QCD1500","QCD2000"};
	double eff1,eff2,eff3,eff4,eff5;
	TH1D* th1[30];
	int aa0[40]={
		1,1,1,1,1,
		1,1,1,1,1,
		1,
		1,1,1,1,1,
		1,1,1,1,1,
		1,
		1,1000,1,1,1,
		1,1,1,1,
	};
	
	int aa[40]={
		2,2,2,2,2,
		2,2,2,2,2,
		2,
		2,2,2,2,2,
		2,2,2,2,2,
		2,
		1000,1904,432,467,449,
		346,120,90,45,
		
	};
	double Xsec[40]={
		0.1,0.1,0.1,0.1,0.1,
		0.1,0.1,0.1,0.1,0.1,
		0.1,
		0.1,0.1,0.1,0.1,0.1,
		0.1,0.1,0.1,0.1,0.1,
		0.1,
		27990000,27990000,1712000,347700,32100,
		6831,1207,119.9,25.24
	};
	bool sigName=0;
	if (a<22)sigName=1;
	
	if(a==100){
	string data1="/data7/syu/NCUGlobalTuples/Run2015D/eec7461/JetHT/crab_JetHT-Run2015D-05Oct2015-v1/160223_142842/0000/NCUGlobalTuples_";
	string data2="/data7/syu/NCUGlobalTuples/Run2015D/eec7461/JetHT/crab_JetHT-Run2015D-PromptReco-v4/160224_140926/0000/NCUGlobalTuples_";
	string data2_2="/data7/syu/NCUGlobalTuples/Run2015D/eec7461/JetHT/crab_JetHT-Run2015D-PromptReco-v4/160224_140926/0001/NCUGlobalTuples_";
	string dataName1="JetHT-Run2015D-05Oct2015-v1";
	string dataName2="JetHT-Run2015D-PromptReco-v4";
	HHbbbbBtagMakeEff(1,412,data1,dataName1,1,2);	
	//HHbbbbBtagMakeEff(1,501,data2,Form("%s1",dataName2.data()),1,2);	
	//HHbbbbBtagMakeEff(501,1000,data2,Form("%s2",dataName2.data()),1,2);	
	HHbbbbBtagMakeEff(1000,1094,data2_2,Form("%s3",dataName2.data()),1,2);	
	//dataPrinter(412,data1,dataName1,1,0);	
	//dataPrinter(856,data2,dataName2,1,0);	
	}
	
	bool sig=0;
	if (a<22)sig=1;
	
	//th1[a]=HHbbbbAnalyzerCompare(aa[a],st1[a],fileName[a],Xsec[a],eff1,eff2,eff3,eff4,eff5,sig);
	//th1[a]=HHbbbbAnalyzerCSVEff(aa[a],st1[a],fileName[a],Xsec[a],eff1,eff2,eff3,eff4,eff5,sig);
	//th1[a]=HHbbbbAnalyzerJetEff(aa[a],st1[a],fileName[a],Xsec[a],eff1,eff2,eff3,eff4,eff5,sig);
	//th1[a]=HHbbbbAnalyzerBaseDCut(aa[a],st1[a],fileName[a],Xsec[a],eff1,eff2,eff3,eff4,eff5,sig);
	//th1[a]=HHbbbbAnalyzerBaseHPHP(aa[a],st1[a],fileName[a],Xsec[a],eff1,eff2,eff3,eff4,eff5,sig);
	if(!(a==22||a==23||a==24||a==25||a==26||a==100))HHbbbbAnalyzerBaseC(aa0[a],aa[a],st1[a],fileName[a],Xsec[a],sigName);
	//dataPrinter(aa[a],st1[a],fileName[a],Xsec[a],sig);
	
}