#include <vector>
#include <iostream>
#include <algorithm>
#include <TH1D.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1.h>
#include <TROOT.h>
#include "TImage.h"
#include "TSystem.h"

void pileupPrinter(string st){
	
	//TFile *f2=TFile::Open("MyDataPileupHistogram_true_down.root");
	TFile *f2=TFile::Open(Form("%s",st.data()));
	
	TH1D* sf2d[10];
	sf2d[0]=(TH1D*)f2->FindObjectAny("pileup");
	for(int i=1;i<60;i++){
		cout<<sf2d[0]->GetBinContent(i)<<","<<endl;
	}
	
}