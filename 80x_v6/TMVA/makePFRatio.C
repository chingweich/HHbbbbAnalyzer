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
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TROOT.h>
#include "TImage.h"
#include "TSystem.h"
#include "TStyle.h"
#include <TClonesArray.h>
#include <fstream>
#include <cmath>
#include <TSystem.h>
#include <TF1.h>
#include <string>
#include <sstream>
#include "../../setNCUStyle.C"

#define  nWidth 5
#define  nBmin 11

#define nMass 3
#define nCat 4

TCanvas* c1;

void makePFRatio(){
	setNCUStyle(true);
	c1 = new TCanvas("c1","",800,600);
	TFile* tf1=TFile::Open("PFRatio/data.root");
	
	string massName[nMass]={"Thea","HCorr","Reg"};
	string catName[nCat]={"PP","PF","FP","FF"};
	string tau21Name[2]={"withTau21","woTau21"};
	string catNameShort[nCat]={"P","F"};
	string looseTight[2]={"loose","tight"};
	
	for(int i=0;i<nMass;i++){
		for(int j=0;j<2;j++){
			for(int k=0;k<2;k++){
				for(int w=0;w<2;w++){
					TH1D* th1,* th2;
				th1= (TH1D*) tf1->Get(Form("%s_%s_P%s_%s",looseTight[w].data(),massName[i].data(),catNameShort[j].data(),tau21Name[k].data()));
				th2= (TH1D*) tf1->Get(Form("%s_%s_F%s_%s",looseTight[w].data(),massName[i].data(),catNameShort[j].data(),tau21Name[k].data()));
				
				TF1 *fa = new TF1("fa","[0]+[1]*x+[2]*x*x+[3]*pow(x,3)",-3,3);
				
				th1->Divide(th2);
				th1->Fit(fa);
				
				ofstream  file ; 
				file.open(Form("PFRatio/%s_%s_%s_%s.txt",looseTight[w].data(),massName[i].data(),catNameShort[j].data(),tau21Name[k].data()));
				file<<fa->GetParameter(0)<<endl;
				file<<fa->GetParameter(1)<<endl;
				file<<fa->GetParameter(2)<<endl;
				file<<fa->GetParameter(3)<<endl;
				file.close();
				
				th1->Draw();
				c1->Print(Form("PFRatio/%s_%s_%s_%s.pdf",looseTight[w].data(),massName[i].data(),catNameShort[j].data(),tau21Name[k].data()));
				}
			}
		}
	}
	
	
}
