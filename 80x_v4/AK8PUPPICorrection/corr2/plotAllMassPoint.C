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
#include <string>
#include <sstream>
#include"../../setNCUStyle.C"
#include "TF1.h"
TCanvas* c1;

void plotAllMassPoint(){
	setNCUStyle(true);
	c1 = new TCanvas("c1","",1360,768);
	string  masspoint[11]={"1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500"};
	int  masspointy[11]={1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500};
	TFile *f[40];
	for(int i=0;i<11;i++){
		f[i]=TFile::Open(Form("B%s.root",masspoint[i].data()));
		//f[i+11]=TFile::Open(Form("R%s.root",masspoint[i].data()));
	}
	
	TLegend *leg = new TLegend(0.28, 0.65, 0.48, 0.90);
  
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
	for(int i=0;i<11;i++){
		if(i==0||i==8)continue;
		TH1D* th1=(TH1D*)f[i]->Get("mass");
		th1->SetLineColor(97-4*i);
		if(i==1)th1->SetMaximum(th1->GetMaximum()*1.8);
		th1->Scale(1/th1->Integral());
		th1->Draw("same");
		th1->SetLineWidth(2);
		leg->AddEntry(th1,Form("M=%s",masspoint[i].data()));
	}
	leg->Draw("same");
	c1->Print("SDMass.pdf");
	
}
