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
#include"../../../setNCUStyle.C"

TCanvas* c1;

#define  nWidth 4
#define  nBmin 4

void Limit(){
	setNCUStyle(true);
	c1 = new TCanvas("c1","",800,600);
	
	int width [nWidth]={25,30,35,40};
	int bmin[nBmin]={100,105,110,115};
	string  masspoint[11]={"1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500"};
	
	TLegend *leg = new TLegend(0.75, 0.68, 0.96, 0.95);
  
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  
  string categoryName[7]={"0a","0c","1a","1c","2a","2b","2d"};
	
	for(int h=0;h<7;h++){
		TFile* tf1=TFile::Open(Form("%s/MassPlotFineBins_subtr_Moriond_Silver105to135.root",categoryName[h].data()));
		if(h==2) tf1=TFile::Open(Form("%s/MassPlotFineBins_subtr_Moriond_Silver110to145.root",categoryName[h].data()));
		if(h>=3) tf1=TFile::Open(Form("%s/MassPlotFineBins_subtr_Moriond_Silver115to145.root",categoryName[h].data()));
			 if (!tf1 || !tf1->IsOpen())continue;
			 TGraphAsymmErrors* tg1=(TGraphAsymmErrors*)tf1->Get("LimitExpectedCLs");
		tg1->GetYaxis()->SetTitle("95% CLs on #sigma(X#rightarrowHH)#timesBR(HH#rightarrowb#bar{b}b#bar{b})[fb]");
			 tg1->SetLineStyle(1);
			 tg1->SetFillColor(0);
			 tg1->SetLineColor(h+1);
			 if (h<2)leg->AddEntry(tg1,Form("%s 105-135",categoryName[h].data()));
			 else if(h==2)leg->AddEntry(tg1,Form("%s 110-145",categoryName[h].data()));
			 else if(h>=3)leg->AddEntry(tg1,Form("%s 115-145",categoryName[h].data()));
			 tg1->SetMaximum(12);
			 tg1->SetMinimum(0.5);
			// cout<<k<<","<<m<<endl;
			 //c1->SetLogy(1);
			 if(h==0)tg1->Draw("APL");
			 else tg1->Draw("PL same");
	}
  
	
	leg->Draw("same");
	c1->Print("Limit.pdf");
}
