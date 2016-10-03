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

TCanvas* c1;

void BDTLimit_noCom(){
	setNCUStyle(true);
	c1 = new TCanvas("c1","",1360,768);
	TFile *f[40];
	f[0]=TFile::Open("Graviton_subtr_3btag_cat1TMe.root");
	f[1]=TFile::Open("Graviton_subtr_3btag_cat1TLe.root");
	f[2]=TFile::Open("Graviton_subtr_3btag_cat1MMe.root");
	f[3]=TFile::Open("Graviton_subtr_3btag_cat1MLe.root");
	f[4]=TFile::Open("Graviton_subtr_3btag_cat1LLe.root");
	
	f[5]=TFile::Open("../BDTLimit/Graviton_subtr_4btag_cat0TT.root");
	f[6]=TFile::Open("../BDTLimitTT/Graviton_subtrcsv.root");
	
	TGraphAsymmErrors* tg1[8];
	
	TLegend *leg = new TLegend(0.65, 0.58, 0.96, 0.92);
  
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);

string st[7]={"TM","TL","MM","ML","LL","TT","3b+4b"};
  
	for(int i=0;i<7;i++){
		tg1[i]=(TGraphAsymmErrors*)f[i]->Get("LimitExpectedCLs");
		tg1[i]->SetLineColor(i+1);
		tg1[i]->SetLineWidth(2);
		tg1[i]->SetLineStyle(1);
		tg1[i]->SetFillColor(0);
		tg1[i]->SetMinimum(2);
		tg1[i]->SetMaximum(1000);
		tg1[i]->GetYaxis()->SetTitle("95% CLs on #sigma(X#rightarrowHH)#timesBR(HH#rightarrowb#bar{b}b#bar{b})[fb]");
		if(i==0)tg1[i]->Draw("APL");
		else tg1[i]->Draw("PL same");
		leg->AddEntry(tg1[i],Form("%s",st[i].data()));
	}
	leg->Draw("same");
	
	c1->SetLogy();
	c1->Print("BDTLimitwoCombine.pdf");
}
