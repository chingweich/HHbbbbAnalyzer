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

void BDTLimit(){
	setNCUStyle(true);
	c1 = new TCanvas("c1","",1360,768);
	TFile *f[40];
	f[0]=TFile::Open("Graviton_subtr_4btag_cat0LL.root");
	f[1]=TFile::Open("Graviton_subtr_4btag_cat0ML.root");
	f[2]=TFile::Open("Graviton_subtr_4btag_cat0MM.root");
	f[3]=TFile::Open("Graviton_subtr_4btag_cat0TL.root");
	f[4]=TFile::Open("Graviton_subtr_4btag_cat0TM.root");
	f[5]=TFile::Open("Graviton_subtr_4btag_cat0TT.root");
	
	TGraphAsymmErrors* tg1[6];
	
	TLegend *leg = new TLegend(0.75, 0.68, 0.92, 0.87);
  
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);

string st[6]={"LL","ML","MM","TL","TM","TT"};
  
	for(int i=0;i<6;i++){
		tg1[i]=(TGraphAsymmErrors*)f[i]->Get("LimitExpectedCLs");
		tg1[i]->SetLineColor(i+1);
		tg1[i]->SetFillColor(0);
		tg1[i]->SetMinimum(1);
		if(i==0)tg1[i]->Draw("APL");
		else tg1[i]->Draw("PL same");
		leg->AddEntry(tg1[i],Form("%s",st[i].data()));
	}
	leg->Draw("same");
	
	c1->SetLogy();
	c1->Print("BDLimit.pdf");
}
