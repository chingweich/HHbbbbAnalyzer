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
#include "../setNCUStyle.C"

TCanvas* c1;

void drawEffDSV(){
	setNCUStyle(1);
	c1 = new TCanvas("c1","",1360,768);
	TFile *f[40];
	TH1D * th1[40];
	double db[9][9],db2[9][9];
	double Xsec[4]={6831,1207,119.9,25.24};
	
	f[9]=TFile::Open("sf/QCD700.root");
	f[10]=TFile::Open("sf/QCD1000.root");
	f[11]=TFile::Open("sf/QCD1500.root");
	f[12]=TFile::Open("sf/QCD2000.root");
	
	th1[9]=(TH1D*)f[9]->FindObjectAny("cutflow");
	th1[10]=(TH1D*)f[10]->FindObjectAny("cutflow");
	th1[11]=(TH1D*)f[11]->FindObjectAny("cutflow");
	th1[12]=(TH1D*)f[12]->FindObjectAny("cutflow");
	double fixNumber=20991/30584.3;
	for (int i=0;i<4;i++)th1[i+9]->Scale(fixNumber* 12883.846147301*Xsec[i]/th1[i+9]->GetBinContent(1));
	
	th1[9]->Add(th1[10]);
	th1[9]->Add(th1[11]);
	th1[9]->Add(th1[12]);
	
	string  masspoint[11]={"1200","1400","1600","1800","2000","2500","3000","4000","4500"};
	double   masspointx[9]={1200,1400,1600,1800,2000,2500,3000,4000,4500};
	
	for(int i=0;i<9;i++){
		f[i]=TFile::Open(Form("sf/Bulk%s.root",masspoint[i].data()));
		th1[i]=(TH1D*)f[i]->FindObjectAny("cutflow");
		for(int j=0;j<9;j++){
			db[j][i]=th1[i]->GetBinContent(j+1);
			if (j!=0)db[j][i]/=db[0][i];
			db2[j][i]=th1[i]->GetBinContent(j+1)/sqrt(th1[9]->GetBinContent(j+1));
		}
		db[7][i]=db[2][i]+db[3][i];
		db[8][i]=db[4][i]+db[6][i];
		
		db2[7][i]=(th1[i]->GetBinContent(3)+th1[i]->GetBinContent(4))/sqrt(th1[9]->GetBinContent(3)+th1[9]->GetBinContent(4));
		db2[8][i]=(th1[i]->GetBinContent(5)+th1[i]->GetBinContent(7))/sqrt(th1[9]->GetBinContent(5)+th1[9]->GetBinContent(7));
	}
	
	
	
	
	int color[8]={kBlack,kBlue,kBlue+3,kRed,kGreen,kOrange,kGreen,kGreen+3};
	int marker[8]={20,20,21,20,20,21,22,23};
	
	TGraph * tg1[8];
	for (int i=0;i<8;i++){
		tg1[i]=new TGraph(9,masspointx,db[i+1]); 
		tg1[i]->SetLineColor(color[i]);
		tg1[i]->GetXaxis()->SetTitle("M_{jj}[GeV]");
		tg1[i]->GetYaxis()->SetTitle("Efficiency");	
		tg1[i]->SetMarkerColor(color[i]);
		tg1[i]->SetMarkerSize(tg1[i]->GetMarkerSize()*1.8);
		tg1[i]->SetMarkerStyle(marker[i]);
		tg1[i]->SetFillColor(0);
		if(i==0){
			tg1[i]->SetMinimum(0.01);
			tg1[i]->Draw("APL");
		}
		else if (i==4){}
		else tg1[i]->Draw("samePL");
	}
	
	 TLegend *leg = new TLegend(0.75, 0.60, 0.92, 0.87);
  
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);
	leg->SetTextSize(0.04);
	leg->AddEntry(tg1[0], "preselection");
	leg->AddEntry(tg1[1], "3b");
	leg->AddEntry(tg1[2], "4b");
	leg->AddEntry(tg1[3], "2DSV");
	//leg->AddEntry(tg1[4], "exact 1DSV");
	leg->AddEntry(tg1[5], "1DSV+1b||1DSV+2b");
	leg->AddEntry(tg1[6], "3b+4b");
	leg->AddEntry(tg1[7], "2DSV+(1DSV+1b||1DSV+2b)");
	leg->Draw("same");
	
	c1->Print("plot/DSV.pdf");
	
	TGraph * tg2[8];
	for (int i=0;i<8;i++){
		tg2[i]=new TGraph(9,masspointx,db2[i+1]); 
		tg2[i]->SetLineColor(color[i]);
		tg2[i]->GetXaxis()->SetTitle("M_{jj}[GeV]");
		tg2[i]->GetYaxis()->SetTitle("S/#sqrt{B}");	
		tg2[i]->SetMarkerColor(color[i]);
		tg2[i]->SetMarkerStyle(marker[i]);
		tg2[i]->SetFillColor(0);
		if(i==0){
			tg2[i]->SetMaximum(600);
			tg2[i]->SetMinimum(0);
			tg2[i]->Draw("APL");
		}
		else tg2[i]->Draw("samePL");
	}
	leg->Draw("same");
	c1->Print("plot/DSV2.pdf");
}
