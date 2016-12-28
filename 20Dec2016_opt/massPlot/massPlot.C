#include <TLegend.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TH1D.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TObject.h>
#include "TImage.h"
#include "TSystem.h"
#include "TStyle.h"
#include "../../untuplizer.h"
#include <TClonesArray.h>
#include <fstream>
#include <cmath>
#include <TSystem.h>
#include <string>
#include <sstream>
#include "../../setNCUStyle.C"


TH1D* massPlotBase(string input,string variable){
	TH1D* th1=new TH1D("","",100,50,150);
	TFile* tf1;
	tf1=TFile::Open(Form("%s",input.data()));
	TTree* tree;
	tf1->GetObject("mynewTree",tree);
	TreeReader data(tree);
				//data.Print();
	for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
			data.GetEntry(jEntry);
			Float_t  vari = data.GetFloat(Form("%s",variable.data()));
			//Float_t  jet2_puppi_msoftdrop = data.GetFloat("jet2_puppi_msoftdrop");
			th1->Fill(vari);
	}
	return th1;
}

void massPlot(){
	TH1D* th1[3];
	th1[0]=massPlotBase("BulkGrav_M-2000_0.root","jet1_puppi_msoftdrop");
	th1[1]=massPlotBase("BulkGrav_M-2000_0.root","jet1_puppi_msoftdrop_TheaCorr");
	th1[2]=massPlotBase("HCorr.root","jet1_puppi_msoftdrop_TheaCorr");
	
	TCanvas* c1;
	setNCUStyle(true);
	c1 = new TCanvas("c1","",800,600);
	
		TLegend *leg = new TLegend(0.20, 0.33, 0.35, 0.88);
  
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  string st[3]={"raw","Thea","Higgs"};
	for(int i=0;i<3;i++){
		leg->AddEntry(th1[i],Form("%s",st[i].data()));
		th1[i]->SetLineColor(i+1);
		th1[i]->SetLineWidth(3);
		th1[i]->SetMarkerSize(0);
		TF1 *tf1;
		tf1=new TF1("fa1","gaus(25000)",90,150);
		tf1->SetLineColor(i+1);
		th1[i]->Fit(tf1,"","",90,150);
		
		gStyle->SetOptStat(0);
		gStyle->SetOptFit(0);
		leg->AddEntry((TObject*) 0,Form("mean=%f",tf1->GetParameter(1)),"");
		leg->AddEntry((TObject*) 0,Form("#sigma=%f",tf1->GetParameter(2)),"");
		leg->AddEntry((TObject*) 0,Form("#sigma/mean=%f",tf1->GetParameter(2)/tf1->GetParameter(1)),"");
		if(i==0)th1[i]->Draw("");
		else th1[i]->Draw("same");
	}
	for(int i=0;i<3;i++){
		if(i==0)th1[i]->Draw("");
		else th1[i]->Draw("same");
	}
	leg->Draw("same");
	c1->Print("j0.pdf");
	
	th1[0]=massPlotBase("BulkGrav_M-2000_0.root","jet2_puppi_msoftdrop");
	th1[1]=massPlotBase("BulkGrav_M-2000_0.root","jet2_puppi_msoftdrop_TheaCorr");
	th1[2]=massPlotBase("HCorr.root","jet2_puppi_msoftdrop_TheaCorr");
	
leg->Clear();
	for(int i=0;i<3;i++){
		leg->AddEntry(th1[i],Form("%s",st[i].data()));
		th1[i]->SetLineColor(i+1);
		th1[i]->SetLineWidth(3);
		th1[i]->SetMarkerSize(0);
		TF1 *tf1;
		tf1=new TF1("fa1","gaus(25000)",90,150);
		tf1->SetLineColor(i+1);
		th1[i]->Fit(tf1,"","",90,150);
		
		gStyle->SetOptStat(0);
		gStyle->SetOptFit(0);
		leg->AddEntry((TObject*) 0,Form("mean=%f",tf1->GetParameter(1)),"");
		leg->AddEntry((TObject*) 0,Form("#sigma=%f",tf1->GetParameter(2)),"");
		leg->AddEntry((TObject*) 0,Form("#sigma/mean=%f",tf1->GetParameter(2)/tf1->GetParameter(1)),"");
		if(i==0)th1[i]->Draw("");
		else th1[i]->Draw("same");
	}
	for(int i=0;i<3;i++){
		if(i==0)th1[i]->Draw("");
		else th1[i]->Draw("same");
	}
	leg->Draw("same");
	c1->Print("j1.pdf");
	
}