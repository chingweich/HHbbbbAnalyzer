#include <TFile.h>
#include <TLine.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TColor.h>
#include <TLegend.h>
#include <TPaletteAxis.h>
#include <TStyle.h>
#include <TMath.h>
#include <TFrame.h>
#include <TLatex.h>
//#include "../macro/mkPlotsLivia/CMS_lumi.C"
#include <iostream>
#include <vector>
#include "../../untuplizer.h"


double  makePlotBase(string input){
	TFile* tf1;
	tf1=TFile::Open(Form("%s",input.data()));
	TTree* tree;
	tf1->GetObject("mynewTree",tree);
	TreeReader data(tree);
	TH1F* th1=(TH1F*)tf1->Get("CountWeighted");
	double  temp=0;
	for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
			data.GetEntry(jEntry);
			Float_t  vari = data.GetFloat("dijetmass_softdrop_corr");
			if(vari>750)temp++;
			//Float_t  jet2_puppi_msoftdrop = data.GetFloat("jet2_puppi_msoftdrop");
			//th1->Fill(vari);
	}
	return temp/ th1->GetBinContent(1);
}

void makePlots(){
	double eff[10]={0};
	int masspoint[10]={1000,1200,1400,1600,1800,2000,2500,3000,4000,4500};
	double masspointdb[10]={1000,1200,1400,1600,1800,2000,2500,3000,4000,4500};
	for(int i=0;i<10;i++){
		eff[i]=makePlotBase(Form("Summer2016Garviton/BulkGrav_M-%d_0.root",masspoint[i]));
		cout<<eff[i]<<endl;
	}
	TGraph* tg1 =new TGraph(10,masspointdb,eff);
	TCanvas * cboth = new TCanvas("cboth","",550,550);
	cboth->cd();
	gStyle->SetOptStat(0);
	 //gStyle->SetPaintTextFormat("2.1f");
	gStyle->SetPalette(57);
	gStyle->SetFrameLineWidth(3);
	//gStyle->SetPadRightMargin(0.109);
	 //gStyle->SetPadLeftMargin(0.13);
	 //cboth->SetLeftMargin(0.1);
	 //cboth->SetRightMargin(0.1);
	gStyle->SetTitleOffset(2, "Z");
	gStyle->SetNdivisions(605, "XYZ");
	
	tg1->GetXaxis()->SetTitle("m_{X} [GeV]");
	tg1->GetYaxis()->SetTitle("Efficiency");
	//limits->GetZaxis()->SetTitle("#sigma_{95\% CL} / #sigma_{th}");
	tg1->SetTitle("");
	tg1->SetMaximum(0.7);
	tg1->SetMarkerSize(1);
	tg1->SetMarkerStyle(20);
	tg1->SetMinimum(0.4);
	tg1->GetYaxis()->SetTitleOffset(1);
	//tg1->GetZaxis()->SetTitleOffset(0.65);
	// size of axis labels
	tg1->GetXaxis()->SetTitleSize(0.04);
	tg1->GetYaxis()->SetTitleSize(0.04);
	//limits->GetZaxis()->SetTitleSize(0.035);
	tg1->GetXaxis()->SetLabelSize(0.05);
	tg1->GetYaxis()->SetLabelSize(0.04); 
	//limits->GetZaxis()->SetLabelSize(0.025);
	
	tg1->Draw("APL");
	cboth->Print("eff.pdf");
}