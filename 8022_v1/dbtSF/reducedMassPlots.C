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


TH1D* massPlotBase(string input,string variable,double down=100,double up=1000,int category=1){
	TH1D* th1=new TH1D("","",100,down,up);
	TFile* tf1;
	tf1=TFile::Open(Form("%s",input.data()));
	TTree* tree;
	tf1->GetObject("mynewTree",tree);
	TreeReader data(tree);
				//data.Print();
	for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
			data.GetEntry(jEntry);
			Float_t  vari = data.GetFloat(Form("%s",variable.data()));
			Float_t  jet1bbtag = data.GetFloat("jet1bbtag");
			Float_t  jet2bbtag = data.GetFloat("jet2bbtag");
			if((jet1bbtag<0.8 ||jet2bbtag<0.8 )&& category==1)continue;
			if((jet1bbtag<0.3 ||jet2bbtag<0.3 )&&(!(jet1bbtag>0.8 &&jet2bbtag>0.8 ))&& category==2)continue;
			//Float_t  jet2_puppi_msoftdrop = data.GetFloat("jet2_puppi_msoftdrop");
			th1->Fill(vari);
	}
	return th1;
}

TH1D* massPlotBaseSF(string input,string variable,double down=100,double up=1000,int category=1){
	TH1D* th1=new TH1D("","",100,down,up);
	TFile* tf1;
	tf1=TFile::Open(Form("%s",input.data()));
	TTree* tree;
	tf1->GetObject("mynewTree",tree);
	TreeReader data(tree);
				//data.Print();
	for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
			data.GetEntry(jEntry);
			Float_t  vari = data.GetFloat(Form("%s",variable.data()));
			Float_t  jet1bbtag = data.GetFloat("jet1bbtag");
			Float_t  jet2bbtag = data.GetFloat("jet2bbtag");
			if((jet1bbtag<0.8 ||jet2bbtag<0.8 )&& category==1)continue;
			if((jet1bbtag<0.3 ||jet2bbtag<0.3 )&&(!(jet1bbtag>0.8 &&jet2bbtag>0.8 ))&& category==2)continue;
			Float_t  dbtSF = data.GetFloat("dbtSF");
			//Float_t  jet2_puppi_msoftdrop = data.GetFloat("jet2_puppi_msoftdrop");
			th1->Fill(vari,dbtSF);
	}
	return th1;
}

void setTGraph(TH1D* tg1,int i,bool setMax=0){
	
	tg1->SetLineColor(i+1);
	//tg1->SetLineColor(98-4*i);
	tg1->SetMarkerColor(i+1);
	if(i==4)tg1->SetLineColor(kOrange);
	if(i==4)tg1->SetMarkerColor(kOrange);
	//tg1->SetMarkerColor(98-4*i);
	tg1->SetMarkerStyle(20);
	
	tg1->GetXaxis()->SetTitle("m_{jj}[GeV]");
	tg1->GetYaxis()->SetTitle("");
	//limits->GetZaxis()->SetTitle("#sigma_{95\% CL} / #sigma_{th}");
	tg1->SetTitle("");
	tg1->SetLineWidth(2);
	tg1->SetFillColor(0);
	//tg1->SetMaximum(1);
	//tg1->SetMinimum(0.0014);
	tg1->SetMarkerSize(1);
	//tg1->SetMarkerStyle(20);
	//tg1->SetMinimum(0.05);
	tg1->GetYaxis()->SetTitleOffset(1.5);
	//tg1->GetZaxis()->SetTitleOffset(0.65);
	// size of axis labels
	tg1->GetXaxis()->SetTitleSize(0.04);
	tg1->GetYaxis()->SetTitleSize(0.04);
	//limits->GetZaxis()->SetTitleSize(0.035);
	tg1->GetXaxis()->SetLabelSize(0.05);
	tg1->GetYaxis()->SetLabelSize(0.03); 
}


void setLeg(TLegend* leg){
	leg->SetBorderSize(1);
	//leg->SetLineColor(0);                                                                                                                     
	leg->SetLineStyle(1);
	leg->SetLineWidth(1);
	leg->SetFillColor(0);
	leg->SetFillStyle(1001);
}

void reducedMassPlots(){
	int massP[13]={750,800,900,1000,1200,1400,1600,1800,2000,2500,3000,4000,4500};
	TH1D* th1[2];
	TCanvas * c1 = new TCanvas("c1","",600,600);
	c1->cd();
	gStyle->SetOptStat(0);
	gStyle->SetPalette(57);
	gStyle->SetFrameLineWidth(3);
	gStyle->SetPadRightMargin(0.001);
	gStyle->SetPadLeftMargin(0.13);
	c1->SetLeftMargin(0.15);
	c1->SetRightMargin(0.01);
	gStyle->SetNdivisions(605, "XYZ");
	for(int i=0;i<13;i++){
		th1[0]=massPlotBase(Form("BulkGrav_M-%d_0.root",massP[i]),"dijetmass_softdrop_corr",massP[i]*0.6,massP[i]*1.2);
		th1[1]=massPlotBaseSF(Form("BulkGrav_M-%d_0.root",massP[i]),"dijetmass_softdrop_corr",massP[i]*0.6,massP[i]*1.2);
	
		th1[0]->Draw("");
		th1[1]->Draw("same");
		setTGraph(th1[0],0);
		setTGraph(th1[1],1);
		TLegend *leg;
		leg = new TLegend(.66, .76, .88, .88);
		setLeg(leg);
		leg->AddEntry(th1[0],"TT raw");
		leg->AddEntry(th1[1],"TT dbtSF");
		leg->Draw("same");
		//th1[1]->SetLineColor(2);
		c1->Print(Form("reducedMass/TT_grav_%d.pdf",massP[i]));
		
		th1[0]=massPlotBase(Form("BulkGrav_M-%d_0.root",massP[i]),"dijetmass_softdrop_corr",massP[i]*0.7,massP[i]*1.3,2);
		th1[1]=massPlotBaseSF(Form("BulkGrav_M-%d_0.root",massP[i]),"dijetmass_softdrop_corr",massP[i]*0.7,massP[i]*1.3,2);
	
		th1[0]->Draw("");
		th1[1]->Draw("same");
		setTGraph(th1[0],0);
		setTGraph(th1[1],1);
		leg->Clear();
		setLeg(leg);
		leg->AddEntry(th1[0],"LL raw");
		leg->AddEntry(th1[1],"LL dbtSF");
		leg->Draw("same");
		//th1[1]->SetLineColor(2);
		c1->Print(Form("reducedMass/LL_grav_%d.pdf",massP[i]));
	}
	
}