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
#include <fstream>
#include <vector>
#include <TFile.h>

void drawAndMakeTable(){
	double eff[10]={0},eff2[10]={0};
	int masspoint[10]={1000,1200,1400,1600,1800,2000,2500,3000,4000,4500};
	double masspointdb[10]={1000,1200,1400,1600,1800,2000,2500,3000,4000,4500};
	
	TGraph* tg1,* tg2;
	TFile* out=TFile::Open("eff.root");
	tg1=(TGraph*)out->Get("wo");
	tg2=(TGraph*)out->Get("with");
	
	for(int i=0;i<10;i++){
		double temp;
		tg1->GetPoint(i,temp,eff[i]);
		tg2->GetPoint(i,temp,eff2[i]);
		//eff2[i]=tg2->GetBinContent(i+1);
		cout<<eff[i]<<endl;
	}
	ofstream myfile;
	myfile.open ("effTable.txt");
	myfile<<"cut ";
	for(int i=0;i<10;i++)myfile<<" & "<<masspoint[i];
	myfile<<"\\\\"<<endl;
	myfile<<"\\hline"<<endl;
	
	myfile<<"without $\\tau$21 ";
	for(int i=0;i<10;i++)myfile<<" & "<<eff[i];
	myfile<<"\\\\"<<endl;
	myfile<<"with $\\tau$21 ";
	for(int i=0;i<10;i++)myfile<<" & "<<eff2[i];
	myfile<<"\\\\"<<endl;
	
	
	
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
	tg1->SetMaximum(0.15);
	tg1->SetMarkerSize(1);
	tg1->SetMarkerStyle(20);
	tg1->SetMinimum(0.05);
	tg1->GetYaxis()->SetTitleOffset(1);
	//tg1->GetZaxis()->SetTitleOffset(0.65);
	// size of axis labels
	tg1->GetXaxis()->SetTitleSize(0.04);
	tg1->GetYaxis()->SetTitleSize(0.04);
	//limits->GetZaxis()->SetTitleSize(0.035);
	tg1->GetXaxis()->SetLabelSize(0.05);
	tg1->GetYaxis()->SetLabelSize(0.03); 
	//limits->GetZaxis()->SetLabelSize(0.025);
	
	tg1->Draw("APL");
	
	tg2->SetMarkerSize(1);
	tg2->SetMarkerStyle(20);
	tg2->SetMarkerColor(2);
	tg2->SetLineColor(2);
	tg2->SetFillColor(0);
	tg1->SetFillColor(0);
	tg2->Draw("PL same");
	
	TLegend *leg;
  leg = new TLegend(.60, .73, .88, .88);
  leg->SetBorderSize(1);
  //leg->SetLineColor(0);                                                                                                                     
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);


  leg->SetFillColor(0);
  leg->SetLineColor(0);
 
  leg->SetShadowColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->AddEntry(tg1,"without #tau21");
  leg->AddEntry(tg2,"with #tau21");
	leg->Draw("same");
	
	cboth->Print("eff.pdf");
	
}