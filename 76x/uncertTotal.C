#include <TLegend.h>
#include <vector>
#include <iostream>
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
#include "TImage.h"
#include "TSystem.h"
#include "TStyle.h"
//#include "untuplizer.h"
#include <TClonesArray.h>
#include <fstream>
#include <cmath>
#include <TSystem.h>
#include <string>
#include <sstream>
#include "../setNCUStyle.C"
#include<TH2.h>
#include "TLine.h"
#include "TF1.h"

TCanvas* c1;
void uncertBase(string unType,string outputDir,double max,double min,int option=0){
	setNCUStyle();
	c1 = new TCanvas("c1","",1360,768);
	
	TFile* tf1[11],* tf2[11];
	TFile* tf3[11],* tf4[11];
	TFile* tf5[11],* tf6[11];
	string st[11]={"1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500"};
	
	for (int i=0;i<11;i++)tf1[i]=TFile::Open(Form("sf/B%s.root",st[i].data()));
	for (int i=0;i<11;i++)tf3[i]=TFile::Open(Form("sf/B%s_%sUp.root",st[i].data(),unType.data()));
	for (int i=0;i<11;i++)tf5[i]=TFile::Open(Form("sf/B%s_%sDown.root",st[i].data(),unType.data()));
	for (int i=0;i<11;i++)tf2[i]=TFile::Open(Form("sf/R%s.root",st[i].data()));
	double x[11]={1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500};
	double xx[11]={1000,1200,1400,1600,1800,2000,2500,3000,4000,4500};
	double y[6][11];

	for (int i=0;i<11;i++){
		TH1D * th1= (TH1D *)tf1[i]->FindObjectAny("totalMass");
		y[0][i]=th1->GetMean();
		TH1D * th2;
		if(option==0)th2=(TH1D *)tf3[i]->FindObjectAny("totalMass");
		else th2=(TH1D *)tf1[i]->FindObjectAny("totalMass_pileup_up");
		y[1][i]=th2->GetMean();
		TH1D *th3;
		if(option==0)th3=(TH1D *)tf5[i]->FindObjectAny("totalMass");
		else th3=(TH1D *)tf1[i]->FindObjectAny("totalMass_pileup_down");
		y[2][i]=th3->GetMean();
		y[1][i]/=y[0][i];
		y[2][i]/=y[0][i];
		y[0][i]/=y[0][i];
		
		th2->SetLineColor(2);
		th3->SetLineColor(3);
		TLegend* leg ;
		leg=new TLegend(0.751452,0.712447,0.880645,0.863966);
		leg->AddEntry(th1,Form("%s_central",unType.data()));
		leg->AddEntry(th2,Form("%s_up",unType.data()));
		leg->AddEntry(th3,Form("%s_down",unType.data()));
		
		
		int bmin=0,bmax=0;

		for (int k=1;k<25001;k++){
			bmin=k;
			if (th1->GetBinContent(k)/th1->GetMaximum()>0.02) break;
	}

		for (int k=25000;k>0;k--){
			bmax=k;
			if (th1->GetBinContent(k)/th1->GetMaximum()>0.02) break;
	}
	
	cout<<980+bmin*20<<","<<980+bmax*20<<endl;
		th2->GetXaxis()->SetRangeUser(980+bmin*20,980+bmax*20);
		
		th2->SetTitle(Form("M%4.0f",x[i]));
		
		
		th2->Draw("HIST,c");
		th1->Draw("HIST,same,c");
		th3->Draw("HIST,same,c");
		leg->Draw("same");
		if(i==0)c1->Print(Form("%s/plot/mjj.pdf(",outputDir.data()));
			else if (i==10)c1->Print(Form("%s/plot/mjj.pdf)",outputDir.data()));
				else c1->Print(Form("%s/plot/mjj.pdf",outputDir.data()));
	}
	
	double yy[6][10];
	for (int i=0;i<11;i++){
		if(i>=8){
			yy[0][i]=y[0][i+1];
			yy[1][i]=y[1][i+1];
			yy[2][i]=y[2][i+1];
		}
		else {
			yy[0][i]=y[0][i];
			yy[1][i]=y[1][i];
			yy[2][i]=y[2][i];
		}
	}
	
	//TGraph * tg1=new TGraph(11,x,y[0]);
	//TGraph * tg2=new TGraph(11,x,y[1]);
	//TGraph * tg3=new TGraph(11,x,y[2]);
	TGraph * tg1=new TGraph(10,xx,yy[0]);
	TGraph * tg2=new TGraph(10,xx,yy[1]);
	TGraph * tg3=new TGraph(10,xx,yy[2]);
	
	tg1->SetTitle("JES uncertainty");
	tg2->SetTitle("JES_up fit");
	tg3->SetTitle("JES_down fit");
	tg1->SetLineColor(1);
	tg2->SetLineColor(2);
	tg3->SetLineColor(3);
	
	tg1->SetFillColor(kWhite);
	tg2->SetFillColor(kWhite);
	tg3->SetFillColor(kWhite);
	
	tg1->GetXaxis()->SetTitle("M_{X}[GeV]");
	tg1->GetYaxis()->SetTitle("Mean/Mean_{central}");	
	tg1->SetMaximum(max);
	tg1->SetMinimum(min);
	tg1->Draw("APL");
	tg2->Draw("PL,same");
	tg3->Draw("PL,same");
	
	TLegend* leg ;
	leg=new TLegend(0.391452,0.712447,0.590645,0.863966);
	leg->AddEntry(tg1,Form("%s_central",unType.data()));
	leg->AddEntry(tg2,Form("%s_up",unType.data()));
	leg->AddEntry(tg3,Form("%s_down",unType.data()));
	
	leg->Draw("same");
	
	c1->Print(Form("%s/plot/all.pdf",outputDir.data()));
	
	TF1 *fa1 = new TF1("fa1","[0]+[1]*x",1000,4500); 
	TF1 *fa2 = new TF1("fa1","[0]+[1]*x",1000,4500); 
	tg2->Fit(fa1);
	fa1->SetLineColor(4);
	
	tg2->Draw("APL");
	fa1->Draw("same");
	c1->Print(Form("%s/plot/fit_up.pdf",outputDir.data()));
	//c1->Print("plot/fit_up.pdf");
	
	tg3->Fit(fa2);
	fa2->SetLineColor(4);
	tg3->Draw("APL");
	fa2->Draw("same");
	c1->Print(Form("%s/plot/fit_down.pdf",outputDir.data()));
	//c1->Print("plot/fit_down.pdf");

}

void uncertTotal(){
	uncertBase("JES","1_pileup",1.001,0.999,1);
	uncertBase("JES","2_JES",1.02,0.98);
	uncertBase("Btag","3_btag",1.001,0.999);
	uncertBase("tau21","4_tau21",1.002,0.998);
}