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

#define nMassp 10
#define nSt 11

TCanvas* c1;

/*
void HHbbbbDrawerNum(int a,int b,int c,int d){
	
	TFile* tf1[11];
	string st[11]={"1000","1200","1400","1600","2000","2500","3000","3500","4000","4500","QCD"};
	int st2[10]={1000,1200,1400,1600,2000,2500,3000,3500,4000,4500};
	for(int i=0;i<11;i++)tf1[i]=TFile::Open(Form("root_files_Op/%s.root",st[i].data()));
	
	
	
	TH2F *th2f1,*th2f2;
	for(int i=0;i<10;i++){
		
		
		TH2F * th2f1 = (TH2F*)tf1[i]->FindObjectAny(Form("%s_jetPRmassL2L3Corr",st[i].data()));
		TH2F * th2f2 = (TH2F*)tf1[10]->FindObjectAny(Form("%s_jetPRmassL2L3Corr",st[i].data()));
		cout<<st[i].data()<<"="<<th2f1->Integral(a,b,c,d)<<endl;
		cout<<st[i].data()<<"="<<th2f2->Integral(a,b,c,d)<<endl;
		
		
		ofstream myfile;
		myfile.open (Form("opt_%d_%d_%d_%d/cout_%d.txt",a,b,c,d,st2[i]));
		myfile<<"DATA 0"<<endl;
		myfile<<"QCD "<<th2f2->Integral(a,b,c,d)<<endl;
		for(int j=0;j<10;j++)myfile<<"M"<<st2[j]<<" "<<th2f1->Integral(a,b,c,d)<<endl;
		myfile.close();
		
	}
}
*/


void HHbbbbDrawerNum(){
	setNCUStyle();
	c1 = new TCanvas("c1","",1360,768);
	
	TFile* tf1[12],* tf2[12],* tf3[12],* tf4[12];
	string st[12]={"1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500","QCD"};
	string st3[12]={"SIG_M1000","SIG_M1200","SIG_M1400","SIG_M1600","SIG_M1800","SIG_M2000","SIG_M2500","SIG_M3000","SIG_M3500","SIG_M4000","SIG_M4500","QCD"};
	int st2[11]={1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500};
	for(int i=0;i<12;i++){
		if(i<11) tf1[i]=TFile::Open(Form("root_files_cat_newNtuple/%s.root",st[i].data()));
		else tf1[i]=TFile::Open(Form("root_files_cat_newNtuple/%s.root",st[i].data()));
		if(i<11) tf2[i]=TFile::Open(Form("root_files_cat_UpNtuple/%s.root",st[i].data()));
		else tf2[i]=TFile::Open(Form("root_files_cat_newNtuple/%s.root",st[i].data()));
		if(i<11) tf3[i]=TFile::Open(Form("root_files_cat_DownNtuple/%s.root",st[i].data()));
		else tf3[i]=TFile::Open(Form("root_files_cat_newNtuple/%s.root",st[i].data()));
		if(i<11) tf4[i]=TFile::Open(Form("root_files_cat_DNtuple/%s.root",st[i].data()));
		else tf4[i]=TFile::Open(Form("root_files_cat_newNtuple/%s.root",st[i].data()));
	}
	
	TFile* tfM=TFile::Open("/home/nstar/HH4b/BumpSearchLimitCode/input/MassPlotFineBins.root");
	
	double eff1[11],eff2[11],eff3[11],eff4[11];
	double x[11]={1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500};
	TH1D *thSh[12], *thSh2[12], *thSh3[12],*thShM[12],*th2[12];
	
	for(int j=0;j<5;j++){
	
	for(int i=0;i<11;i++){
			thSh[i]=(TH1D*)tf4[i]->FindObjectAny(Form("cat%d",j));
			eff1[i]=thSh[i]->GetEntries();
			eff3[i]=thSh[i]->GetEntries()/50000;
			
			thSh2[i]=(TH1D*)tf1[i]->FindObjectAny(Form("cat%d",j));
			eff2[i]=thSh2[i]->GetEntries();
			eff4[i]=thSh2[i]->GetEntries()/50000;
	}
	TGraph* tgEff1=new TGraph(11,x,eff1);
	TGraph* tgEff2=new TGraph(11,x,eff2);
	TGraph* tgEff3=new TGraph(11,x,eff3);
	TGraph* tgEff4=new TGraph(11,x,eff4);
	
	tgEff2->SetLineColor(2);
	tgEff4->SetLineColor(2);
	
	tgEff1->SetTitle(Form("%db-tag",j));
	tgEff2->SetFillColor(kWhite);
	tgEff1->SetFillColor(kWhite);
	
	tgEff3->SetTitle(Form("%db-tag eff",j));
	tgEff4->SetFillColor(kWhite);
	tgEff3->SetFillColor(kWhite);
	
	TLegend* legEff ;
		legEff=new TLegend(0.791452,0.562447,0.890645,0.783966);
		legEff->AddEntry(tgEff1,"Maxime's cut");
		legEff->AddEntry(tgEff2,"origin cut");
	tgEff1->Draw("APL");
	tgEff2->Draw("PL,same");
	legEff->Draw("same");
	
	if(j==0)c1->Print("0115_eff.pdf(");
	//else if (j==4)c1->Print("0115_eff.pdf)");
		else c1->Print("0115_eff.pdf");
	
	tgEff3->Draw("APL");
	tgEff4->Draw("PL,same");
	legEff->Draw("same");
	
	//if(j==0)c1->Print("0115_eff.pdf(");
	if (j==4)c1->Print("0115_eff.pdf)");
		else c1->Print("0115_eff.pdf");
	}
	
	
	double xx[10]={1200,1400,1600,1800,2000,2500,3000,3500,4000,4500};
	double y1[11],y2[11],y3[11];
	double yy1[10],yy2[10],yy3[10];
	for(int i=0;i<12;i++){
		thSh[i]=(TH1D*)tf1[i]->FindObjectAny("cat4");
		thSh2[i]=(TH1D*)tf2[i]->FindObjectAny("cat4");
		thSh3[i]=(TH1D*)tf3[i]->FindObjectAny("cat4");
		
		
		
		if(i<11){
			
			//y1[i]=thSh[i]->GetMean();
			//y2[i]=thSh2[i]->GetMean();
			//y3[i]=thSh3[i]->GetMean();
			
			y1[i]=thSh[i]->GetMean()/x[i];
			y2[i]=thSh2[i]->GetMean()/x[i];
			y3[i]=thSh3[i]->GetMean()/x[i];
			
			double temp=y1[i];
			y1[i]/=temp;
			y2[i]/=temp;
			y3[i]/=temp;
			
			
		}
		
		
		if(i<11) thSh[i]->SetName(Form("Radion%s_cat0",st[i].data()));
		else {
			thSh[i]->SetName("QCD_cat0");
			thSh[i]->Sumw2();
			thSh[i]->Scale(1.96/2.7);
		}
		th2[i]=(TH1D*)tf1[i]->FindObjectAny("cat3");
		if(i<11) th2[i]->SetName(Form("Radion%s_cat1",st[i].data()));
		else {
			th2[i]->SetName("QCD_cat1");
			th2[i]->Sumw2();
			th2[i]->Scale(1.96/2.7);
		}
		
		
		if(i!=7){
			
		if(i<11)thShM[i]=(TH1D*)tfM->FindObjectAny(Form("Radion%s_cat0",st[i].data()));
		else thShM[i]=(TH1D*)tfM->FindObjectAny("QCD_cat0");
		cout<<Form("Radion%s_cat0",st[i].data())<<endl;
		
		thSh[i]->Rebin(5);
		thShM[i]->Rebin(5);
		
		thSh[i]->Sumw2();
		thShM[i]->Sumw2();
		
		
		thSh[i]->Scale(1/thSh[i]->Integral());
		thShM[i]->Scale(1/thShM[i]->Integral());
		
		thShM[i]->GetXaxis()->SetRangeUser(x[i]*0.8>1000?x[i]*0.8:1000,x[i]*1.25);
		
		thShM[i]->SetLineColor(2);
		TLegend* leg ;
		leg=new TLegend(0.791452,0.562447,0.890645,0.783966);
		leg->AddEntry(thSh[i],"Ching-Wei's");
		leg->AddEntry(thShM[i],"Maxime's");
		thShM[i]->Draw();
		thSh[i]->Draw("same");
		leg->Draw("same");
		
		if(i==0) c1->Print("1228_3.pdf(");
		else if (i==11) c1->Print("1228_3.pdf)");
		else  c1->Print("1228_3.pdf");
		}
		
		
		
	}
	//cout<<"check"<<endl;
	
	for (int i=0;i<10;i++){
		yy1[i]=y1[i+1];
		yy2[i]=y2[i+1];
		yy3[i]=y3[i+1];
		
	}
	
	TGraph * tg1=new TGraph(10,xx,yy1);
	TGraph * tg2=new TGraph(10,xx,yy2);
	TGraph * tg3=new TGraph(10,xx,yy3);
	//TGraph * tg1=new TGraph(11,x,y1);
	//TGraph * tg2=new TGraph(11,x,y2);
	//TGraph * tg3=new TGraph(11,x,y3);
	
	
	tg1->SetLineColor(1);
	tg2->SetLineColor(2);
	tg3->SetLineColor(3);
	
	tg1->SetFillColor(kWhite);
	tg2->SetFillColor(kWhite);
	tg3->SetFillColor(kWhite);
	
	tg1->GetXaxis()->SetTitle("M_{X}[GeV]");
	tg1->GetYaxis()->SetTitle("Mean/M_{X}[GeV](normalized)");	
	tg1->SetMaximum(1.05);
	tg1->Draw("APL");
	tg2->Draw("PL,same");
	tg3->Draw("PL,same");
	
	TLegend* leg ;
	leg=new TLegend(0.191452,0.562447,0.290645,0.783966);
	leg->AddEntry(tg1,"JES");
	leg->AddEntry(tg2,"JES_UP");
	leg->AddEntry(tg3,"JES_DO");
	
	leg->Draw("same");
	
	c1->Print("0115_1.pdf");
	
	TF1 *fa1 = new TF1("fa1","[0]+[1]*x",1200,4500); 
	TF1 *fa2 = new TF1("fa1","[0]+[1]*x",1200,4500); 
	tg2->Fit(fa1);
	tg2->Draw("APL");
	fa1->Draw("same");
	c1->Print("0115_2.pdf");
	
	tg3->Fit(fa2);
	tg3->Draw("APL");
	fa2->Draw("same");
	c1->Print("0115_3.pdf");
	
	TFile* outFile = new TFile("MassPlotFineBins_5.root","recreate");
	
	
	for(int i=0;i<12;i++){
		thSh[i]->Write();
		th2[i]->Write();
	}
	
	outFile->Close();
}

