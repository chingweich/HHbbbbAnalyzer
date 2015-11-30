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
#include "untuplizer.h"
#include <TClonesArray.h>
#include <fstream>
#include <cmath>
#include <TSystem.h>
#include <string>
#include <sstream>
#include "setNCUStyle.C"
#include<TH2.h>
#include "TLine.h"

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

void HHbbbbDrawerNum(int a){
	setNCUStyle();
	c1 = new TCanvas("c1","",1360,768);
	
	TFile* tf1[11];
	string st[11]={"1000","1200","1400","1600","2000","2500","3000","3500","4000","4500","QCD"};
	string st3[11]={"SIG_M1000","SIG_M1200","SIG_M1400","SIG_M1600","SIG_M2000","SIG_M2500","SIG_M3000","SIG_M3500","SIG_M4000","SIG_M4500","QCD"};
	int st2[10]={1000,1200,1400,1600,2000,2500,3000,3500,4000,4500};
	for(int i=0;i<11;i++)tf1[i]=TFile::Open(Form("root_files_Op_DeltaEta_shape/%s.root",st[i].data()));
	
	
	
	//TH2F *th2f1,*th2f2;
	for(int i=0;i<10;i++){
		//TH1D * th1 = (TH1D*)tf1[i]->FindObjectAny(Form("%s_dEta",st[i].data()));
		//TH1D * th2 = (TH1D*)tf1[10]->FindObjectAny(Form("%s_dEta",st[i].data()));
		/*
		TLegend *leg = new TLegend(.75, .25, .95, .35);
		leg->AddEntry(th1,"Signal");
		leg->AddEntry(th2,"Bkg.");
		th2->SetLineColor(2);
		th1->Rebin(20);
		th1->SetTitle(Form("M_{X}=%s[GeV]",st[i].data()));
		th1->SetXTitle("deltaEta");
		th2->Rebin(20);
		th1->Scale(1/th1->Integral());
		th2->Scale(1/th2->Integral());
		th1->Draw();
		th2->Draw("same");
		leg->Draw("same");
		if(i==0)c1->Print("pdf/dEta.pdf(");
		else if (i==9)c1->Print("pdf/dEta.pdf)");
		else c1->Print("pdf/dEta.pdf");
		//cout<<st[i].data()<<"="<<th2f1->Integral(a,b,c,d)<<endl;
		//cout<<st[i].data()<<"="<<th2f2->Integral(a,b,c,d)<<endl;
		
		
		ofstream myfile;
		myfile.open (Form("opt_%d/cout_%d.txt",a,st2[i]));
		myfile<<"DATA 0"<<endl;
		myfile<<"QCD "<<th2->Integral(1,a)<<endl;
		for(int j=0;j<10;j++)myfile<<"M"<<st2[j]<<" "<<th1->Integral(1,a)<<endl;
		myfile.close();
		*/
	}
	//int numSh[8]={8,9,10,11,12,13,14,15};
	
	TH1D * thSh[11];
	for(int i=0;i<11;i++){
		thSh[i]=(TH1D*)tf1[i]->FindObjectAny(Form("Mass_dEta%d",a));
		thSh[i]->SetName(Form("%s",st3[i].data()));
	}
	//cout<<"check"<<endl;
	
	
	ofstream myfile;
	myfile.open (Form("opt_%d/cout.txt",a));
	myfile<<"DATA 0"<<endl;
	myfile<<"QCD "<<thSh[10]->Integral()<<endl;
	for(int j=0;j<10;j++)myfile<<"M"<<st2[j]<<" "<<thSh[j]->Integral()<<endl;
	myfile.close();
	
	TFile* outFile = new TFile(Form("opt_%d/shape.root",a),"recreate");
	
	TH1D * thD=new TH1D("data_obs","Mass_dEta8",1000,0,5000);
	thD->Write();
	for(int i=0;i<11;i++)thSh[i]->Write();
	outFile->Close();
	
	thSh[0]->SetLineColor(1);
	thSh[3]->SetLineColor(2);
	thSh[6]->SetLineColor(3);
	thSh[9]->SetLineColor(4);
	thSh[10]->SetLineColor(5);
	
	for(int i=0;i<11;i++){
		thSh[i]->Scale(1/thSh[i]->Integral());
		thSh[i]->Rebin(5);
	}
	
	TLegend *leg = new TLegend(.75, .75, .95, .95);
		leg->AddEntry(thSh[0],"M1000");
		leg->AddEntry(thSh[3],"M1600");
		leg->AddEntry(thSh[6],"M3000");
		leg->AddEntry(thSh[9],"M4500");
		leg->AddEntry(thSh[10],"QCD");
	
	thSh[0]->Draw();
	thSh[3]->Draw("same");
	thSh[6]->Draw("same");
	thSh[9]->Draw("same");
	thSh[10]->Draw("same");
	leg->Draw("same");
	c1->SetLogy();
	c1->Print(Form("opt_%d/Mjj.png",a));
	
	
}