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

void draw(){
	setNCUStyle();
	c1 = new TCanvas("c1","",1360,768);
	string  masspoint[11]={"1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500"};
	TFile *f[40];
	for(int i=0;i<11;i++){
		f[i]=TFile::Open(Form("B%s_2.root",masspoint[i].data()));
		f[i+11]=TFile::Open(Form("R%s_2.root",masspoint[i].data()));
	}
	
	f[22]=TFile::Open("QCD700_2.root");
	f[23]=TFile::Open("QCD1000_2.root");
	f[24]=TFile::Open("QCD1500_2.root");
	f[25]=TFile::Open("QCD2000_2.root");
	f[26]=TFile::Open("../root_files_btaggedScaleFactor/data.root");
	f[27]=TFile::Open("data.root");
	
	TFile *f2[40];
	for(int i=0;i<11;i++){
		f2[i]=TFile::Open(Form("B%s.root",masspoint[i].data()));
		f2[i+11]=TFile::Open(Form("R%s.root",masspoint[i].data()));
	}
	
	f2[22]=TFile::Open("QCD700.root");
	f2[23]=TFile::Open("QCD1000.root");
	f2[24]=TFile::Open("QCD1500.root");
	f2[25]=TFile::Open("QCD2000.root");
	f2[26]=TFile::Open("data.root");
	f2[27]=TFile::Open("test.root");
	
	
	TH1D* sf2d[10];
	sf2d[0]=(TH1D*)f2[22]->FindObjectAny("effD_l_1d");
	sf2d[1]=(TH1D*)f2[23]->FindObjectAny("effD_l_1d");
	sf2d[2]=(TH1D*)f2[24]->FindObjectAny("effD_l_1d");
	sf2d[3]=(TH1D*)f2[25]->FindObjectAny("effD_l_1d");
	sf2d[4]=(TH1D*)f2[26]->FindObjectAny("effD_l_1d");
	sf2d[5]=(TH1D*)f2[22]->FindObjectAny("effN_l_1d");
	sf2d[6]=(TH1D*)f2[23]->FindObjectAny("effN_l_1d");
	sf2d[7]=(TH1D*)f2[24]->FindObjectAny("effN_l_1d");
	sf2d[8]=(TH1D*)f2[25]->FindObjectAny("effN_l_1d");
	sf2d[9]=(TH1D*)f2[26]->FindObjectAny("effN_l_1d");
	
	
	TH1D* sf2d2[4];
	sf2d2[0]=(TH1D*)sf2d[0]->Clone();
	sf2d2[0]->Add(sf2d[1]);
	sf2d2[0]->Add(sf2d[2]);
	sf2d2[0]->Add(sf2d[3]);
	
	sf2d2[1]=(TH1D*)sf2d[5]->Clone();
	sf2d2[1]->Add(sf2d[6]);
	sf2d2[1]->Add(sf2d[7]);
	sf2d2[1]->Add(sf2d[8]);
	
	TH1D* thh=(TH1D*)f2[27]->FindObjectAny("sf1d");
	
	
	//sf2d2[1]->RebinX(5);
	//sf2d2[0]->RebinX(5);
	//sf2d[9]->RebinX(5);
	//sf2d[4]->RebinX(5);
	
	//sf2d2[1]->RebinY(2);
	//sf2d2[0]->RebinY(2);
	//sf2d[9]->RebinY(2);
	//sf2d[4]->RebinY(2);
	
	double err1[200],err2[200];
	
	for (int i=0;i<200;i++){
		//cout<<5+i*10<<","<<pow(sf2d2[1]->GetBinError(i+1)/sf2d2[1]->GetBinContent(i+1),2)<<","<<pow(sf2d2[0]->GetBinError(i+1)/sf2d2[0]->GetBinContent(i+1),2);
		if      (sf2d2[1]->GetBinContent(i+1)==0 && sf2d2[0]->GetBinContent(i+1)!=0)err1[i]=(sqrt(pow(sf2d2[0]->GetBinError(i+1)/sf2d2[0]->GetBinContent(i+1),2)));
		else if (sf2d2[0]->GetBinContent(i+1)==0 && sf2d2[1]->GetBinContent(i+1)!=0)err1[i]=(sqrt(pow(sf2d2[1]->GetBinError(i+1)/sf2d2[1]->GetBinContent(i+1),2)));
		else if (sf2d2[1]->GetBinContent(i+1)==0 && sf2d2[0]->GetBinContent(i+1)==0)err1[i]=0;
		else err1[i]=(sqrt(pow(sf2d2[1]->GetBinError(i+1)/sf2d2[1]->GetBinContent(i+1),2)+pow(sf2d2[0]->GetBinError(i+1)/sf2d2[0]->GetBinContent(i+1),2)));
		
		if      (sf2d[9]->GetBinContent(i+1)==0 && sf2d[4]->GetBinContent(i+1)!=0)err2[i]=(sqrt(pow(sf2d[4]->GetBinError(i+1)/sf2d[4]->GetBinContent(i+1),2)));
		else if (sf2d[9]->GetBinContent(i+1)!=0 && sf2d[4]->GetBinContent(i+1)==0)err2[i]=(sqrt(pow(sf2d[9]->GetBinError(i+1)/sf2d[9]->GetBinContent(i+1),2)));
		else if (sf2d[9]->GetBinContent(i+1)==0 && sf2d[4]->GetBinContent(i+1)==0)err2[i]=0;
		else err2[i]=(sqrt(pow(sf2d[9]->GetBinError(i+1)/sf2d[9]->GetBinContent(i+1),2)+pow(sf2d[4]->GetBinError(i+1)/sf2d[4]->GetBinContent(i+1),2)));
		cout<<"1="<<err1[i]<<",2="<<err2[i]<<",,";
		if (sf2d2[0]->GetBinContent(i+1)!=0)err1[i]*=(sf2d2[1]->GetBinContent(i+1)/sf2d2[0]->GetBinContent(i+1));
		else err1[i]=0;
		if (sf2d[4]->GetBinContent(i+1)!=0)err2[i]*=(sf2d[9]->GetBinContent(i+1)/sf2d[4]->GetBinContent(i+1));
		else err2[i]=0;
		cout<<"1="<<err1[i]<<",2="<<err2[i]<<endl;
		//sf2d2[1]->SetBinError(i+1,err1[i]);
		//sf2d[9]->SetBinError(i+1,err2[i]);
	}
	
	sf2d2[1]->Divide(sf2d2[0]);
	sf2d[9]->Divide(sf2d[4]);
	
	for (int i=0;i<200;i++){
	
		sf2d2[1]->SetBinError(i+1,err1[i]);
		sf2d[9]->SetBinError(i+1,err2[i]);
	}
	
	sf2d[9]->SetXTitle("Pt(GeV)");
	sf2d[1]->SetXTitle("Pt(GeV)");
	
	sf2d[9]->GetXaxis()->SetRangeUser(0,1000);
	sf2d[9]->Draw("e");
	c1->Print("sf/eff_data_1d.pdf");
	
	sf2d2[1]->GetXaxis()->SetRangeUser(0,1000);
	sf2d2[1]->Draw("e");
	c1->Print("sf/eff_QCD_1d.pdf");
	
	
	for (int i=0;i<200;i++){
		//cout<<5+i*10<<","<<pow(sf2d2[1]->GetBinError(i+1)/sf2d2[1]->GetBinContent(i+1),2)<<","<<pow(sf2d2[0]->GetBinError(i+1)/sf2d2[0]->GetBinContent(i+1),2);
		if      (sf2d2[1]->GetBinContent(i+1)!=0 && sf2d[9]->GetBinContent(i+1)==0)err1[i]=(sqrt(pow(sf2d2[1]->GetBinError(i+1)/sf2d2[1]->GetBinContent(i+1),2)));
		else if (sf2d2[1]->GetBinContent(i+1)==0 && sf2d[9]->GetBinContent(i+1)!=0)err1[i]=(sqrt(pow(sf2d[9]->GetBinError(i+1)/sf2d[9]->GetBinContent(i+1),2)));
		else if (sf2d2[1]->GetBinContent(i+1)==0 && sf2d[9]->GetBinContent(i+1)==0)err1[i]=0;
		else err1[i]=(sqrt(pow(sf2d2[1]->GetBinError(i+1)/sf2d2[1]->GetBinContent(i+1),2)+pow(sf2d[9]->GetBinError(i+1)/sf2d[9]->GetBinContent(i+1),2)));
		
		if (sf2d2[0]->GetBinContent(i+1)!=0)err1[i]*=(sf2d2[1]->GetBinContent(i+1)/sf2d2[0]->GetBinContent(i+1));
		else err1[i]=0;
		//sf2d[9]->SetBinError(i+1,err1[i]);
		
	}
	
	
	sf2d[9]->Divide(sf2d2[1]);
	for (int i=0;i<200;i++)sf2d[9]->SetBinError(i+1,err1[i]);
	
	TLegend* leg5 ;																																						
	leg5=new TLegend(0.451452,0.662447,0.670645,0.883966);
	leg5->AddEntry(sf2d[9],"derived");
	leg5->AddEntry(thh,"csv file");
	
	sf2d[9]->SetYTitle("scale factor");
	sf2d[9]->Draw("e");
	
	thh->SetLineColor(2);
	thh->Draw("same");
	leg5->Draw("same");
	c1->Print("sf/sf_1d.pdf");
	
	string QCDName[4]={"700","1000","1500","2000"};
	

	
	for(int i=0;i<26;i++){
		TH1D * thz,* thx;
		thz=(TH1D*)f2[i]->FindObjectAny("effD_b_1d");
		thx=(TH1D*)f2[i]->FindObjectAny("effN_b_1d");
		
		thz->Rebin(10);
		thx->Rebin(10);
		TGraphAsymmErrors* tg=new TGraphAsymmErrors(thx,thz);
		tg->Draw("APL");
		if(i<11)c1->Print(Form("plot/B%s/eff_b_1d.pdf",masspoint[i].data()));
		else if (i<22 && i>10)c1->Print(Form("plot/R%s/eff_b_1d.pdf",masspoint[i-11].data()));
		else c1->Print(Form("plot/QCD%s/eff_b_1d.pdf",QCDName[i-22].data()));
		
		thz=(TH1D*)f2[i]->FindObjectAny("effD_c_1d");
		thx=(TH1D*)f2[i]->FindObjectAny("effN_c_1d");
		thz->Rebin(10);
		thx->Rebin(10);
		TGraphAsymmErrors* tg2=new TGraphAsymmErrors(thx,thz);
		tg2->Draw("APL");
		if(i<11)c1->Print(Form("plot/B%s/eff_c_1d.pdf",masspoint[i].data()));
		else if (i<22 && i>10)c1->Print(Form("plot/R%s/eff_c_1d.pdf",masspoint[i-11].data()));
		else c1->Print(Form("plot/QCD%s/eff_c_1d.pdf",QCDName[i-22].data()));
		
		thz=(TH1D*)f2[i]->FindObjectAny("effD_l_1d");
		thx=(TH1D*)f2[i]->FindObjectAny("effN_l_1d");
		thz->Rebin(10);
		thx->Rebin(10);
		TGraphAsymmErrors* tg3=new TGraphAsymmErrors(thx,thz);
		tg3->Draw("APL");
		if(i<11)c1->Print(Form("plot/B%s/eff_l_1d.pdf",masspoint[i].data()));
		else if (i<22 && i>10)c1->Print(Form("plot/R%s/eff_l_1d.pdf",masspoint[i-11].data()));
		else c1->Print(Form("plot/QCD%s/eff_l_1d.pdf",QCDName[i-22].data()));
		
	}
	string comparingName[51]={
		"Pt_j0_sj0_0b","Pt_j0_sj1_0b","Pt_j1_sj0_0b","Pt_j1_sj1_0b",
		"Pt_j0_sj0_1b","Pt_j0_sj1_1b","Pt_j1_sj0_1b","Pt_j1_sj1_1b",
		"Pt_j0_sj0_2b","Pt_j0_sj1_2b","Pt_j1_sj0_2b","Pt_j1_sj1_2b",
		"deltaR0_0b","deltaR1_0b","deltaR0_1b","deltaR1_1b","deltaR0_2b","deltaR1_2b",
		"Eta_j0_sj0_0b","Eta_j0_sj1_0b","Eta_j1_sj0_0b","Eta_j1_sj1_0b",
		"Eta_j0_sj0_1b","Eta_j0_sj1_1b","Eta_j1_sj0_1b","Eta_j1_sj1_1b",
		"Eta_j0_sj0_2b","Eta_j0_sj1_2b","Eta_j1_sj0_2b","Eta_j1_sj1_2b",
		"Pt_j0_0b","Pt_j1_0b","Pt_j0_1b","Pt_j1_1b","Pt_j0_2b","Pt_j1_2b",
		"Eta_j0_0b","Eta_j1_0b","Eta_j0_1b","Eta_j1_1b","Eta_j0_2b","Eta_j1_2b",
		"totalMass_0b","totalMass_1b","totalMass_2b",
		"prMass_j0_0b","prMass_j1_0b","prMass_j0_1b","prMass_j1_1b","prMass_j0_2b","prMass_j1_2b"
	};
	
	for(int i=0;i<26;i++){
		
		TH2D* th2d=(TH2D*)f[i]->FindObjectAny("SF_vs_Pt_b");
		th2d->Draw("colz");
		if(i<11)c1->Print(Form("plot/B%s/SF_vs_Pt_b.pdf",masspoint[i].data()));
		else if (i<22 && i>10)c1->Print(Form("plot/R%s/SF_vs_Pt_b.pdf",masspoint[i-11].data()));
		else c1->Print(Form("plot/QCD%s/SF_vs_Pt_b.pdf",QCDName[i-22].data()));
		
		th2d=(TH2D*)f[i]->FindObjectAny("SF_vs_Pt_c");
		th2d->Draw("colz");
		if(i<11)c1->Print(Form("plot/B%s/SF_vs_Pt_c.pdf",masspoint[i].data()));
		else if (i<22 && i>10)c1->Print(Form("plot/R%s/SF_vs_Pt_c.pdf",masspoint[i-11].data()));
		else c1->Print(Form("plot/QCD%s/SF_vs_Pt_c.pdf",QCDName[i-22].data()));
		
		th2d=(TH2D*)f[i]->FindObjectAny("SF_vs_Pt_l");
		th2d->Draw("colz");
		if(i<11)c1->Print(Form("plot/B%s/SF_vs_Pt_l.pdf",masspoint[i].data()));
		else if (i<22 && i>10)c1->Print(Form("plot/R%s/SF_vs_Pt_l.pdf",masspoint[i-11].data()));
		else c1->Print(Form("plot/QCD%s/SF_vs_Pt_l.pdf",QCDName[i-22].data()));
		
		TH1D* th3;
		/*
		for(int j=0;j<51;j++){
			TH1D* th5[2];
			th5[0]=(TH1D*)f[i]->FindObjectAny(Form("%s",comparingName[j].data()));
			th5[1]=(TH1D*)f[i]->FindObjectAny(Form("%ss",comparingName[j].data()));
			th5[1]->SetLineColor(2);
			th5[0]->Draw();
			th5[1]->Draw("same");
			TLegend* leg3 ;
			leg3=new TLegend(0.691452,0.662447,0.890645,0.883966);
			leg3->AddEntry(th5[0],"withoutSF");
			leg3->AddEntry(th5[1],"withSF");
			leg3->Draw("same");
			if(i<11)c1->Print(Form("plot/B%s/%s.pdf",masspoint[i].data(),comparingName[j].data()));
			else if (i<22 && i>10)c1->Print(Form("plot/R%s/%s.pdf",masspoint[i-11].data(),comparingName[j].data()));
			else c1->Print(Form("plot/QCD%s/%s.pdf",QCDName[i-22].data(),comparingName[j].data()));
		}
		*/
		
		
		c1->SetLogy(1);
		th3=(TH1D*)f[i]->FindObjectAny("SF_jet0_sub0_pass");
		th3->Draw();
		if(i<11)c1->Print(Form("plot/B%s/.pdf",masspoint[i].data()));
		else if (i<22 && i>10)c1->Print(Form("plot/R%s/SF_jet0_sub0_pass.pdf",masspoint[i-11].data()));
		else c1->Print(Form("plot/QCD%s/SF_jet0_sub0_pass.pdf",QCDName[i-22].data()));
		th3=(TH1D*)f[i]->FindObjectAny("SF_jet0_sub1_pass");
		th3->Draw();
		if(i<11)c1->Print(Form("plot/B%s/SF_jet0_sub1_pass.pdf",masspoint[i].data()));
		else if (i<22&& i>10)c1->Print(Form("plot/R%s/SF_jet0_sub1_pass.pdf",masspoint[i-11].data()));
		else c1->Print(Form("plot/QCD%s/SF_jet0_sub1_pass.pdf",QCDName[i-22].data()));
		th3=(TH1D*)f[i]->FindObjectAny("SF_jet1_sub0_pass");
		th3->Draw();
		if(i<11)c1->Print(Form("plot/B%s/SF_jet1_sub0_pass.pdf",masspoint[i].data()));
		else if (i<22&& i>10)c1->Print(Form("plot/R%s/SF_jet1_sub0_pass.pdf",masspoint[i-11].data()));
		else c1->Print(Form("plot/QCD%s/SF_jet1_sub0_pass.pdf",QCDName[i-22].data()));
		th3=(TH1D*)f[i]->FindObjectAny("SF_jet1_sub1_pass");
		th3->Draw();
		if(i<11)c1->Print(Form("plot/B%s/SF_jet1_sub1_pass.pdf",masspoint[i].data()));
		else if (i<22&& i>10)c1->Print(Form("plot/R%s/SF_jet1_sub1_pass.pdf",masspoint[i-11].data()));
		else c1->Print(Form("plot/QCD%s/SF_jet1_sub1_pass.pdf",QCDName[i-22].data()));
		th3=(TH1D*)f[i]->FindObjectAny("SF_jet0_sub0_fail");
		th3->Draw();
		if(i<11)c1->Print(Form("plot/B%s/SF_jet0_sub0_fail.pdf",masspoint[i].data()));
		else if (i<22&& i>10)c1->Print(Form("plot/R%s/SF_jet0_sub0_fail.pdf",masspoint[i-11].data()));
		else c1->Print(Form("plot/QCD%s/SF_jet0_sub0_fail.pdf",QCDName[i-22].data()));
		th3=(TH1D*)f[i]->FindObjectAny("SF_jet0_sub1_fail");
		th3->Draw();
		if(i<11)c1->Print(Form("plot/B%s/SF_jet0_sub1_fail.pdf",masspoint[i].data()));
		else if (i<22&& i>10)c1->Print(Form("plot/R%s/SF_jet0_sub1_fail.pdf",masspoint[i-11].data()));
		else c1->Print(Form("plot/QCD%s/SF_jet0_sub1_fail.pdf",QCDName[i-22].data()));
		th3=(TH1D*)f[i]->FindObjectAny("SF_jet1_sub0_fail");
		th3->Draw();
		if(i<11)c1->Print(Form("plot/B%s/SF_jet1_sub0_fail.pdf",masspoint[i].data()));
		else if (i<22&& i>10)c1->Print(Form("plot/R%s/SF_jet1_sub0_fail.pdf",masspoint[i-11].data()));
		else c1->Print(Form("plot/QCD%s/SF_jet1_sub0_fail.pdf",QCDName[i-22].data()));
		th3=(TH1D*)f[i]->FindObjectAny("SF_jet1_sub1_fail");
		th3->Draw();
		if(i<11)c1->Print(Form("plot/B%s/SF_jet1_sub1_fail.pdf",masspoint[i].data()));
		else if (i<22&& i>10)c1->Print(Form("plot/R%s/SF_jet1_sub1_fail.pdf",masspoint[i-11].data()));
		else c1->Print(Form("plot/QCD%s/SF_jet1_sub1_fail.pdf",QCDName[i-22].data()));
		th3=(TH1D*)f[i]->FindObjectAny("weight");
		th3->Draw();
		if(i<11)c1->Print(Form("plot/B%s/weight.pdf",masspoint[i].data()));
		else if (i<22&& i>10)c1->Print(Form("plot/R%s/weight.pdf",masspoint[i-11].data()));
		else c1->Print(Form("plot/QCD%s/weight.pdf",QCDName[i-22].data()));
		
		th3=(TH1D*)f[i]->FindObjectAny("weight_0b");
		th3->Draw();
		if(i<11)c1->Print(Form("plot/B%s/weight_0b.pdf",masspoint[i].data()));
		else if (i<22&& i>10)c1->Print(Form("plot/R%s/weight_0b.pdf",masspoint[i-11].data()));
		else c1->Print(Form("plot/QCD%s/weight_0b.pdf",QCDName[i-22].data()));
		
		th3=(TH1D*)f[i]->FindObjectAny("weight_1b");
		th3->Draw();
		if(i<11)c1->Print(Form("plot/B%s/weight_1b.pdf",masspoint[i].data()));
		else if (i<22&& i>10)c1->Print(Form("plot/R%s/weight_1b.pdf",masspoint[i-11].data()));
		else c1->Print(Form("plot/QCD%s/weight_1b.pdf",QCDName[i-22].data()));
		
		th3=(TH1D*)f[i]->FindObjectAny("weight_2b");
		th3->Draw();
		if(i<11)c1->Print(Form("plot/B%s/weight_2b.pdf",masspoint[i].data()));
		else if (i<22&& i>10)c1->Print(Form("plot/R%s/weight_2b.pdf",masspoint[i-11].data()));
		else c1->Print(Form("plot/QCD%s/weight_2b.pdf",QCDName[i-22].data()));
		
		TH1D* th4=(TH1D*)f[i]->FindObjectAny("weight_2b_ll");
		TH1D* th5=(TH1D*)f[i]->FindObjectAny("weight_2b_onel");
		th4->SetLineColor(2);
		th5->SetLineColor(3);
		th4->Draw("same");
		th5->Draw("same");
		TLegend* leg ;
		leg=new TLegend(0.691452,0.662447,0.890645,0.883966);
		leg->AddEntry(th4,"twoLightJet");
		leg->AddEntry(th5,"oneLightJet");	
		leg->Draw("same");	
		if(i<11)c1->Print(Form("plot/B%s/weight_2b.pdf",masspoint[i].data()));
		else if (i<22&& i>10)c1->Print(Form("plot/R%s/weight_2b.pdf",masspoint[i-11].data()));
		else c1->Print(Form("plot/QCD%s/weight_2b.pdf",QCDName[i-22].data()));		
		
		c1->SetLogy(0);
	}
	
	
	string st[52]={
		"Nbtagjet","Pt_j0_sj0_0b","Pt_j0_sj1_0b","Pt_j1_sj0_0b","Pt_j1_sj1_0b",
		"Pt_j0_sj0_1b","Pt_j0_sj1_1b","Pt_j1_sj0_1b","Pt_j1_sj1_1b",
		"Pt_j0_sj0_2b","Pt_j0_sj1_2b","Pt_j1_sj0_2b","Pt_j1_sj1_2b",
		"deltaR0_0b","deltaR1_0b","deltaR0_1b","deltaR1_1b","deltaR0_2b","deltaR1_2b",
		"Eta_j0_sj0_0b","Eta_j0_sj1_0b","Eta_j1_sj0_0b","Eta_j1_sj1_0b",
		"Eta_j0_sj0_1b","Eta_j0_sj1_1b","Eta_j1_sj0_1b","Eta_j1_sj1_1b",
		"Eta_j0_sj0_2b","Eta_j0_sj1_2b","Eta_j1_sj0_2b","Eta_j1_sj1_2b",
		"Pt_j0_0b","Pt_j1_0b","Pt_j0_1b","Pt_j1_1b","Pt_j0_2b","Pt_j1_2b",
		"Eta_j0_0b","Eta_j1_0b","Eta_j0_1b","Eta_j1_1b","Eta_j0_2b","Eta_j1_2b",
		"totalMass_0b","totalMass_1b","totalMass_2b",
		"prMass_j0_0b","prMass_j1_0b","prMass_j0_1b","prMass_j1_1b","prMass_j0_2b","prMass_j1_2b"
	};
	
	TH1D* th1[9][36];
	
	
	for(int j=0;j<52;j++){
		th1[0][j]=(TH1D*)f[1]->FindObjectAny((Form("%s",st[j].data())));
		th1[1][j]=(TH1D*)f[4]->FindObjectAny((Form("%s",st[j].data())));
		th1[2][j]=(TH1D*)f[6]->FindObjectAny((Form("%s",st[j].data())));
		th1[3][j]=(TH1D*)f[15]->FindObjectAny((Form("%s",st[j].data())));
		//th1[4][j]=(TH1D*)f[22]->FindObjectAny((Form("%s",st[j].data())));
		//th1[5][j]=(TH1D*)f[23]->FindObjectAny((Form("%s",st[j].data())));
		//th1[6][j]=(TH1D*)f[24]->FindObjectAny((Form("%s",st[j].data())));
		//th1[7][j]=(TH1D*)f[25]->FindObjectAny((Form("%s",st[j].data())));
		if(j==0)th1[8][j]=(TH1D*)f[26]->FindObjectAny((Form("%s",st[j].data())));
		else th1[8][j]=(TH1D*)f[27]->FindObjectAny((Form("%s",st[j].data())));
		
		
		
		TH1D* th3[4];
		th3[0]=(TH1D*)f[22]->FindObjectAny((Form("%s",st[j].data())));
		th3[1]=(TH1D*)f[23]->FindObjectAny((Form("%s",st[j].data())));
		th3[2]=(TH1D*)f[24]->FindObjectAny((Form("%s",st[j].data())));
		th3[3]=(TH1D*)f[25]->FindObjectAny((Form("%s",st[j].data())));
		
		if(j==0){
		th3[0]=(TH1D*)f[22]->FindObjectAny((Form("%sS",st[j].data())));
		th3[1]=(TH1D*)f[23]->FindObjectAny((Form("%sS",st[j].data())));
		th3[2]=(TH1D*)f[24]->FindObjectAny((Form("%sS",st[j].data())));
		th3[3]=(TH1D*)f[25]->FindObjectAny((Form("%sS",st[j].data())));
		}
		
		for(int ii=0;ii<4;ii++){
			th3[ii]->Sumw2();
			th3[ii]->Scale(7902/11079.8);
		}
		
		
		TH1D* th6[4];
		TH1D* th7[4];
		if(j!=0){
		th6[0]=(TH1D*)f[22]->FindObjectAny((Form("%ss",st[j].data())));
		th6[1]=(TH1D*)f[23]->FindObjectAny((Form("%ss",st[j].data())));
		th6[2]=(TH1D*)f[24]->FindObjectAny((Form("%ss",st[j].data())));
		th6[3]=(TH1D*)f[25]->FindObjectAny((Form("%ss",st[j].data())));
		for(int ii=0;ii<4;ii++){
			th6[ii]->Sumw2();
			th6[ii]->Scale(7902/11079.8);
		}
		th7[3]=(TH1D*)th6[3]->Clone();
		th7[3]->Add(th6[0]);
		th7[3]->Add(th6[1]);
		th7[3]->Add(th6[2]);
		
		}
		
		
		
		
		cout<<"check"<<endl;
		th1[4][j]=(TH1D*)th3[0]->Clone();
		th1[5][j]=(TH1D*)th3[1]->Clone();
		th1[6][j]=(TH1D*)th3[2]->Clone();
		th1[7][j]=(TH1D*)th3[3]->Clone();
		
		th1[5][j]->Add(th3[0]);
		th1[6][j]->Add(th3[0]);
		th1[7][j]->Add(th3[0]);
		
		th1[6][j]->Add(th3[1]);
		th1[7][j]->Add(th3[1]);
		
		th1[7][j]->Add(th3[2]);
		
		th1[4][j]->SetLineColor(5);
		th1[5][j]->SetLineColor(41);
		th1[6][j]->SetLineColor(42);
		th1[7][j]->SetLineColor(43);
		th1[0][j]->SetLineColor(7);
		th1[1][j]->SetLineColor(3);
		th1[2][j]->SetLineColor(4);
		th1[3][j]->SetLineColor(6);
		th1[8][j]->SetLineColor(1);
		
		th1[1][j]->SetLineStyle(1);
		th1[2][j]->SetLineStyle(2);
		th1[3][j]->SetLineStyle(3);
		th1[0][j]->SetLineStyle(4);
		
		
		th1[1][j]->SetLineWidth(th1[0][j]->GetLineWidth()*3);
		th1[2][j]->SetLineWidth(th1[0][j]->GetLineWidth()*3);
		th1[3][j]->SetLineWidth(th1[0][j]->GetLineWidth()*3);
		th1[0][j]->SetLineWidth(th1[0][j]->GetLineWidth()*3);
		th1[8][j]->SetLineWidth(th1[8][j]->GetLineWidth()*3);
		th1[8][j]->SetMarkerSize(th1[8][j]->GetMarkerSize()*2);
		
		th1[4][j]->SetFillColor(5);
		th1[5][j]->SetFillColor(41);
		th1[6][j]->SetFillColor(42);
		th1[7][j]->SetFillColor(43);
		
		if(j!=0)th7[3]->SetLineColor(2);
		if(j!=0)th7[3]->SetLineWidth(th7[3]->GetLineWidth()*2);
		TLegend* leg ;
		leg=new TLegend(0.691452,0.662447,0.890645,0.883966);
		leg->AddEntry(th1[8][j],"Data");
		leg->AddEntry(th1[4][j],"QCD700-1000");																																																												
		leg->AddEntry(th1[5][j],"QCD1000-1500");
		leg->AddEntry(th1[6][j],"QCD1500-2000");
		leg->AddEntry(th1[7][j],"QCD2000-inf");
		if(j!=0)leg->AddEntry(th7[3],"QCDTotalWithSF");
		
		TLegend* leg2 ;																																						
		leg2=new TLegend(0.451452,0.662447,0.670645,0.883966);
		leg2->AddEntry(th1[0][j],"M_{G}=1.2TeV");
		leg2->AddEntry(th1[1][j],"M_{G}=1.8TeV");
		leg2->AddEntry(th1[2][j],"M_{G}=2.5TeV");
		leg2->AddEntry(th1[3][j],"M_{R}=1.8TeV");
		
		if(j==18||j==19)th1[7][j]->GetXaxis()->SetRangeUser(300,1200);
		
		if(j==0)th1[7][j]->SetMaximum(th1[7][j]->GetMaximum()*2);
		else th1[7][j]->SetMaximum(th1[7][j]->GetMaximum()*2);
		
		//c1->SetLogy(1);
		//if(j==20)th1[7][j]->SetMaximum(th1[7][j]->GetMaximum()*50);
		//else th1[7][j]->SetMaximum(th1[7][j]->GetMaximum()*100);
		
		th1[7][j]->SetXTitle(Form("%s",st[j].data()));
		th1[7][j]->Draw("hist");
		th1[6][j]->Draw("hist,ssame");
		th1[5][j]->Draw("hist,ssame");
		th1[4][j]->Draw("hist,ssame");
		th1[0][j]->Draw("hist,ssame");
		th1[1][j]->Draw("hist,ssame");
		th1[2][j]->Draw("hist,ssame");
		th1[3][j]->Draw("hist,ssame");
		
		if(j!=0)th7[3]->Draw("same");
		th1[8][j]->Draw("e,same");
		leg->Draw("same");
		leg2->Draw("same");
		//cout<<"c"<<endl;
		
		c1->Print(Form("plot/%s_afterSF.pdf",st[j].data()));
		
		if(j==0){
			th1[8][j]->Divide(th1[7][j]);
			th1[8][j]->Draw();
			th1[8][j]->SetYTitle("Data/MC");
			c1->Print(Form("plot/%s_ratio_afterSF.pdf",st[j].data()));
		}
		
	}
	
}