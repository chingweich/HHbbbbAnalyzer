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
#include "TF1.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "RooTFnBinding.h" 
#include "RooLandau.h"
#include "RooFFTConvPdf.h"
#include "TVirtualFFT.h"
#include "RooBifurGauss.h"
#include "RooHistPdf.h"
#include "RooBreitWigner.h"
#include "RooNovosibirsk.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "TText.h"
#include "RooChi2Var.h"
#include "RooCBShape.h"
#include "RooVoigtian.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooBifurGauss.h"
#define nMasspoint 11

TCanvas* c1;

using namespace RooFit ;

void makeCorrRoofit2(){
	setNCUStyle(true);
	c1 = new TCanvas("c1","",1360,768);
	
	
	
	TFile* tf1[nMasspoint];
	int masspoint[nMasspoint]={700,800,900,1000,1200,1400,1600,1800,2000,2500,3000};
	for(int i=0;i<nMasspoint;i++){
		tf1[i]=TFile::Open(Form("corr/B%d.root",masspoint[i]));
		//f[i+11]=TFile::Open(Form("R%s.root",masspoint[i].data()));
	}
	double xsec[9]={1.90,0.763,0.33,0.155,7.65e-2,1.58e-2,3.73e-3,2.08e-4,4.91e-5};
		
	double ptBins[14]={200,300,400,500,600,700,800,900,1000,1250,1500,1750,2000,2500};
	double ptBinsCenter[16]={350,450,550,650,750,850,950,1125,1375,1625,1875,2250,2750};
	double ptBinsCenterE[16]={350,450,550,650,750,850,950,1125,1375,1625,1875,2250,2750};
	double ptBinsError[16]={0};
	double ptBinsErrorE[16]={0};
	
	double mean[2][16];
	double sigma[2][16];
	
	RooRealVar xpt("x","x",0,3200) ;
	
	for(int i=0;i<nMasspoint;i++){
		TH1D* th1=(TH1D	*)tf1[i]->Get("ptBarrel");
		ptBinsCenter[i]=th1->GetMean();
		ptBinsError[i]=20;
		th1->SetTitle(Form("%d",masspoint[i]));
		th1->Draw("");
		if(i==0)c1->Print("plots/ptBarrel.pdf(");
		else if(i==nMasspoint-1)c1->Print("plots/ptBarrel.pdf)");
		else  c1->Print("plots/ptBarrel.pdf");
		delete th1;
	}

	for(int i=0;i<nMasspoint;i++){
		TH1D* th1=(TH1D	*)tf1[i]->Get("ptEndcap");
		ptBinsCenterE[i]=th1->GetMean();
		ptBinsErrorE[i]=20;
		th1->SetTitle(Form("%d",masspoint[i]));
		th1->Draw("");
		if(i==0)c1->Print("plots/ptEndcap.pdf(");
		else if(i==nMasspoint-1)c1->Print("plots/ptEndcap.pdf)");
		else  c1->Print("plots/ptEndcap.pdf");
		delete th1;
	}
	
	
	RooRealVar x("x","x",40,150) ;
	
	for(int i=0;i<nMasspoint;i++){
		TH1D* th1=(TH1D	*)tf1[i]->Get("recoBarrelMass");
		th1->SetTitle(Form("%d",masspoint[i]));
		th1->SetMaximum(th1->GetMaximum()*1.3);
		cout<<i<<"="<<Form("%d",masspoint[i])<<endl;

		RooDataHist dh("dh","dh",x,Import(*th1)) ;
		/*RooRealVar mean("mean","mean",125,50,150) ;
		RooRealVar sigmaR("sigmaR","sigmaR",4,0.1,9) ;
		RooRealVar sigmaL("sigmaL","sigmaL",10,0.1,18) ;
		RooBifurGauss Bi("bi","bi",x,mean,sigmaL,sigmaR);
		*/
		
		RooRealVar m0("m0","m0",117,50,150);
		RooRealVar sig("sigma","sigma",12,0,100);
		RooRealVar alpha("alpha","alpha",1,0,10);
		RooRealVar n("n","n",1,0,50);
		RooCBShape CB("CB","CB",x,m0,sig,alpha,n);
		RooPlot* frame=x.frame(Title(Form("%d",masspoint[i])));
		dh.plotOn(frame);
		CB.fitTo(dh) ; 
		CB.plotOn(frame,LineColor(kRed)) ;
		CB.paramOn(frame,Layout(0.55)) ;
		mean[0][i]=125/m0.getValV();
		sigma[0][i]=m0.getError()*(125/m0.getValV())/m0.getValV();
		//mean[4][i]=125/th1->GetMean();
		//sigma[4][i]=th1->GetMeanError()/th1->GetMean();
		//th1->Draw();
		frame->SetMaximum(th1->GetMaximum()*1.3);
		frame->Draw() ;
		
		if(i==0)c1->Print("plots/recoBarrel.pdf(");
		else if(i==nMasspoint-1)c1->Print("plots/recoBarrel.pdf)");
		else  c1->Print("plots/recoBarrel.pdf");
		delete th1;
	}
	
	for(int i=0;i<nMasspoint;i++){
		TH1D* th1=(TH1D	*)tf1[i]->Get("recoEndcapMass");
		th1->SetTitle(Form("%d",masspoint[i]));
		th1->SetMaximum(th1->GetMaximum()*1.3);
		
		RooDataHist dh("dh","dh",x,Import(*th1)) ;
		/*RooRealVar mean("mean","mean",125,50,150) ;
		RooRealVar sigmaR("sigmaR","sigmaR",4,0.1,9) ;
		RooRealVar sigmaL("sigmaL","sigmaL",10,0.1,18) ;
		RooBifurGauss Bi("bi","bi",x,mean,sigmaL,sigmaR);
		*/
		
		RooRealVar m0("m0","m0",117,50,150);
		RooRealVar sig("sigma","sigma",12,0,100);
		RooRealVar alpha("alpha","alpha",1,0,10);
		RooRealVar n("n","n",1,0,50);
		RooCBShape CB("CB","CB",x,m0,sig,alpha,n);
		RooPlot* frame=x.frame(Title(Form("%d",masspoint[i])));
		dh.plotOn(frame);
		CB.fitTo(dh) ; 
		CB.plotOn(frame,LineColor(kRed)) ;
		CB.paramOn(frame,Layout(0.55)) ;
		mean[1][i]=125/m0.getValV();
		sigma[1][i]=m0.getError()*(125/m0.getValV())/m0.getValV();
		//mean[4][i]=125/th1->GetMean();
		//sigma[4][i]=th1->GetMeanError()/th1->GetMean();
		//th1->Draw();
		frame->SetMaximum(th1->GetMaximum()*1.3);
		frame->Draw() ;
		//cout<<i<<"="<<mean[5][i]<<endl;
		if(i==0)c1->Print("plots/recoEndcap.pdf(");
		else if(i==nMasspoint-1)c1->Print("plots/recoEndcap.pdf)");
		else  c1->Print("plots/recoEndcap.pdf");
		delete th1;
	}
	
	
	
	
	TGraphErrors* tg1[6];
	
	tg1[4]=new TGraphErrors(nMasspoint,ptBinsCenter,mean[0],ptBinsError,sigma[0]);
	tg1[5]=new TGraphErrors(nMasspoint,ptBinsCenterE,mean[1],ptBinsErrorE,sigma[1]);
	
	for(int i=0;i<14;i++)cout<<i<<"="<<mean[0][i]<<endl;
	for(int i=0;i<14;i++)cout<<i<<"="<<mean[1][i]<<endl;
	for(int i=0;i<14;i++)cout<<i<<","<<ptBinsCenter[i]<<","<<ptBinsError[i]<<endl;
	
	for(int i=0;i<nMasspoint;i++)cout<<ptBinsCenter[i]<<",";
	cout<<endl;
	for(int i=0;i<nMasspoint;i++)cout<<mean[0][i]<<",";
	cout<<endl;
	for(int i=0;i<nMasspoint;i++)cout<<ptBinsCenterE[i]<<",";
	cout<<endl;
	for(int i=0;i<nMasspoint;i++)cout<<mean[1][i]<<",";
	cout<<endl;
	
	TLegend *leg = new TLegend(0.68, 0.65, 0.94, 0.90);
  
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
	
	
	tg1[4]->GetXaxis()->SetTitle("jet Pt");
	tg1[4]->GetYaxis()->SetTitle("M_{PDG}/M_{Reco}");
	tg1[4]->SetTitle("Gen Correction");
	tg1[4]->SetMinimum(1);
	tg1[4]->SetMaximum(1.4);
	tg1[4]->Draw("APL");
	tg1[4]->SetFillColor(0);
	tg1[5]->SetFillColor(0);
	tg1[5]->SetLineColor(2);
	tg1[5]->SetMarkerColor(2);
	tg1[5]->Draw("PLsame");
	
	tg1[4]->SetTitle("Higgs mass correction");
	/*
	TF1* recoOneBarel = new TF1("genBarel","[0]+[1]*pow(x*[2],-[3])");
	  recoOneBarel->SetParameters(
   				 1.00626,
   				 -1.06161,
   				 0.07999,
   				 1.20454
   				 );
	TF1* recoOneEndcap = new TF1("genEndcap","[0]+[1]*pow(x*[2],-[3])");
	  recoOneEndcap->SetParameters(
   				 1.00626,
   				 -1.06161,
   				 0.07999,
   				 1.20454
   				 );
	*/
	
	
leg->Clear();


 TF1* puppisd_corrGEN      = new TF1("puppisd_corrGEN","[0]+[1]*pow(x*[2],-[3])");
  puppisd_corrGEN->SetParameters(
   				 1.00626,
   				 -1.06161,
   				 0.07999,
   				 1.20454
   				 );
  TF1* puppisd_corrRECO_cen = new TF1("puppisd_corrRECO_cen","([0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5))*([6]+[7]*pow(x*[8],-[9]))",200,2000);
  puppisd_corrRECO_cen->SetParameters(
   				      1.05807,
   				      -5.91971e-05,
   				      2.296e-07,
   				      -1.98795e-10,
   				      6.67382e-14,
   				      -7.80604e-18,
					 1.00626,
   				 -1.06161,
   				 0.07999,
   				 1.20454
   				      );

  TF1* puppisd_corrRECO_for = new TF1("puppisd_corrRECO_for","([0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5))*([6]+[7]*pow(x*[8],-[9]))",200,2000);
  puppisd_corrRECO_for->SetParameters(
   				      1.26638,
   				      -0.000658496,
   				      9.73779e-07,
   				      -5.93843e-10,
   				      1.61619e-13,
   				      -1.6272e-17,
						 1.00626,
   				 -1.06161,
   				 0.07999,
   				 1.20454);
					
					
  
  leg->AddEntry(tg1[4],"H corr barrel");
  leg->AddEntry(tg1[5],"H corr endcap");
  leg->AddEntry(puppisd_corrRECO_cen,"Thea barrel");
  leg->AddEntry(puppisd_corrRECO_for,"Thea endcap");

	leg->Draw("same");
	//genBarel->Draw("same");
//	genEndcap->Draw("same");
	puppisd_corrRECO_cen->SetLineColor(3);
	puppisd_corrRECO_for->SetLineColor(4);
	//tg1[2]->Draw("APL");
	puppisd_corrRECO_cen->Draw("same");
	puppisd_corrRECO_for->Draw("same");
	
	c1->Print("plots/Correction.pdf");
	
}
