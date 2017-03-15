#include <cstring>
#include <cerrno>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <unistd.h>
#include <errno.h>
#include <iomanip>
// ROOT headers
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TIterator.h"
#include "TMath.h"

#include "TLatex.h"
#include "TString.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"

// RooFit headers
#include "RooAbsPdf.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooHistPdf.h"
#include "RooMsgService.h"
#include "RooNLLVar.h"
#include "RooPlot.h"
#include "RooRandom.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "TStyle.h"

// RooStats headers
#include "RooStats/HLFactory.h"

#include "RooAbsPdf.h"
#include "RooHist.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAbsData.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooProduct.h"
#include "RooExtendPdf.h"
#include "RooBernstein.h"
#include "RooFitResult.h"
#include "RooMinimizer.h"
#include "RooCmdArg.h"
#include "RooConstVar.h"
#include "RooRealVar.h"

using namespace RooFit;
using namespace RooStats ;
using namespace std ;

double WSMakerForRSS(int j=0){
	
	string DBT[2]={"TT","LL"};
	string SB="_sb";
	bool isSB=1;
	double rss0[2]={0};
	//for(int j=0;j<2;j++){
		TFile * tf1 = TFile::Open(Form("w_data_%s.root",DBT[j].data()));
		RooWorkspace *ws1 = (RooWorkspace*)tf1->Get("HH4b");
		//ws1->Print();
		RooDataSet* data=(RooDataSet*)ws1->data("data_obs");
		if(isSB)data=(RooDataSet*)ws1->data("data_obs_sb");
		TFile * tf2 = TFile::Open(Form("w_background_%s.root",DBT[j].data()));
		RooWorkspace *ws2 = (RooWorkspace*)tf2->Get("HH4b");
		//ws2->Print();
		RooAbsPdf* bkg=ws2->pdf("bg_");
		if(isSB)bkg=ws2->pdf("bgSB_");
		RooRealVar *mass =  ws2->var("x");
		
		RooPlot *plot = mass->frame();
		
		data->plotOn(plot,RooFit::LineColor(kBlack));
		bkg->plotOn(plot,RooFit::LineColor(kGreen));
		plot->SetTitle("PDF fits to toy data");
		
		plot->Draw();
		
		RooArgSet* set = new RooArgSet(*mass);
		
		
		double normalisation =  data->sumEntries();
		cout<<normalisation<<endl;
		for (int i = 0; i < (plot->getHist(Form("h_data_obs%s",SB.data())))->GetN()  ; i++){
			double x0, y0;
			(plot->getHist(Form("h_data_obs%s",SB.data())))->GetPoint(i, x0, y0);
			double xc, yc;
			double xlowErr =(plot->getHist(Form("h_data_obs%s",SB.data())))->GetErrorXlow(i);
			double xhighErr = (plot->getHist(Form("h_data_obs%s",SB.data())))->GetErrorXhigh(i);
			(plot->getHist(Form("h_data_obs%s",SB.data())))->GetPoint(i, xc, yc);
			double xmin = xc-fabs(xlowErr), xmax = xc+fabs(xhighErr);
			mass->setRange("A",xmin,xmax);
			RooAbsReal* intBin0 =bkg->createIntegral(*set,*set,"A") ;
			double dintBin0 = intBin0->getVal();
			rss0[j]+=TMath::Power((dintBin0*normalisation - y0),2);
			if(dintBin0*normalisation - y0>0.1)cout<<dintBin0<<","<<normalisation<<","<<dintBin0*normalisation <<","<<y0<<","<<dintBin0*normalisation - y0<<endl;
		}
	//}
	
	
	
	//for(int j=0;j<2;j++)cout<<DBT[j]<<",RSS="<<rss0[j]<<endl;
	//cout<<"RSS="<<rss0[1]<<endl;
	return rss0[j];
}

double WSMakerForRSS3(int j,int & dataN){
	
	string DBT[2]={"TT","LL"};
	string SB="_sb";
	bool isSB=1;
	double rss0[2]={0};
	//for(int j=0;j<2;j++){
		TFile * tf1 = TFile::Open(Form("w_data_%s.root",DBT[j].data()));
		RooWorkspace *ws1 = (RooWorkspace*)tf1->Get("HH4b");
		//ws1->Print();
		RooDataSet* data=(RooDataSet*)ws1->data("data_obs");
		if(isSB)data=(RooDataSet*)ws1->data("data_obs_sb");
		TFile * tf2 = TFile::Open(Form("w_background_%s.root",DBT[j].data()));
		RooWorkspace *ws2 = (RooWorkspace*)tf2->Get("HH4b");
		//ws2->Print();
		RooAbsPdf* bkg=ws2->pdf("bg_");
		if(isSB)bkg=ws2->pdf("bgSB1c_");
		RooRealVar *mass =  ws2->var("x");
		
		RooPlot *plot = mass->frame();
		
		data->plotOn(plot,RooFit::LineColor(kBlack));
		bkg->plotOn(plot,RooFit::LineColor(kGreen));
		plot->SetTitle("PDF fits to toy data");
		
	
		plot->Draw();
	
		RooArgSet* set = new RooArgSet(*mass);
		
		
		double normalisation =  data->sumEntries();
		cout<<normalisation<<endl;
		for (int i = 0; i < (plot->getHist(Form("h_data_obs%s",SB.data())))->GetN()  ; i++){
			double x0, y0;
			(plot->getHist(Form("h_data_obs%s",SB.data())))->GetPoint(i, x0, y0);
			double xc, yc;
			double xlowErr =(plot->getHist(Form("h_data_obs%s",SB.data())))->GetErrorXlow(i);
			double xhighErr = (plot->getHist(Form("h_data_obs%s",SB.data())))->GetErrorXhigh(i);
			(plot->getHist(Form("h_data_obs%s",SB.data())))->GetPoint(i, xc, yc);
			double xmin = xc-fabs(xlowErr), xmax = xc+fabs(xhighErr);
			mass->setRange("A",xmin,xmax);
			RooAbsReal* intBin0 =bkg->createIntegral(*set,*set,"A") ;
			double dintBin0 = intBin0->getVal();
			rss0[j]+=TMath::Power((dintBin0*normalisation - y0),2);
			if(dintBin0*normalisation - y0>0.1)cout<<dintBin0<<","<<normalisation<<","<<dintBin0*normalisation <<","<<y0<<","<<dintBin0*normalisation - y0<<endl;
		}
	//}
	
	
	dataN= (plot->getHist(Form("h_data_obs%s",SB.data())))->GetN()  ;
	//for(int j=0;j<2;j++)cout<<DBT[j]<<",RSS="<<rss0[j]<<endl;
	//cout<<"RSS="<<rss0[1]<<endl;
	return rss0[j];
}

double WSMakerForRSS1(int j,int & dataN){
	
	string DBT[2]={"TT","LL"};
	string SB="_sb";
	bool isSB=1;
	double rss0[2]={0};
	//for(int j=0;j<2;j++){
		TFile * tf1 = TFile::Open(Form("w_data_%s.root",DBT[j].data()));
		RooWorkspace *ws1 = (RooWorkspace*)tf1->Get("HH4b");
		//ws1->Print();
		RooDataSet* data=(RooDataSet*)ws1->data("data_obs");
		if(isSB)data=(RooDataSet*)ws1->data("data_obs_sb");
		TFile * tf2 = TFile::Open(Form("w_background_%s.root",DBT[j].data()));
		RooWorkspace *ws2 = (RooWorkspace*)tf2->Get("HH4b");
		//ws2->Print();
		RooAbsPdf* bkg=ws2->pdf("bg_");
		if(isSB)bkg=ws2->pdf("bgSB1a_");
		RooRealVar *mass =  ws2->var("x");
		
		RooPlot *plot = mass->frame();
		
		data->plotOn(plot,RooFit::LineColor(kBlack));
		bkg->plotOn(plot,RooFit::LineColor(kGreen));
		plot->SetTitle("PDF fits to toy data");
		
		
		plot->Draw();
		
		RooArgSet* set = new RooArgSet(*mass);
		
		
		double normalisation =  data->sumEntries();
		cout<<normalisation<<endl;
		for (int i = 0; i < (plot->getHist(Form("h_data_obs%s",SB.data())))->GetN()  ; i++){
			double x0, y0;
			(plot->getHist(Form("h_data_obs%s",SB.data())))->GetPoint(i, x0, y0);
			double xc, yc;
			double xlowErr =(plot->getHist(Form("h_data_obs%s",SB.data())))->GetErrorXlow(i);
			double xhighErr = (plot->getHist(Form("h_data_obs%s",SB.data())))->GetErrorXhigh(i);
			(plot->getHist(Form("h_data_obs%s",SB.data())))->GetPoint(i, xc, yc);
			double xmin = xc-fabs(xlowErr), xmax = xc+fabs(xhighErr);
			mass->setRange("A",xmin,xmax);
			RooAbsReal* intBin0 =bkg->createIntegral(*set,*set,"A") ;
			double dintBin0 = intBin0->getVal();
			rss0[j]+=TMath::Power((dintBin0*normalisation - y0),2);
			if(dintBin0*normalisation - y0>0.1)cout<<dintBin0<<","<<normalisation<<","<<dintBin0*normalisation <<","<<y0<<","<<dintBin0*normalisation - y0<<endl;
		}
	//}
	
	
	dataN= (plot->getHist(Form("h_data_obs%s",SB.data())))->GetN()  ;
	//for(int j=0;j<2;j++)cout<<DBT[j]<<",RSS="<<rss0[j]<<endl;
	//cout<<"RSS="<<rss0[1]<<endl;
	return rss0[j];
}

void ftest(){
	double rss0[2],rss3[2],rss1[2];
	int dataN=0;
	rss0[0]=WSMakerForRSS();
	rss3[0]=WSMakerForRSS3(0,dataN);
	rss1[0]=WSMakerForRSS1(0,dataN);
	
	double fTT3=(rss0[0]-rss3[0])/(rss3[0]/(dataN-3));
	double fTT1=(rss1[0]-rss0[0])/(rss0[0]/(dataN-2));
	//cout<<"fTT="<<fTT3<<endl;
	
	rss0[1]=WSMakerForRSS(1);
	rss3[1]=WSMakerForRSS3(1,dataN);
	rss1[1]=WSMakerForRSS1(1,dataN);
	//cout<<rss0[1]<<","<<rss3[1]<<endl;
	double fLL3=(rss0[1]-rss3[1])/(rss3[1]/(dataN-3));
	double fLL1=(rss1[1]-rss0[1])/(rss0[1]/(dataN-2));
	
	cout<<"rss 1par="<<rss1[0]<<" , rss 2par="<<rss0[0]<<" , rss 3par="<<rss3[0]<<endl;
	cout<<"fTT2to3="<<fTT3<<",1to2="<<fTT1<<endl;
	cout<<"rss 1par="<<rss1[1]<<" , rss 2par="<<rss0[1]<<" , rss 3par="<<rss3[1]<<endl;
	cout<<"fLL2to3="<<fLL3<<",1to2="<<fLL1<<endl;
	cout<<"P(TT)3="<<1-TMath::FDistI(fTT3,1,40)<<endl;
	cout<<"P(TT)2="<<1-TMath::FDistI(fTT1,1,41)<<endl;
	cout<<"P(LL)3="<<1-TMath::FDistI(fLL3,1,40)<<endl;
	cout<<"P(LL)2="<<1-TMath::FDistI(fLL1,1,41)<<endl;
	//cout<<"fLL="<<fLL3<<endl;
}