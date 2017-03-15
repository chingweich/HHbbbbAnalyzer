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

void WSMakerForRSS(){
	
	string DBT[2]={"TT","LL"};
	string SB="_sb";
	bool isSB=1;
	double rss0[2]={0};
	for(int j=0;j<2;j++){
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
		
		TCanvas* c1;
		c1 = new TCanvas("c1","", 600, 600);
		plot->Draw();
		if(isSB) c1->Print(Form("%s_SB.pdf",DBT[j].data()));
		else c1->Print(Form("%s.pdf",DBT[j].data()));
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
			cout<<dintBin0<<","<<normalisation<<","<<dintBin0*normalisation <<","<<y0<<","<<dintBin0*normalisation - y0<<endl;
		}
	}
	
	
	
	for(int j=0;j<2;j++)cout<<DBT[j]<<",RSS="<<rss0[j]<<endl;
	//cout<<"RSS="<<rss0[1]<<endl;
}