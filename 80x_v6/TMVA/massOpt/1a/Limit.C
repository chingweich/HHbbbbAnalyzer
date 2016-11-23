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
#include"../../../../setNCUStyle.C"

TCanvas* c1;

#define  nWidth 4
#define  nBmin 4

void Limit(){
	setNCUStyle(true);
	c1 = new TCanvas("c1","",800,600);
	
	int width [nWidth]={25,30,35,40};
	int bmin[nBmin]={100,105,110,115};
	string  masspoint[11]={"1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500"};
	
	TLegend *leg = new TLegend(0.75, 0.68, 0.96, 0.95);
  
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
	
	for(int k=1;k<nWidth;k++){
		 for(int m=0;m<nBmin;m++){
			 
			 if(width[k]+bmin[m]>166)continue;
			 
			  if( 
			 // //(k==4 && m==6)||
			  (k==4 && m==5)||
			     (k==4 && m==4)||
			   (k==4 && m==3)||
			   // (k==4 && m==2)||
			   // (k==4 && m==1)||
			    (k==4 && m==0)||
			 
			 // // // // (k==3 && m==7)||
			 // // // // // (k==3 && m==6)||
			 // (k==3 && m==5)||
			    // // (k==3 && m==4)||
			   // (k==3 && m==3)||
			   // (k==3 && m==2)||
			   // (k==3 && m==1)||
			     // (k==3 && m==0)||
			 
			      (k==2 && m==0)||
			      (k==2 && m==1)||
			// //(k==2 && m==2)||
			   // (k==2 && m==3)||
			    // // (k==2 && m==4)||
		 // (k==2 && m==5)||
			 // // // // // (k==2 && m==6)||
			 
			     (k==1 && m==0)||
			     // (k==1 && m==1)||
			    (k==1 && m==2)||
			// (k==1 && m==3)||
			     // (k==1 && m==4)||
			    // (k==1 && m==5)||
			   
			     // (k==0 && m==0)||
			     // (k==0 && m==1)||
			    // (k==0 && m==2)||
			     // (k==0 && m==3)||
			    // // (k==0 && m==4)||
			 
			  // (k==0 && m==6)
			  (k==10 && m==6)
				 )continue;
			 
			 TFile* tf1=TFile::Open(Form("MassPlotFineBins_subtr_Moriond_Silver%dto%d.root",bmin[m],bmin[m]+width[k]));
			 if (!tf1 || !tf1->IsOpen())continue;
			 TGraphAsymmErrors* tg1=(TGraphAsymmErrors*)tf1->Get("LimitExpectedCLs");
			 
			 tg1->GetYaxis()->SetTitle("95% CLs on #sigma(X#rightarrowHH)#timesBR(HH#rightarrowb#bar{b}b#bar{b})[fb]");
			 tg1->SetLineStyle(k+1);
			 tg1->SetFillColor(0);
			 tg1->SetLineColor(m+1);
			 leg->AddEntry(tg1,Form("%dto%d",bmin[m],bmin[m]+width[k]));
			 tg1->SetMaximum(12);
			 tg1->SetMinimum(1);
			 cout<<k<<","<<m<<endl;
			 //c1->SetLogy(1);
			 if(k==1 && m==1)tg1->Draw("APL");
			 else tg1->Draw("PL same");
		}
	}
	leg->Draw("same");
	c1->Print("Limit.pdf");
}
