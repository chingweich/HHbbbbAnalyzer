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
void massDrawer(){
	setNCUStyle();
	c1 = new TCanvas("c1","",1360,768);
	
	TFile* tf1[2];
	tf1[0]=TFile::Open("output/Zprime1000.root");
	tf1[1]=TFile::Open("output/Zprime2500.root");
	
	TH1D * th1[8];
	th1[0]= (TH1D *)tf1[0]->FindObjectAny("FATjetPRmass");
	th1[1]= (TH1D *)tf1[0]->FindObjectAny("FATjetPRmassL2L3Corr");
	th1[2]= (TH1D *)tf1[0]->FindObjectAny("FATjetPuppiSDmass");
	th1[3]= (TH1D *)tf1[0]->FindObjectAny("FATjetPuppiSDmassL2L3Corr");
	th1[4]= (TH1D *)tf1[1]->FindObjectAny("FATjetPRmass");
	th1[5]= (TH1D *)tf1[1]->FindObjectAny("FATjetPRmassL2L3Corr");
	th1[6]= (TH1D *)tf1[1]->FindObjectAny("FATjetPuppiSDmass");
	th1[7]= (TH1D *)tf1[1]->FindObjectAny("FATjetPuppiSDmassL2L3Corr");
	
	for(int i=0;i<8;i++)th1[i]->Rebin(5);
	
	for(int i=0;i<8;i++)th1[i]->Scale(1/th1[i]->Integral());
	
	th1[0]->SetMaximum(th1[0]->GetMaximum()*1.3);
	th1[1]->SetMaximum(th1[1]->GetMaximum()*1.3);
	th1[0]->SetLineColor(1);
	th1[2]->SetLineColor(2);
	th1[4]->SetLineColor(3);
	th1[6]->SetLineColor(4);
	
	th1[1]->SetLineColor(1);
	th1[3]->SetLineColor(2);
	th1[5]->SetLineColor(3);
	th1[7]->SetLineColor(4);
	
	TLegend* leg ;
	leg=new TLegend(0.141452,0.552447,0.450645,0.863966);
	leg->AddEntry(th1[0],"1TeV,PRmass");
	leg->AddEntry(th1[2],"1TeV,PuppiSDmass");
	leg->AddEntry(th1[4],"2.5TeV,PRmass");
	leg->AddEntry(th1[6],"2.5TeV,PuppiSDmass");
	
	th1[0]->SetTitle("without L2L3Corr");
	th1[1]->SetTitle("with L2L3Corr");
	th1[0]->SetXTitle("Mass[GeV]");
	th1[1]->SetXTitle("Mass[GeV]");
	th1[0]->Draw("c");
	th1[2]->Draw("same,c");
	th1[4]->Draw("same,c");
	th1[6]->Draw("same,c");
	leg->Draw("same");
	c1->Print("output/wo.pdf");
	
	leg->Clear();
	leg->AddEntry(th1[1],"1TeV,PRmassL2L3Corr");
	leg->AddEntry(th1[3],"1TeV,PuppiSDmassL2L3Corr");
	leg->AddEntry(th1[5],"2.5TeV,PRmassL2L3Corr");
	leg->AddEntry(th1[7],"2.5TeV,PuppiSDmassL2L3Corr");
	
	
	th1[1]->Draw("c");
	th1[3]->Draw("same,c");
	th1[5]->Draw("same,c");
	th1[7]->Draw("same,c");
	leg->Draw("same");
	c1->Print("output/with.pdf");

}

