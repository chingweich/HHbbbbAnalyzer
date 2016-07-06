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
#include"TGraphAsymmErrors.h"

#define nMassp 10
#define nSt 11

TCanvas* c1;
void massDrawerBulk(){
	setNCUStyle();
	c1 = new TCanvas("c1","",1360,768);
	
	TFile* tf1[2];
	tf1[0]=TFile::Open("output/Bulk1200.root");
	tf1[1]=TFile::Open("output/Bulk4000.root");
	
	TH1D * th1[16];
	th1[0]= (TH1D *)tf1[0]->FindObjectAny("FATjetPRmass");
	th1[1]= (TH1D *)tf1[0]->FindObjectAny("FATjetPRmassL2L3Corr");
	th1[2]= (TH1D *)tf1[0]->FindObjectAny("FATjetPuppiSDmass");
	th1[3]= (TH1D *)tf1[0]->FindObjectAny("FATjetPuppiSDmassL2L3Corr");
	th1[4]= (TH1D *)tf1[1]->FindObjectAny("FATjetPRmass");
	th1[5]= (TH1D *)tf1[1]->FindObjectAny("FATjetPRmassL2L3Corr");
	th1[6]= (TH1D *)tf1[1]->FindObjectAny("FATjetPuppiSDmass");
	th1[7]= (TH1D *)tf1[1]->FindObjectAny("FATjetPuppiSDmassL2L3Corr");
	th1[8]= (TH1D *)tf1[0]->FindObjectAny("FATjetPRmassS");
	th1[9]= (TH1D *)tf1[0]->FindObjectAny("FATjetPRmassL2L3CorrS");
	th1[10]= (TH1D *)tf1[0]->FindObjectAny("FATjetPuppiSDmassS");
	th1[11]= (TH1D *)tf1[0]->FindObjectAny("FATjetPuppiSDmassL2L3CorrS");
	th1[12]= (TH1D *)tf1[1]->FindObjectAny("FATjetPRmassS");
	th1[13]= (TH1D *)tf1[1]->FindObjectAny("FATjetPRmassL2L3CorrS");
	th1[14]= (TH1D *)tf1[1]->FindObjectAny("FATjetPuppiSDmassS");
	th1[15]= (TH1D *)tf1[1]->FindObjectAny("FATjetPuppiSDmassL2L3CorrS");
	
	for(int i=0;i<16;i++)th1[i]->Rebin(2);
	
	cout<<th1[0]->GetBinCenter(th1[0]->GetMaximumBin())<<","<<th1[4]->GetBinCenter(th1[4]->GetMaximumBin())<<endl;
	cout<<th1[2]->GetBinCenter(th1[2]->GetMaximumBin())<<","<<th1[6]->GetBinCenter(th1[6]->GetMaximumBin())<<endl;
	cout<<th1[1]->GetBinCenter(th1[1]->GetMaximumBin())<<","<<th1[5]->GetBinCenter(th1[5]->GetMaximumBin())<<endl;
	cout<<th1[3]->GetBinCenter(th1[3]->GetMaximumBin())<<","<<th1[7]->GetBinCenter(th1[7]->GetMaximumBin())<<endl;
	
	//for(int i=0;i<16;i++)th1[i]->Rebin(2.5);
	
	
	
	for(int i=0;i<16;i++)th1[i]->Scale(1/th1[i]->Integral());
	
	th1[0]->SetMaximum(th1[0]->GetMaximum()*1.3);
	th1[1]->SetMaximum(th1[1]->GetMaximum()*1.3);
	th1[8]->SetMaximum(th1[8]->GetMaximum()*1.3);
	th1[9]->SetMaximum(th1[9]->GetMaximum()*1.3);
	
	th1[0]->SetLineColor(1);
	th1[2]->SetLineColor(2);
	th1[4]->SetLineColor(3);
	th1[6]->SetLineColor(4);
	
	th1[1]->SetLineColor(1);
	th1[3]->SetLineColor(2);
	th1[5]->SetLineColor(3);
	th1[7]->SetLineColor(4);
	
	th1[8]->SetLineColor(1);
	th1[10]->SetLineColor(2);
	th1[12]->SetLineColor(3);
	th1[14]->SetLineColor(4);
	
	th1[9]->SetLineColor(1);
	th1[11]->SetLineColor(2);
	th1[13]->SetLineColor(3);
	th1[15]->SetLineColor(4);
	
	TLegend* leg ;
	leg=new TLegend(0.141452,0.552447,0.450645,0.863966);
	leg->AddEntry(th1[0],"1.2TeV,PRmass");
	leg->AddEntry(th1[2],"1.2TeV,PuppiSDmass");
	leg->AddEntry(th1[4],"4TeV,PRmass");
	leg->AddEntry(th1[6],"4TeV,PuppiSDmass");
	
	
	th1[0]->SetTitle("without L2L3Corr");
	th1[1]->SetTitle("with L2L3Corr");
	th1[0]->SetXTitle("Mass[GeV]");
	th1[1]->SetXTitle("Mass[GeV]");
	th1[0]->Draw("c");
	th1[2]->Draw("same,c");
	th1[4]->Draw("same,c");
	th1[6]->Draw("same,c");
	leg->Draw("same");
	c1->Print("output/bulk_leading_wo.pdf");
	
	th1[8]->Draw("");
	th1[10]->Draw("same");
	th1[12]->Draw("same");
	th1[14]->Draw("same");
	leg->Draw("same");
	c1->Print("output/bulk_subLeading_wo.pdf");
	
	th1[0]->Add(th1[8]);
	th1[2]->Add(th1[10]);
	th1[4]->Add(th1[12]);
	th1[6]->Add(th1[14]);
	th1[0]->SetMaximum(th1[0]->GetMaximum()*1.3);
	th1[0]->Draw("c");
	th1[2]->Draw("same,c");
	th1[4]->Draw("same,c");
	th1[6]->Draw("same,c");
	leg->Draw("same");
	c1->Print("output/bulk_L+S_wo.pdf");
	
	leg->Clear();
	leg->AddEntry(th1[1],"1.2TeV,PRmassL2L3Corr");
	leg->AddEntry(th1[3],"1.2TeV,PuppiSDmassL2L3Corr");
	leg->AddEntry(th1[5],"4TeV,PRmassL2L3Corr");
	leg->AddEntry(th1[7],"4TeV,PuppiSDmassL2L3Corr");
	
	th1[1]->Draw("c");
	th1[3]->Draw("same,c");
	th1[5]->Draw("same,c");
	th1[7]->Draw("same,c");
	leg->Draw("same");
	c1->Print("output/bulk_leading_with.pdf");
	
	th1[9]->Draw("");
	th1[11]->Draw("same");
	th1[13]->Draw("same");
	th1[15]->Draw("same");
	leg->Draw("same");
	c1->Print("output/bulk_subLeading_with.pdf");
	
	
	
	
	th1[1]->Add(th1[9]);
	th1[3]->Add(th1[11]);
	th1[5]->Add(th1[13]);
	th1[7]->Add(th1[15]);
	th1[1]->SetMaximum(th1[1]->GetMaximum()*1.3);
	th1[1]->Draw("c");
	th1[3]->Draw("same,c");
	th1[5]->Draw("same,c");
	th1[7]->Draw("same,c");
	leg->Draw("same");
	c1->Print("output/bulk_L+S_with.pdf");
	
	cout<<th1[0]->GetBinCenter(th1[0]->GetMaximumBin())<<","<<th1[4]->GetBinCenter(th1[4]->GetMaximumBin())<<endl;
	cout<<th1[2]->GetBinCenter(th1[2]->GetMaximumBin())<<","<<th1[6]->GetBinCenter(th1[6]->GetMaximumBin())<<endl;
	cout<<th1[1]->GetBinCenter(th1[1]->GetMaximumBin())<<","<<th1[5]->GetBinCenter(th1[5]->GetMaximumBin())<<endl;
	cout<<th1[3]->GetBinCenter(th1[3]->GetMaximumBin())<<","<<th1[7]->GetBinCenter(th1[7]->GetMaximumBin())<<endl;
	
	TH1D * th2[20];
	th2[0]= (TH1D *)tf1[0]->FindObjectAny("PV");
	th2[1]= (TH1D *)tf1[0]->FindObjectAny("PVFATjetPRmass");
	th2[2]= (TH1D *)tf1[0]->FindObjectAny("PVFATjetPRmassL2L3Corr");
	th2[3]= (TH1D *)tf1[0]->FindObjectAny("PVFATjetPuppiSDmass");
	th2[4]= (TH1D *)tf1[0]->FindObjectAny("PVFATjetPuppiSDmassL2L3Corr");
	
	th2[5]= (TH1D *)tf1[1]->FindObjectAny("PV");
	th2[6]= (TH1D *)tf1[1]->FindObjectAny("PVFATjetPRmass");
	th2[7]= (TH1D *)tf1[1]->FindObjectAny("PVFATjetPRmassL2L3Corr");
	th2[8]= (TH1D *)tf1[1]->FindObjectAny("PVFATjetPuppiSDmass");
	th2[9]= (TH1D *)tf1[1]->FindObjectAny("PVFATjetPuppiSDmassL2L3Corr");
	
	th2[10]= (TH1D *)tf1[0]->FindObjectAny("PVFATjetPRmassTau21");
	th2[11]= (TH1D *)tf1[0]->FindObjectAny("PVFATjetPuppiSDmassPuppiTau21");
	th2[12]= (TH1D *)tf1[0]->FindObjectAny("PVFATjetPRmassL2L3CorrTau21");
	th2[13]= (TH1D *)tf1[0]->FindObjectAny("PVFATjetPuppiSDmassL2L3CorrPuppiTau21");
	
	
	TGraphAsymmErrors * tg1[20];
	tg1[0]= new TGraphAsymmErrors(th2[1],th2[0]);
	tg1[1]= new TGraphAsymmErrors(th2[2],th2[0]);
	tg1[2]= new TGraphAsymmErrors(th2[3],th2[0]);
	tg1[3]= new TGraphAsymmErrors(th2[4],th2[0]);
	tg1[4]= new TGraphAsymmErrors(th2[6],th2[5]);
	tg1[5]= new TGraphAsymmErrors(th2[7],th2[5]);
	tg1[6]= new TGraphAsymmErrors(th2[8],th2[5]);
	tg1[7]= new TGraphAsymmErrors(th2[9],th2[5]);
	tg1[8]= new TGraphAsymmErrors(th2[10],th2[0]);
	tg1[9]= new TGraphAsymmErrors(th2[11],th2[0]);
	tg1[10]= new TGraphAsymmErrors(th2[12],th2[0]);
	tg1[11]= new TGraphAsymmErrors(th2[13],th2[0]);
	
	tg1[0]->SetLineColor(1);
	tg1[2]->SetLineColor(2);
	tg1[4]->SetLineColor(3);
	tg1[6]->SetLineColor(4);
	
	tg1[1]->SetLineColor(1);
	tg1[3]->SetLineColor(2);
	tg1[5]->SetLineColor(3);
	tg1[7]->SetLineColor(4);
	
	tg1[8]->SetLineColor(3);
	tg1[9]->SetLineColor(4);
	tg1[10]->SetLineColor(3);
	tg1[11]->SetLineColor(4);
	
	TLegend* leg2 ;
	leg2=new TLegend(0.241452,0.16447,0.550645,0.43966);
	
	leg2->AddEntry(tg1[0],"1.2TeV,PRmass");
	leg2->AddEntry(tg1[2],"1.2TeV,PuppiSDmass");
	leg2->AddEntry(tg1[4],"4TeV,PRmass");
	leg2->AddEntry(tg1[6],"4TeV,PuppiSDmass");
	
	tg1[0]->Draw("APL");
	tg1[2]->Draw("PL,same");
	tg1[4]->Draw("PL,same");
	tg1[6]->Draw("PL,same");
	leg2->Draw("same");
	c1->Print("output/bulk_PV_wo.pdf");
	
	leg2->Clear();
	leg2->AddEntry(tg1[1],"1.2TeV,PRmassL2L3Corr");
	leg2->AddEntry(tg1[3],"1.2TeV,PuppiSDmassL2L3Corr");
	leg2->AddEntry(tg1[5],"4TeV,PRmassL2L3Corr");
	leg2->AddEntry(tg1[7],"4TeV,PuppiSDmassL2L3Corr");
	
	
	tg1[0]->SetTitle("without L2L3Corr");
	tg1[1]->SetTitle("with L2L3Corr");
	tg1[1]->GetXaxis()->SetTitle("PV");
	tg1[0]->GetXaxis()->SetTitle("PV");
	tg1[1]->GetYaxis()->SetTitle("efficiency");
	tg1[0]->GetYaxis()->SetTitle("efficiency");
	
	tg1[1]->Draw("APL");
	tg1[3]->Draw("PL,same");
	tg1[5]->Draw("PL,same");
	tg1[7]->Draw("PL,same");
	leg2->Draw("same");
	c1->Print("output/bulk_PV_with.pdf");
	
	
	TLegend* leg3 ;
	leg3=new TLegend(0.141452,0.66447,0.550645,0.93966);
	leg3->Clear();
	leg3->AddEntry(tg1[0],"1.2TeV,PRmass");
	leg3->AddEntry(tg1[2],"1.2TeV,PuppiSDmass");
	leg3->AddEntry(tg1[8],"1.2TeV,PRmass+tau21");
	leg3->AddEntry(tg1[9],"1.2TeV,PuppiSDmass+puppitau21");
	
	tg1[0]->SetMaximum(1.5);
	tg1[1]->SetMaximum(1.5);
	
	tg1[0]->Draw("APL");
	tg1[2]->Draw("PL,same");
	tg1[8]->Draw("PL,same");
	tg1[9]->Draw("PL,same");
	leg3->Draw("same");
	c1->Print("output/bulk1200_PV_tau21_wo.pdf");
	
	leg3->Clear();
	leg3->AddEntry(tg1[1],"1.2TeV,PRmassL2L3Corr");
	leg3->AddEntry(tg1[3],"1.2TeV,PuppiSDmassL2L3Corr");
	leg3->AddEntry(tg1[10],"1.2TeV,PRmassL2L3Corr+tau21");
	leg3->AddEntry(tg1[11],"1.2TeV,PuppiSDmassL2L3Corr+puppitau21");
	
	tg1[1]->Draw("APL");
	tg1[3]->Draw("PL,same");
	tg1[10]->Draw("PL,same");
	tg1[11]->Draw("PL,same");
	leg3->Draw("same");
	c1->Print("output/bulk1200_PV_tau21_with.pdf");
	
	TH1D * th3[16];
	th3[0]= (TH1D *)tf1[0]->FindObjectAny("tau21");
	th3[1]= (TH1D *)tf1[0]->FindObjectAny("tau21Puppi");
	th3[0]->SetLineColor(1);
	th3[1]->SetLineColor(2);
	
	
	TLegend* leg4 ;
	leg4=new TLegend(0.641452,0.76447,0.850645,0.87966);
	leg4->Clear();
	leg4->AddEntry(th3[0],"tau21");
	leg4->AddEntry(th3[1],"tau21Puppi");
	th3[0]->Draw("");
	th3[1]->Draw("same");
	leg4->Draw("same");
	c1->Print("output/tau21/HH.pdfj");
}

