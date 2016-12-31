#include <TLegend.h>
#include <vector>
#include <iostream>
#include <fstream>
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
#include <TObject.h>
#include "TImage.h"
#include "TSystem.h"
#include "TStyle.h"
#include "../../untuplizer.h"
#include <TClonesArray.h>
#include <fstream>
#include <cmath>
#include <TSystem.h>
#include <string>
#include <sstream>
#include "../../setNCUStyle.C"


TH1D* massPlotBase(string input,string variable,double down=100,double up=1000){
	TH1D* th1=new TH1D("","",100,down,up);
	TFile* tf1;
	tf1=TFile::Open(Form("%s",input.data()));
	TTree* tree;
	tf1->GetObject("mynewTree",tree);
	TreeReader data(tree);
				//data.Print();
	for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
			data.GetEntry(jEntry);
			Float_t  vari = data.GetFloat(Form("%s",variable.data()));
			//Float_t  jet2_puppi_msoftdrop = data.GetFloat("jet2_puppi_msoftdrop");
			th1->Fill(vari);
	}
	return th1;
}

void variableComparison(){
	TH1D* th1[3];
	//th1[0]=massPlotBase("BulkGrav_M-2000_0.root","jet1_puppi_msoftdrop");
	//th1[1]=massPlotBase("BulkGrav_M-2000_0.root","jet1_puppi_msoftdrop_TheaCorr");
	//th1[2]=massPlotBase("HCorr.root","jet1_puppi_msoftdrop_TheaCorr");
	
	TCanvas* c1;
	setNCUStyle(true);
	c1 = new TCanvas("c1","",800,600);
	
	TLegend *leg = new TLegend(0.70, 0.63, 0.85, 0.88);
  
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);
	leg->SetTextSize(0.04);
  
	vector<string> st;
	vector<double> up;
	vector<double> down;
	st.push_back("jet1pt"); up.push_back(800);  down.push_back(100);  
	st.push_back("jet2pt"); up.push_back(800);  down.push_back(100);  
	st.push_back("jet1eta"); up.push_back(3);  down.push_back(-3);  
	st.push_back("jet2eta"); up.push_back(3);  down.push_back(-3);  
	st.push_back("etadiff"); up.push_back(1.5);  down.push_back(0);  
	st.push_back("dijetmass_softdrop_corr"); up.push_back(1500);  down.push_back(700);  
	st.push_back("jet1bbtag"); up.push_back(1);  down.push_back(-1);  
	st.push_back("jet2bbtag"); up.push_back(1);  down.push_back(-1);  
	st.push_back("jet1_puppi_tau21"); up.push_back(1);  down.push_back(0);  
	st.push_back("jet2_puppi_tau21"); up.push_back(1);  down.push_back(0);  
	st.push_back("jet1_puppi_msoftdrop"); up.push_back(150);  down.push_back(40);  
	st.push_back("jet2_puppi_msoftdrop"); up.push_back(150);  down.push_back(40);  
	st.push_back("jet1_puppi_msoftdrop_TheaCorr"); up.push_back(150);  down.push_back(40);  
	st.push_back("jet2_puppi_msoftdrop_TheaCorr"); up.push_back(150);  down.push_back(40);  
	st.push_back("");
	
	string legendName[2]={"prompt","reReco"};
	
	for(unsigned int j=0;j<st.size()-1;j++){
		th1[0]=massPlotBase("JetHT.root",st[j],down[j],up[j]);
		th1[1]=massPlotBase("reRecoMinitree.root",st[j],down[j],up[j]);
		
		leg->Clear();
		for(int i=0;i<2;i++){
			th1[i]->SetTitle(Form("%s",st[j].data()));
			th1[i]->Scale(1/th1[i]->Integral());
			th1[i]->SetLineColor(i+1);
			th1[i]->SetLineWidth(3);
			th1[i]->SetMarkerSize(0);
			
		}
		th1[0]->SetMaximum(th1[0]->GetMaximum()>th1[1]->GetMaximum()?th1[0]->GetMaximum()*1.2:th1[1]->GetMaximum()*1.2);
		for(int i=0;i<2;i++){
			if(i==0)th1[i]->Draw("");
			else th1[i]->Draw("same");
			leg->AddEntry(th1[i],Form("%s",legendName[i].data()));
		}
		
		leg->Draw("same");
		
		if(j==0)c1->Print("c.pdf(");
		else if(j==st.size()-2)c1->Print("c.pdf)");
		else c1->Print("c.pdf");
	}
  
  

	
}