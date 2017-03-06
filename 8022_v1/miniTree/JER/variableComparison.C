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
#include "../../../untuplizer.h"
#include <TClonesArray.h>
#include <fstream>
#include <cmath>
#include <TSystem.h>
#include <string>
#include <sstream>
#include "../../../setNCUStyle.C"


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
			Float_t  etadiff = data.GetFloat("etadiff");
			Float_t  dijetmass_softdrop_corr = data.GetFloat("dijetmass_softdrop_corr");
			Float_t  jet1bbtag = data.GetFloat("jet1bbtag");
			Float_t  jet2bbtag = data.GetFloat("jet2bbtag");
			Float_t  jet1_puppi_msoftdrop_TheaCorr = data.GetFloat("jet1_puppi_msoftdrop_TheaCorr");
			Float_t  jet2_puppi_msoftdrop_TheaCorr = data.GetFloat("jet2_puppi_msoftdrop_TheaCorr");
			if(jet1_puppi_msoftdrop_TheaCorr<105 ||jet1_puppi_msoftdrop_TheaCorr>135)continue;
			if(jet2_puppi_msoftdrop_TheaCorr<105 ||jet2_puppi_msoftdrop_TheaCorr>135)continue;
			//if(dijetmass_softdrop_corr<mjj*0.8||dijetmass_softdrop_corr>mjj*1.2)continue;
			if(jet1bbtag<0.3 ||jet2bbtag<0.3)continue;
			//if(jet1bbtag>0.8 &&jet2bbtag>0.8)continue;
			if(etadiff>1.3)continue;
			//Float_t  jet2_puppi_msoftdrop = data.GetFloat("jet2_puppi_msoftdrop");
			th1->Fill(vari);
	}
	return th1;
}

void setTGraph(TGraph* tg1,int i,bool setMax=0){
	
	tg1->SetLineColor(i+1);
	//tg1->SetLineColor(98-4*i);
	tg1->SetMarkerColor(i+1);
	if(i==4)tg1->SetLineColor(kOrange);
	if(i==4)tg1->SetMarkerColor(kOrange);
	//tg1->SetMarkerColor(98-4*i);
	tg1->SetMarkerStyle(20);
	
	tg1->GetXaxis()->SetTitle("m_{jj}[GeV]");
	tg1->GetYaxis()->SetTitle("");
	//limits->GetZaxis()->SetTitle("#sigma_{95\% CL} / #sigma_{th}");
	tg1->SetTitle("");
	tg1->SetLineWidth(2);
	tg1->SetFillColor(0);
	//tg1->SetMaximum(1);
	//tg1->SetMinimum(0.0014);
	tg1->SetMarkerSize(1);
	//tg1->SetMarkerStyle(20);
	//tg1->SetMinimum(0.05);
	tg1->GetYaxis()->SetTitleOffset(1.5);
	//tg1->GetZaxis()->SetTitleOffset(0.65);
	// size of axis labels
	tg1->GetXaxis()->SetTitleSize(0.04);
	tg1->GetYaxis()->SetTitleSize(0.04);
	//limits->GetZaxis()->SetTitleSize(0.035);
	tg1->GetXaxis()->SetLabelSize(0.05);
	tg1->GetYaxis()->SetLabelSize(0.03); 
}

void setLeg(TLegend* leg){
	leg->SetBorderSize(1);
	//leg->SetLineColor(0);                                                                                                                     
	leg->SetLineStyle(1);
	leg->SetLineWidth(1);
	leg->SetFillColor(0);
	leg->SetFillStyle(1001);
}

void variableComparison(){
	//TH1D* th1[3];
	//th1[0]=massPlotBase("BulkGrav_M-2000_0.root","jet1_puppi_msoftdrop");
	//th1[1]=massPlotBase("BulkGrav_M-2000_0.root","jet1_puppi_msoftdrop_TheaCorr");
	//th1[2]=massPlotBase("HCorr.root","jet1_puppi_msoftdrop_TheaCorr");
	
	TCanvas* c1;
	//setNCUStyle(true);
	c1 = new TCanvas("c1","",600,600);

	TLegend *leg = new TLegend(0.50, 0.73, 0.85, 0.88);
  
	setLeg(leg);
  
	vector<string> st;
	vector<double> up;
	vector<double> down;
	
	st.push_back("dijetmass_softdrop_corr"); up.push_back(1500);  down.push_back(700);  
	st.push_back("dijetmass_softdrop_corr_JERup"); up.push_back(1500);  down.push_back(700);  
	st.push_back("dijetmass_softdrop_corr_JERdown"); up.push_back(1500);  down.push_back(700);  
	
	int masspoint[10]={1000,1200,1400,1600,1800,2000,2500,3000,4000,4500};
	double masspointd[10]={1000,1200,1400,1600,1800,2000,2500,3000,4000,4500};
	TH1D* th1[10],*th2[10],*th3[10];
	double JERup[10],JERdown[10];
	double JERWidthup[10],JERWidthdown[10];
	
	for(unsigned int i=0;i<10;i++){
		//TFile* tf1;
		//tf1=TFile::Open();
		//th1[0]=massPlotBase("JetHT.root",st[j],down[j],up[j]);
		//th1[1]=massPlotBase("reRecoMinitree_v2.root",st[j],down[j],up[j]);
		
		th1[i]=massPlotBase(Form("BulkGrav_M-%d_0.root.root",masspoint[i]),st[0],masspoint[i]*0.7,masspoint[i]*1.3);
		th2[i]=massPlotBase(Form("BulkGrav_M-%d_0.root.root",masspoint[i]),st[1],masspoint[i]*0.7,masspoint[i]*1.3);
		th3[i]=massPlotBase(Form("BulkGrav_M-%d_0.root.root",masspoint[i]),st[2],masspoint[i]*0.7,masspoint[i]*1.3);
		
		JERup[i]=th2[i]->GetMean()/th1[i]->GetMean();
		JERWidthup[i]=th2[i]->GetRMS()/th1[i]->GetRMS();
		JERdown[i]=th3[i]->GetMean()/th1[i]->GetMean();
		JERWidthdown[i]=th3[i]->GetRMS()/th1[i]->GetRMS();
		
		cout<<"i="<<i<<","<<JERup[i]<<","<<JERdown[i]<<endl;
		cout<<"Width : i="<<i<<","<<JERWidthup[i]<<","<<JERWidthdown[i]<<endl;
		leg->Clear();
		
		th1[i]->Draw("");
		th2[i]->SetLineColor(2);
		th3[i]->SetLineColor(3);
		th2[i]->Draw("same");
		th3[i]->Draw("same");
		c1->Print(Form("Mjj/%d.pdf",masspoint[i]));
	}
	TGraph* tg1[4];
	tg1[0]=new TGraph(10,masspointd,JERup);
	tg1[1]=new TGraph(10,masspointd,JERdown);
	tg1[2]=new TGraph(10,masspointd,JERWidthup);
	tg1[3]=new TGraph(10,masspointd,JERWidthdown);
	for(int i=0;i<4;i++)setTGraph(tg1[i],0);
	
	leg->AddEntry(tg1[0],"up");
	leg->AddEntry(tg1[1],"down");
	tg1[0]->Draw("APL");
	tg1[1]->SetLineColor(2);
	tg1[0]->SetMinimum(0.999);
	tg1[0]->SetMaximum(1.0015);
	tg1[1]->Draw("samePL");
	leg->Draw("same");
	
	c1->Print("mean.pdf");
	
	tg1[2]->Draw("APL");
	tg1[2]->SetMinimum(0.985);
	tg1[2]->SetMaximum(1.02);
	tg1[3]->SetLineColor(2);
	tg1[3]->Draw("samePL");
	leg->Draw("same");
	c1->Print("width.pdf");
}