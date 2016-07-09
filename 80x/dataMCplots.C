#include <string>
#include <iostream>
#include <TPad.h>
#include <TH1D.h>
#include <TFile.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <THStack.h>
#include <TLegend.h>
#include <TTree.h>
#include <TKey.h>
#include <TSystemDirectory.h>
#include "../setNCUStyle.C"
#include <TSystem.h>
//#include "../readHists.h"

double rangeUserUp=0,rangeUserDown=0;
int isSetRange=0;

void myPlot(vector< TH1D*> h_Z,
           vector< TH1D*> v_data,
	    TH1D* h_data, TH1D* h_bkg){

  h_data->Reset();
  for(unsigned int i=0;i<v_data.size();i++)h_data->Add(v_data[i]);
  
  TLegend *leg = new TLegend(0.73, 0.60, 0.90, 0.87);
  
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);

   h_bkg->Reset();
    THStack *h_stack = new THStack("h_stack", "");
  for(unsigned int i=0;i<h_Z.size();i++){
	  h_Z[i]->SetFillColor(97-3*i);
	  h_Z[i]->SetLineColor(kBlack);
	  h_bkg->Add(h_Z[i]);
	   h_stack->Add(h_Z[i]);
	  
  }
   leg->AddEntry(h_Z[0], "QCD700", "f");
   leg->AddEntry(h_Z[1], "QCD1000", "f");
   leg->AddEntry(h_Z[2], "QCD1500", "f");
   leg->AddEntry(h_Z[3], "QCD2000", "f");
  /*
leg->AddEntry(h_Zjets, "Z+Jets", "f");
  leg->AddEntry(h_TT, "t#bar{t}", "f");
  leg->AddEntry(h_WW, "WW", "f");
  leg->AddEntry(h_WZ, "WZ", "f");
  leg->AddEntry(h_ZZ, "ZZ", "f");
  leg->AddEntry(h_ZH, "ZH", "f");
*/

  h_data->SetLineColor(kBlack);
  h_data->SetMarkerStyle(8);
  h_data->SetMarkerSize(1.5);
  h_data->GetYaxis()->SetTitleOffset(1.3);
  h_data->GetXaxis()->SetTitle("");
  h_data->GetXaxis()->SetLabelOffset(999);
  h_data->GetXaxis()->SetLabelSize(0);
  
  int bmin=0,bmax=0;

		for (int k=1;k<25001;k++){
			bmin=k;
			if (h_data->GetBinContent(k)/h_data->GetMaximum()>0.02) break;
	}

		for (int k=25000;k>0;k--){
			bmax=k;
			if (h_data->GetBinContent(k)/h_data->GetMaximum()>0.02) break;
	}
	double width=h_data->GetBinWidth(1);
	rangeUserUp=(bmax-0.5)*width+h_data->GetBinCenter(1);
	rangeUserDown=(bmin-0.5)*width+h_data->GetBinCenter(1);
	if(isSetRange)h_data->GetXaxis()->SetRangeUser(rangeUserDown,rangeUserUp);
	//if(isSetRange)h_stack->GetXaxis()->SetRangeUser(rangeUserDown,rangeUserUp);
	//h_stack->GetXaxis()->SetRangeUser((bmin-0.5)*width+h_data->GetBinCenter(1),(bmax-0.5)*width+h_data->GetBinCenter(1));

  if( h_data->GetMaximum() < h_stack->GetMaximum() ){
	h_stack->SetMaximum(h_data->GetMaximum()*1.3);
    h_stack->Draw("histe");
    h_stack->GetHistogram()->GetYaxis()->SetTitle("Event Numbers");
    h_stack->GetHistogram()->GetYaxis()->SetTitleSize(h_data->GetYaxis()->GetTitleSize());
    h_stack->GetHistogram()->GetYaxis()->SetLabelSize(h_data->GetYaxis()->GetLabelSize());
    h_stack->GetHistogram()->GetYaxis()->SetTitleOffset(1.3);
    h_stack->GetHistogram()->GetXaxis()->SetTickLength(0);
    h_stack->GetHistogram()->GetXaxis()->SetLabelOffset(999);
    h_data->Draw("elsame");
  
  }
    
  else{

    h_data->GetYaxis()->SetTitle("Event Numbers");
    h_data->Draw("el");
    h_stack->Draw("histesame");
    h_data->Draw("elsame");

  }

    
  
  leg->AddEntry(h_data, "Data", "lp");
  leg->Draw();

  TLatex *lar = new TLatex();

  lar->SetNDC(kTRUE);
  lar->SetTextSize(0.04);
  lar->SetLineWidth(5);
  //lar->DrawLatex(0.14, 0.94, "CMS #it{#bf{2015}}");
  lar->DrawLatex(0.60, 0.94, "L = 4.513 fb^{-1} at #sqrt{s} = 13 TeV");

}

void myRatio(TH1D* h_data, TH1D *h_bkg){
cout<<h_bkg->Integral()/h_data->Integral()<<endl;
  TH1D* h_ratio = (TH1D*)h_bkg->Clone("h_ratio");
	
  h_ratio->Reset();

  Int_t nbin = h_ratio->GetNbinsX();
  Float_t ratio[nbin];
  Float_t error[nbin];
  Float_t numer_nbincontent[nbin];
  Float_t denom_nbincontent[nbin];
  Float_t numer_binerror[nbin];
  Float_t denom_binerror[nbin];

  for(Int_t i = 1; i <= nbin; i++){
	
	
    numer_nbincontent[i] = h_data->GetBinContent(i);
    denom_nbincontent[i] = h_bkg ->GetBinContent(i);
    numer_binerror[i]    = h_data->GetBinError(i);
    denom_binerror[i]    = h_bkg ->GetBinError(i);

    if( denom_nbincontent[i] <= 0 || numer_nbincontent[i] <= 0 ) continue;
    if( denom_binerror[i] <= 0 || numer_binerror[i] <= 0 ) continue;

    ratio[i] = (Float_t)numer_nbincontent[i]/denom_nbincontent[i];
    error[i] = (ratio[i])*sqrt(pow(numer_binerror[i]/numer_nbincontent[i],2)+pow(denom_binerror[i]/denom_nbincontent[i],2));

    h_ratio->SetBinContent(i,ratio[i]);
    h_ratio->SetBinError(i,error[i]);

  }

  h_ratio->SetLineColor(kBlack);
  h_ratio->SetMarkerStyle(8);
  h_ratio->SetMarkerSize(1.5);
  h_ratio->SetTitle("");
  h_ratio->GetYaxis()->SetTitle("Data/MC");
  h_ratio->GetYaxis()->SetTitleOffset(0.45);
  h_ratio->GetXaxis()->SetLabelSize(0.1);
  h_ratio->GetXaxis()->SetLabelOffset(0.005);
  h_ratio->GetXaxis()->SetTitleSize(0.125);
  h_ratio->GetXaxis()->SetTitleOffset(0.8);
  h_ratio->GetYaxis()->SetLabelSize(0.1);
  h_ratio->GetYaxis()->SetTitleSize(0.1);
  h_ratio->GetYaxis()->SetNdivisions(505);
  h_ratio->GetYaxis()->SetRangeUser(0,2);
  h_ratio->Draw();

  Float_t x0 = h_bkg->GetXaxis()->GetXmin();
  Float_t x1 = h_bkg->GetXaxis()->GetXmax();
  Float_t y0 = 1.;
  Float_t y1 = 1.;

  TLine* one = new TLine(x0,y0,x1,y1);

  one->SetLineColor(2);
  one->SetLineStyle(1);
  one->SetLineWidth(2);
  one->Draw("same");
 if(isSetRange)h_ratio->GetXaxis()->SetRangeUser(rangeUserDown,rangeUserUp);
  h_ratio->Draw("same");

}

void dataMCplots(){

  setNCUStyle(true);
 
  Float_t up_height     = 0.8;
  Float_t dw_correction = 1.455;
  Float_t dw_height     = (1-up_height)*dw_correction;

  TCanvas c("c","",0,0,1000,900);
  c.Divide(1,2);

  TPad* c_up = (TPad*) c.GetListOfPrimitives()->FindObject("c_1");
  TPad* c_dw = (TPad*) c.GetListOfPrimitives()->FindObject("c_2"); 

  c_up->SetPad(0,1-up_height,1,1);
  c_dw->SetPad(0,0,1,dw_height);
  c_dw->SetBottomMargin(0.25);

  // To get the name of histograms
  
  TFile* tf1[7];
  tf1[0]=TFile::Open("sf/QCD700.root");
  tf1[1]=TFile::Open("sf/QCD1000.root");
  tf1[2]=TFile::Open("sf/QCD1500.root");
  tf1[3]=TFile::Open("sf/QCD2000.root");
  tf1[4]=TFile::Open("sf/data.root");
  
  vector<std::string> h_name;
  
  for(int i=0;i<2;i++){
		for(int j=0;j<2;j++){
			for(int k=0;k<5;k++){
				h_name.push_back(Form("Pt_j%d_sj%d_%db",i,j,k));  
				h_name.push_back(Form("Eta_j%d_sj%d_%db",i,j,k));  
				h_name.push_back(Form("subCSV_j%d_sj%d_%db",i,j,k));  
				
			}
		}
		for(int k=0;k<5;k++){
			h_name.push_back(Form("deltaR_j%d_%db",i,k));  
			h_name.push_back(Form("Pt_j%d_%db",i,k));  
			h_name.push_back(Form("Eta_j%d_%db",i,k));  
			h_name.push_back(Form("prMassL2L3_j%d_%db",i,k));  
			h_name.push_back(Form("tau21_j%d_%db",i,k));  
			h_name.push_back(Form("SDMassL2L3_j%d_%db",i,k));  
			h_name.push_back(Form("puppiTau21_j%d_%db",i,k));  
		
		}
	}
	for(int k=0;k<5;k++){
		h_name.push_back(Form("totalMass_%db",k));  
		h_name.push_back(Form("deltaEta_%db",k));  
	}
  
 // h_name.push_back("cutflow");  
h_name.push_back("Nbtagjet");  
  for(unsigned int i = 0; i < h_name.size(); i++){
	 
	  
	//cout<<h_name[i]<<endl;
	TH1D *th1[10];
	for(int k=0;k<5;k++)th1[k]=(TH1D* )tf1[k]->FindObjectAny(Form("%s",h_name[i].data()));
	if(h_name[i].find("3b")!= std::string::npos)continue;	
	TString endfix;
	if(h_name[i].find("4b")!= std::string::npos){
		endfix=gSystem->GetFromPipe(Form("file=%s; test=${file%%*_4b}; echo \"${test}\"",h_name[i].data()));
		//cout<<endfix<<endl;
		for(int j=0;j<5;j++){
				//cout<<th1[j]->Integral()<<endl;
				TH1D *th2=(TH1D* )tf1[j]->FindObjectAny(Form("%s_0b",endfix.Data()));
				th1[j]->Add(th2);
				//cout<<th1[j]->Integral()<<"1,"<<endl;
				th2=(TH1D* )tf1[j]->FindObjectAny(Form("%s_1b",endfix.Data()));
				th1[j]->Add(th2);
				//cout<<th1[j]->Integral()<<"2,"<<endl;
				th2=(TH1D* )tf1[j]->FindObjectAny(Form("%s_2b",endfix.Data()));
				th1[j]->Add(th2);
				//cout<<th1[j]->Integral()<<"3,"<<endl;
				th2=(TH1D* )tf1[j]->FindObjectAny(Form("%s_3b",endfix.Data()));
				th1[j]->Add(th2);
				//cout<<"4btest="<<th1[j]->Integral()<<endl;
		}
	}
	
	
	
	
    TH1D *h_data = (TH1D* )th1[0]->Clone("h_data");
    TH1D *h_bkg  = (TH1D* )th1[0]->Clone("h_bkg");
	TH1D *temp = (TH1D* )th1[0]->Clone("h_bkg");
	 c_up->cd();
	// if(h_name[i].find("cutflow")!= std::string::npos)c_up->SetLogy();	
	  if(h_name[i].find("Pt")!= std::string::npos)isSetRange=1;	
	vector<TH1D* > v2;
	vector<TH1D* > vd;
	double fixNumber=13445/22434.4;//7306/11857.3;
	double Xsec[4]={6831,1207,119.9,25.24};
	double scaleTemp[4];
	for(int k=0;k<4;k++){
		TH1D *th2=(TH1D* )tf1[k]->FindObjectAny("cutflow");
		
		th1[k]->Scale(fixNumber*4513.309419004*Xsec[k]/th2->GetBinContent(1));
		scaleTemp[k]=fixNumber*4513.309419004*Xsec[k]/th2->GetBinContent(1);
		v2.push_back(th1[k]);
		
	}	
	vd.push_back(th1[4]);
	//cout<<h_name[i]<<endl;
    myPlot(v2,
vd,
	   h_data, h_bkg);

    c_up->RedrawAxis();
    c_dw->cd();

    myRatio(h_data, h_bkg);
	
	
	if(h_name[i].find("Pt")!= std::string::npos)isSetRange=0;	
   // c.Draw();
	if(h_name[i].find("0b")!= std::string::npos)c.Print(Form("dataMC/0b/%s.pdf",h_name[i].data()));
	else if(h_name[i].find("1b")!= std::string::npos)c.Print(Form("dataMC/1b/%s.pdf",h_name[i].data()));
	else if(h_name[i].find("2b")!= std::string::npos)c.Print(Form("dataMC/2b/%s.pdf",h_name[i].data()));
    else if(h_name[i].find("4b")!= std::string::npos)c.Print(Form("dataMC/all/%s.pdf",endfix.Data()));
	else c.Print(Form("dataMC/all/%s.pdf",h_name[i].data()));
	for(int k=0;k<4;k++)th1[k]->Scale(1/scaleTemp[k]);
  }

}