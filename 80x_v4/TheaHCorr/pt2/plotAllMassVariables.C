// example code to run Bulk Graviton->ZZ->ZlepZhad selections on electron-channel

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TString.h>
#include <map>
#include <TH1D.h>
#include <TFile.h>
#include <TF1.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TObject.h>
#include "../../setNCUStyle.C"


using namespace std;

float FWHM(TH1F* hist)
{
  int bin1 = hist->FindFirstBinAbove(hist->GetMaximum()/2);
  int bin2 = hist->FindLastBinAbove(hist->GetMaximum()/2);
  float fwhm = hist->GetBinCenter(bin2) - hist->GetBinCenter(bin1);
  // return (fwhm/2.36);
  return (fwhm     );
  //  return fwhm;
}


void plotAllMassVariablesBase(std::string inputFile){

  setNCUStyle();
  TString outputFile;
  outputFile=gSystem->GetFromPipe(Form("file=%s; test=${file%%.root}; echo \"${test}.pdf\"",inputFile.data()));
  cout << "output file name = " << outputFile.Data() << endl;

  TString HeaderName;
  HeaderName=gSystem->GetFromPipe(Form("file=%s; test=${file%%.root}; echo \"${test}\"",inputFile.data()));

  const int NTYPES=5;
  const int NHISTOS=3;
  TFile *inf = new TFile(inputFile.data());
  TH1F* hmass[NTYPES][NHISTOS];
  TH1F* hdiffmass[NTYPES][NHISTOS];
  //int COLORS[] ={1,2,4,8,kOrange-3,kCyan+2,kGreen+3};
  int COLORS[] ={kOrange,kGreen+2,1,2,4};


  std::string prefix[]={"leading","subleading","both"};
  std::string name[]={"PR","PRCorr","AK8SD","AK8SDCorrThea","AK8SDHCorr"};

  float max[3]={-9999,-9999,-9999};
  float maxdiff[3]={-9999,-9999,-9999};
  
  double hmassFit[4][NTYPES][NHISTOS];
  double hdiffmassFit[4][NTYPES][NHISTOS];
  gStyle->SetOptStat(0);
		gStyle->SetOptFit(0);
  for(int i=0; i < NTYPES;i++){
    for(int k=0; k < NHISTOS; k++){
	 TF1 *tf1[4];
	 	
	//hmass[i][k]->Sumw2();
	    
      hmass[i][k] = (TH1F*)inf->FindObjectAny(Form("h_%s_%s",name[i].data(),prefix[k].data()));
	hmass[i][k]->Sumw2();
      hmass[i][k]->Scale(1.0/hmass[i][k]->Integral());
      if( hmass[i][k]->GetMaximum()>max[k])
	max[k]=hmass[i][k]->GetMaximum();
tf1[0]=new TF1("fa1","gaus(25000)", hmass[i][k]->GetBinCenter(hmass[i][k]->GetMaximumBin())-20, hmass[i][k]->GetBinCenter(hmass[i][k]->GetMaximumBin())+20);
tf1[0] ->SetLineColor(COLORS[i]);
 hmass[i][k]->Fit(tf1[0],"","", hmass[i][k]->GetBinCenter(hmass[i][k]->GetMaximumBin())-20, hmass[i][k]->GetBinCenter(hmass[i][k]->GetMaximumBin())+20);
 hmass[i][k]->Fit(tf1[0],"","", hmass[i][k]->GetBinCenter(hmass[i][k]->GetMaximumBin())-20, hmass[i][k]->GetBinCenter(hmass[i][k]->GetMaximumBin())+20);
//tf1[0] ->SetLineColor(COLORS[i]);


 hmassFit[0][i][k]=tf1[0]->GetParameter(1);
 hmassFit[1][i][k]=tf1[0]->GetParError(1);
 hmassFit[2][i][k]=tf1[0]->GetParameter(2);
 hmassFit[3][i][k]=tf1[0]->GetParError(2);

	//hdiffmass[i][k]->Sumw2();

      hdiffmass[i][k] = (TH1F*)inf->FindObjectAny(Form("h_diff_%s_%s",name[i].data(),prefix[k].data()));			
hdiffmass[i][k]->Sumw2();	
      hdiffmass[i][k]->Scale(1.0/hdiffmass[i][k]->Integral());
      if( hdiffmass[i][k]->GetMaximum()>maxdiff[k])
	maxdiff[k]=hdiffmass[i][k]->GetMaximum();
tf1[1]=new TF1("fa1","gaus(25000)", hdiffmass[i][k]->GetBinCenter(hdiffmass[i][k]->GetMaximumBin())-0.2, hdiffmass[i][k]->GetBinCenter(hdiffmass[i][k]->GetMaximumBin())+0.2);
tf1[1]=new TF1("fa1","gaus(25000)", hdiffmass[i][k]->GetBinCenter(hdiffmass[i][k]->GetMaximumBin())-0.2, hdiffmass[i][k]->GetBinCenter(hdiffmass[i][k]->GetMaximumBin())+0.2);
		tf1[1] ->SetLineColor(COLORS[i]);
hdiffmass[i][k]->Fit(tf1[1],"","", hdiffmass[i][k]->GetBinCenter(hdiffmass[i][k]->GetMaximumBin())-0.2, hdiffmass[i][k]->GetBinCenter(hdiffmass[i][k]->GetMaximumBin())+0.2);

hdiffmassFit[0][i][k]=tf1[1]->GetParameter(1);
 hdiffmassFit[1][i][k]=tf1[1]->GetParError(1);
 hdiffmassFit[2][i][k]=tf1[1]->GetParameter(2);
 hdiffmassFit[3][i][k]=tf1[1]->GetParError(2);
    }
  }
  
  TCanvas* c1 = new TCanvas("c1","",500,500);

  TLegend* leg = new TLegend(0.179,0.692,0.326,0.882);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);



  int nPage=0;
  for(int k=0; k < NHISTOS; k++){
    for(int i=0; i < NTYPES;i++){

      hmass[i][k]->SetMaximum(max[k]*1.1);
      hmass[i][k]->SetLineWidth(3);
      hmass[i][k]->SetLineColor(COLORS[i]);
      hmass[i][k]->SetXTitle("Mass [GeV]");
      hmass[i][k]->SetTitle(k==2? Form("%s, %s jets",HeaderName.Data(),prefix[k].data()):
			    Form("%s, %s jet",HeaderName.Data(),prefix[k].data())
			    );
      if(i==0)
	hmass[i][k] ->Draw();
      else
	hmass[i][k] ->Draw("same");
      if(k==0)
	leg->AddEntry(hmass[i][k], name[i].data(),"l");
      if(i==NTYPES-1)
	leg->Draw("same");
      if(i==NTYPES-1 && nPage==0)
	{
	  c1->Print(Form("%s(",outputFile.Data()),"pdf");
	  nPage++;
	}
      else if(i==NTYPES-1)
	{
	  c1->Print(Form("%s",outputFile.Data()),"pdf");
	  nPage++;
	}
    }

    TLegend* leg2= new TLegend(0.657,0.263,0.904,0.889);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.035);
    leg2->SetBorderSize(0);

    for(int i=0; i < NTYPES;i++){

      hdiffmass[i][k]->SetMaximum(maxdiff[k]*1.5);
      hdiffmass[i][k]->SetLineWidth(3);
      hdiffmass[i][k]->SetLineColor(COLORS[i]);
      hdiffmass[i][k]->SetXTitle("(Mass-125)/125 [GeV]");
      hdiffmass[i][k]->SetTitle(k==2? Form("%s, %s jets",HeaderName.Data(),prefix[k].data()):
			    Form("%s, %s jet",HeaderName.Data(),prefix[k].data())
			    );
      if(i==0)
     	hdiffmass[i][k] ->Draw();
      else
     	hdiffmass[i][k] ->Draw("same");
      ofstream fout;
      fout.open(Form("%s_%s.dat",prefix[k].data(),name[i].data()),ios::out | ios::app);

      string tagname="Sigma = ";
      tagname += Form("%.3f",hdiffmassFit[2][i][k]);
      string tagname2="FWHM = ";
      tagname2 += Form("%.3f",FWHM(hdiffmass[i][k])); 
      string tagname3="Mean = ";
      tagname3 += Form("%.3f",hdiffmassFit[0][i][k]);

      //fout << hdiffmass[i][k]->GetMean() << " " << hdiffmass[i][k]->GetMeanError() << " " << hdiffmass[i][k]->GetRMS() << " " << hdiffmass[i][k]->GetRMSError()  << endl;
      fout << hdiffmassFit[0][i][k] << " " << hdiffmassFit[1][i][k] << " " << hdiffmassFit[2][i][k] << " " << hdiffmassFit[3][i][k]  << endl;
      fout.close();

      ofstream fout2;
      fout2.open(Form("rel_%s_%s.dat",prefix[k].data(),name[i].data()),ios::out | ios::app);
	 fout2 << hmassFit[0][i][k] << " " << hmassFit[1][i][k] << " " << hmassFit[2][i][k] << " " << hmassFit[3][i][k]  << endl;
      //fout2 << hmass[i][k]->GetMean() << " " << hmass[i][k]->GetMeanError() << " " << hmass[i][k]->GetRMS() << " " << hmass[i][k]->GetRMSError()  << endl;
      //fout2 << hmass[i][k]->GetMean() << " " << hmass[i][k]->GetMeanError() << " " << hmass[i][k]->GetRMS() << " " << hmass[i][k]->GetRMSError()  << endl;
      fout2.close();
 
      leg2->AddEntry(hdiffmass[i][k], name[i].data(),"l");
      leg2->AddEntry((TObject*)0, tagname3.data(),"");
      leg2->AddEntry((TObject*)0, tagname.data(),"");
      // leg2->AddEntry((TObject*)0, tagname2.data(),"");
      leg2->Draw("same");
      if(i==NTYPES-1)
	{
	  leg2->Draw("same");
	  if(k==2)
	    c1->Print(Form("%s)",outputFile.Data()),"pdf");
	  else
	    c1->Print(Form("%s",outputFile.Data()),"pdf");
	  nPage++;
	}
     }
  } // end loop over  jet type
 
 


}

void plotMultiGraphs(string masspoint,bool jetPtFix=0){

  std::string prefix[]={"leading","subleading","both"};
  std::string name[]={"PR","PRCorr","AK8SD","AK8SDCorrThea","AK8SDHCorr"};

  const int NTYPES=5;
  float mass[10]={300,400,500,600,700,800,900,1000,1250,1500};
 float mass2 [10]={452.286,540.273,627.553,714.669,799.949,883.213,1097.31,1310.44,1728.48,1937.44};
 if(jetPtFix){
	 for(int i=0;i<10;i++)mass[i]=mass2[i];
 }
 // int MARKERS[7]={20,21,22,23,34,29,24};
  int MARKERS[7]={23,34,20,22,21};
  //int COLORS[NTYPES]={1,4,2,kOrange,kGreen+2};
  // int COLORS[] ={1,2,4,8,kOrange-3,kCyan+2,kGreen+3};
   int COLORS[] ={kOrange,kGreen+2,1,2,4};
  setNCUStyle(true);

  
  
  TCanvas* c1 = new TCanvas("c1","",700,500);

  for(int i=0; i<3; i++){
	  
	 

    TGraphErrors* graph_mean[NTYPES];
    TGraphErrors* graph_RMS[NTYPES];
    TGraphErrors* graph_RMSMean[NTYPES];

    TMultiGraph *mg = new TMultiGraph();
    TMultiGraph *mg_h = new TMultiGraph();
    TMultiGraph *mg_a = new TMultiGraph();

    for(int j=0; j<NTYPES; j++)
      {
		 cout<<"here"<<i<<","<<j<<endl;
		
	ifstream fin;
	fin.open(Form("%s_%s.dat",prefix[i].data(),name[j].data()));
	
	float mean[9],meanerr[9],RMS[9],RMSerr[9];
	float rel[9], relerr[9];
	float masserr[9];
	for(int il=0; il<9; il++)
	  {
	    fin >> mean[il] >> meanerr[il] >> RMS[il] >> RMSerr[il];
	    masserr[il]=2;
	    //meanerr[il]/=mean[il];
	    //RMSerr[il]/=RMS[il];
	  }
	fin.close();

	int nPoint=9;
	if(jetPtFix)nPoint=8;
	graph_mean[j] = new TGraphErrors(nPoint,mass,mean,masserr,meanerr);
	graph_mean[j]->SetName(Form("gr_Mean_%d",j));
	graph_mean[j]->SetMarkerStyle(MARKERS[j]);
	graph_mean[j]->SetMarkerColor(COLORS[j]);
	graph_mean[j]->SetLineColor(COLORS[j]);
	graph_mean[j]->SetMarkerSize(1.1);
	mg->Add(graph_mean[j]);

	
	
	graph_RMS[j] = new TGraphErrors(nPoint,mass,RMS,masserr,RMSerr);
	graph_RMS[j]->SetName(Form("gr_RMS_%d",j));
	graph_RMS[j]->SetMarkerStyle(MARKERS[j]);
	graph_RMS[j]->SetMarkerColor(COLORS[j]);
	graph_RMS[j]->SetLineColor(COLORS[j]);
	graph_RMS[j]->SetMarkerSize(1.2);
	mg_h->Add(graph_RMS[j]);

	ifstream fin2;
	fin2.open(Form("rel_%s_%s.dat",prefix[i].data(),name[j].data()));
	
	for(int il=0; il<9; il++)
	  {
	    fin2 >> mean[il] >> meanerr[il] >> RMS[il] >> RMSerr[il];
	    float yield = RMS[il]/mean[il];
	    float err = yield*sqrt(pow(RMSerr[il]/RMS[il],2)+
				   pow(meanerr[il]/mean[il],2));
	    rel[il] = yield;
	    relerr[il] = err;
	    cout<<yield<<","<<err<<endl;
	  }
	fin2.close();

	
	graph_RMSMean[j] = new TGraphErrors(9,mass,rel,masserr,relerr);
	graph_RMSMean[j]->SetName(Form("gr_RMSMean_%d",j));
	graph_RMSMean[j]->SetMarkerStyle(MARKERS[j]);
	graph_RMSMean[j]->SetMarkerColor(COLORS[j]);
	graph_RMSMean[j]->SetLineColor(COLORS[j]);
	graph_RMSMean[j]->SetMarkerSize(1.1);
	mg_a->Add(graph_RMSMean[j]);

	
      } // end loop of mass types

    mg->SetTitle(i==2? Form("%s jets",prefix[i].data()):
		 Form("%s jet",prefix[i].data())
		 );


    mg->Draw("AP");
    mg->GetXaxis()->SetTitle("Jet Pt[GeV]");
    mg->GetYaxis()->SetTitleOffset(1.1);
    mg->GetYaxis()->SetTitle("Mean of (Mass-125)/125");
    mg->GetYaxis()->SetRangeUser(-0.3,0.2);
    
    TLegend* leg = new TLegend(0.648,0.704,0.397,0.877);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.05);
    leg->SetBorderSize(0);
    for(int itype=0; itype<NTYPES; itype++)
      {
	leg->AddEntry(Form("gr_Mean_%d",itype), Form("%s",name[itype].data()),"p");  
      }
    
    leg->Draw("same");

    string output = Form("MassMean_%s_%s",prefix[i].data(),masspoint.data());
    string final;
    final = output + ".gif";
    c1->Print(final.data());
    final = output + ".pdf";
    c1->Print(final.data());

    mg_h->SetTitle(i==2? Form("%s jets",prefix[i].data()):
		 Form("%s jet",prefix[i].data())
		 );
    
    mg_h->Draw("AP");
    mg_h->GetXaxis()->SetTitle("Jet Pt[GeV]");
    mg_h->GetYaxis()->SetTitleOffset(1.1);
    mg_h->GetYaxis()->SetTitle("#sigma of (Mass-125)/125");
    mg_h->GetYaxis()->SetRangeUser(0.05,0.22);
    
    leg->Draw("same");
    
    output = Form("MassRMS_%s_%s",prefix[i].data(),masspoint.data());
    final = output + ".gif";
    c1->Print(final.data());
    final = output + ".pdf";
    c1->Print(final.data());



    mg_a->SetTitle(i==2? Form("%s jets",prefix[i].data()):
		 Form("%s jet",prefix[i].data())
		 );
    
    mg_a->Draw("AP");
    mg_a->GetXaxis()->SetTitle("Jet Pt[GeV]");
    mg_a->GetYaxis()->SetTitleOffset(1.1);
    mg_a->GetYaxis()->SetTitle("#sigma/Mean of mass");
    mg_a->GetYaxis()->SetRangeUser(0.07,0.25);
    
    leg->Draw("same");
    
    output = Form("MassRMSMean_%s_%s",prefix[i].data(),masspoint.data());
    final = output + ".gif";
    c1->Print(final.data());
    final = output + ".pdf";
    c1->Print(final.data());

    
  } // end loop of jet types



}

void  plotAllMassVariables(string masspoint="",bool plotJetPtFix=0){
	if(plotJetPtFix){
		plotAllMassVariablesBase("0/B1000.root");
		plotAllMassVariablesBase("0/B1200.root");
		plotAllMassVariablesBase("0/B1400.root");
		plotAllMassVariablesBase("0/B1600.root");
		plotAllMassVariablesBase("0/B1800.root");
		plotAllMassVariablesBase("0/B2000.root");
		plotAllMassVariablesBase("0/B2500.root");
		plotAllMassVariablesBase("0/B3000.root");
		plotMultiGraphs("FixJetPt",1);
	}
	else{ 
	plotAllMassVariablesBase(Form("300/B%s.root",masspoint.data()));
	plotAllMassVariablesBase(Form("400/B%s.root",masspoint.data()));
	plotAllMassVariablesBase(Form("500/B%s.root",masspoint.data()));
	plotAllMassVariablesBase(Form("600/B%s.root",masspoint.data()));
	plotAllMassVariablesBase(Form("700/B%s.root",masspoint.data()));
	plotAllMassVariablesBase(Form("800/B%s.root",masspoint.data()));
	plotAllMassVariablesBase(Form("900/B%s.root",masspoint.data()));
	plotAllMassVariablesBase(Form("1000/B%s.root",masspoint.data()));
	plotAllMassVariablesBase(Form("1250/B%s.root",masspoint.data()));
	//plotAllMassVariablesBase(Form("1500/B%s.root",masspoint.data()));
	plotMultiGraphs(masspoint);
	
	}
}
/*
void  plotAllMassVariables(){
	plotAllMassVariablesmp("1000");
	plotAllMassVariablesmp("1200");
	plotAllMassVariablesmp("1400");
	plotAllMassVariablesmp("1600");
	plotAllMassVariablesmp("1800");
	plotAllMassVariablesmp("2000");
	plotAllMassVariablesmp("2500");
	plotAllMassVariablesmp("3000");
	plotAllMassVariablesmp("merge");
}
*/
