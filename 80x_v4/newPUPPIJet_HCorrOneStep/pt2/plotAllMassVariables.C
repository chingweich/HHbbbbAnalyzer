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
#include <TCanvas.h>
#include <TObject.h>
#include "setNCUStyle.C"

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


void plotAllMassVariables(std::string inputFile){

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
