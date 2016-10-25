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
#include <string>
#include "TGaxis.h"
#include "TPad.h"
#include <fstream>

#define  nWidth 9
#define  nBmin 5

#define masspointConst 9

double  getEff(TH1F *h,int bmin,int bmax);
double getSign(double d1,double d2);
double  getErr(TH1F *h,int bmin,int bmax);
double  getErr2(TH1F *h,int bmin,int bmax);

using namespace std;

TFile *f,*f2,*f3;
TCanvas* c1,*c2;
double x1NDC = 0.6522;
  double y1NDC = 0.7835;
  double x2NDC = 0.9122;
  double y2NDC = 0.9860;

void setLeg(TLegend *leg){
leg->SetFillColor(18);
leg->SetFillStyle(0);
leg->SetTextSize(0.03);
leg->SetBorderSize(2);
}
/*
void makeEff(TString fin){

  c1 = new TCanvas("c1","",1360,768);
  f= TFile::Open(fin.Data());
  TH1F * th1 = (TH1F*)f->FindObjectAny("HMass");
  int width [nWidth]={30,35,50,45,50};
  int bmin[nBmin]={91,93,95,97,99,101,103,105,107,109};
  double eff[nWidth][nBmin];
  for(int i=0;i<nWidth;i++){
    for(int j=0;j<nBmin;j++){
      eff[i][j]=getEff(th1,bmin[j],bmin[j]+width[i]);
      cout<<"range={"<<bmin[j]<<","<<bmin[j]+width[i]<<"} and efficiency="<<eff[i][j]<<endl;
    }
  }
}
*/
void makeEff0929(){
  c1 = new TCanvas("c1","",1360,768);
  c2 = new TCanvas("c2","",1360,768);
  c1->cd();
  string output="FATjetATLASmassL2L3Corr";
  int twikiSign[masspointConst][nWidth][nBmin];
  int twikiSignNum[masspointConst][nWidth][nBmin];
  double twikiSignValue[masspointConst][nWidth][nBmin];
  double twikiSignNumValue[masspointConst][nWidth][nBmin];
  double twikiEffValue[masspointConst][nWidth][nBmin];
  double tgraphMass[masspointConst]={1200,1400,1600,1800,2000,2500,3000,4000,4500};
  double tgraphMassError[masspointConst]={0};
  double tgraphSigEff[masspointConst]={0};
  double tgraphSigEffError[masspointConst]={0};
  double tgraphBkgEff[masspointConst]={0};
  double tgraphBkgEffError[masspointConst]={0};
  double tgraphSgnEff[masspointConst]={0};
  double tgraphSgnEffError[masspointConst]={0};
  double tgraphSigEff145[masspointConst]={0};
  double tgraphSigEffError145[masspointConst]={0};
  double tgraphBkgEff145[masspointConst]={0};
  double tgraphBkgNum145[masspointConst]={0};
  double tgraphBkgNum[masspointConst]={0};
  double tgraphBkgEffError145[masspointConst]={0};
  double tgraphSgnEff145[masspointConst]={0};
  double tgraphSgnEffError145[masspointConst]={0};
  //double twikiEffNumValue[13][nWidth][nBmin];
  //double twikiWidth[13][nWidth][nBmin];
  //double twikiBmin[13][nWidth][nBmin];
  //string  masspoint[4]={"1000","2000","3000","4500"};
  string  masspoint[masspointConst]={"1200","1400","1600","1800","2000","2500","3000","4000","4500"};
  int width [nWidth]={20,25,30,35,40,45,50,55,60};
  int bmin[nBmin]={90,95,100,105,110};
  //int bmin[nBmin]={0,5,10,15,20};
  for (int massP=0;massP<masspointConst;massP++){
  TString fin =Form("category/B%s.root",masspoint[massP].data()),
    fin2 = "category/QCD700_1.root";
  // fin2=Form("root_files/DYHTall-%s.root",masspoint[massP].data());
  f= TFile::Open(fin.Data());
  f2= TFile::Open(fin.Data());
  TH1F * th1 = (TH1F*)f->FindObjectAny("reg2Mass_j0");
  TH1F * th2 = (TH1F*)f2->FindObjectAny("reg2Mass_j0");
  
  //TString fin3 = Form("root_files/BulkGravitonZlepZqq-%s.root",masspoint[massP].data());
  //f3= TFile::Open(fin3.Data());
  //TH1F * th3 = (TH1F*)f3->FindObjectAny("HMass");
  //cout<<"before"<<th2->GetEntries()<<endl;
  th2->Sumw2();
  //cout<<"here";
  //th2->Add(th3);
  //cout<<"after"<<th2->GetEntries()<<endl;
  double eff[nWidth][nBmin],eff2[nWidth][nBmin],err[nWidth][nBmin],err2[nWidth][nBmin],err2Num[nWidth][nBmin];
  double sign[nWidth][nBmin],signCP[nWidth][nBmin],signErr[nWidth][nBmin];
  double signNum[nWidth][nBmin],signNumCP[nWidth][nBmin],signNumErr[nWidth][nBmin];
  for(int i=0;i<nWidth;i++){
    for(int j=0;j<nBmin;j++){
      eff[i][j]=getEff(th1,bmin[j],bmin[j]+width[i]);
      eff2[i][j]=getEff(th2,bmin[j],bmin[j]+width[i]);
	  //if(eff2[i][j]==0)cout<<massP<<","<<eff2[i][j]<<endl;
	  err[i][j]=getErr(th1,bmin[j],bmin[j]+width[i]);
      err2[i][j]=getErr(th2,bmin[j],bmin[j]+width[i]);
	  err2Num[i][j]=getErr2(th2,bmin[j],bmin[j]+width[i]);
	  //if(massP==0)cout<<"method1="<<err2[i][j]<<",method2="<<err2Num[i][j]/th2->Integral()<<endl;
      sign[i][j]=getSign(eff[i][j],eff2[i][j]);
	  signNum[i][j]=getSign(eff[i][j],eff2[i][j]*th2->Integral());
      signCP[i][j]=sign[i][j];
      signNumCP[i][j]=signNum[i][j];
	  signErr[i][j]=(1/(1+sqrt(eff2[i][j])))*sqrt(err[i][j]*err[i][j]+ (eff[i][j]*eff[i][j]*err2[i][j]*err2[i][j])/(4*eff2[i][j]*(1+sqrt(eff2[i][j]))*(1+sqrt(eff2[i][j]))));
      //signNumErr[i][j]=(1/(1+sqrt(eff2[i][j]*th2->Integral())))*sqrt(err[i][j]*err[i][j]+ (eff[i][j]*eff[i][j]*err2[i][j]*err2[i][j]*th2->Integral())/(4*eff2[i][j]*(1+sqrt(eff2[i][j]*th2->Integral()))*(1+sqrt(eff2[i][j]*th2->Integral()))));
	  signNumErr[i][j]=(1/(1+sqrt(eff2[i][j]*th2->Integral())))*sqrt(err[i][j]*err[i][j]+ (eff[i][j]*eff[i][j]*err2Num[i][j]*err2Num[i][j])/(4*eff2[i][j]*th2->Integral()*(1+sqrt(eff2[i][j]*th2->Integral()))*(1+sqrt(eff2[i][j]*th2->Integral()))));
		if(eff2[i][j]<1e-7){
			signErr[i][j]=err[i][j]/(1+sqrt(eff2[i][j]));
			signNumErr[i][j]=err[i][j]/(1+sqrt(eff2[i][j]*th2->Integral()));
		}
		if(massP==2 &&i==3 && j==3){
			tgraphBkgNum145[massP]=eff2[i][j]*th2->Integral();
		tgraphSigEff145[massP]=eff[i][j];
		tgraphSigEffError145[massP]=err[i][j];
		tgraphBkgEff145[massP]=eff2[i][j];
		tgraphBkgEffError145[massP]=err2[i][j];
		tgraphSgnEff145[massP]=signNum[i][j];
		tgraphSgnEffError145[massP]=signNumErr[i][j];
		}
		else if (massP==7 &&i==2 && j==2){
			tgraphBkgNum145[massP]=eff2[i][j]*th2->Integral();
		tgraphSigEff145[massP]=eff[i][j];
		tgraphSigEffError145[massP]=err[i][j];
		tgraphBkgEff145[massP]=eff2[i][j];
		tgraphBkgEffError145[massP]=err2[i][j];
		tgraphSgnEff145[massP]=signNum[i][j];
		tgraphSgnEffError145[massP]=signNumErr[i][j];
		}
		else if (massP==9 &&i==2 && j==2){
			tgraphBkgNum145[massP]=eff2[i][j]*th2->Integral();
		tgraphSigEff145[massP]=eff[i][j];
		tgraphSigEffError145[massP]=err[i][j];
		tgraphBkgEff145[massP]=eff2[i][j];
		tgraphBkgEffError145[massP]=err2[i][j];
		tgraphSgnEff145[massP]=signNum[i][j];
		tgraphSgnEffError145[massP]=signNumErr[i][j];
			
		}
		else if (massP==11 &&i==4 && j==1){
			
			tgraphBkgNum145[massP]=eff2[i][j]*th2->Integral();
		tgraphSigEff145[massP]=eff[i][j];
		tgraphSigEffError145[massP]=err[i][j];
		tgraphBkgEff145[massP]=eff2[i][j];
		tgraphBkgEffError145[massP]=err2[i][j];
		tgraphSgnEff145[massP]=signNum[i][j];
		tgraphSgnEffError145[massP]=signNumErr[i][j];
		}
		
	    else {}
		
    }
  }
  
  vector<int> signI,signJ;
  for(int ij=0;ij<nWidth*nBmin;ij++){
    double tempSign=0;
    int tempI=0,tempJ=0;
    for(int i=0;i<nWidth;i++){
      for(int j=0;j<nBmin;j++){
        if(signCP[i][j]>tempSign){
          cout<<"i="<<i<<",j="<<j<<",signcp="<<signCP[i][j]<<",temp="<<tempSign<<endl;
          tempSign=signCP[i][j];
	  tempI=i;
	  tempJ=j;
	}
      }
    }

    //cout<<signCP[tempI][tempJ]<<endl;
    signCP[tempI][tempJ]=0;
	if(width[tempI]+bmin[tempJ]>150){
		
		continue;
	}
	//if(massP==0)cout<<"tempI"<<tempI<<",tempJ"<<tempJ<<endl;
	signI.push_back(tempI);
    signJ.push_back(tempJ);
  }
  
  vector<int> signINum,signJNum;
  for(int ij=0;ij<nWidth*nBmin;ij++){
    double tempSign=0;
    int tempI=0,tempJ=0;
    for(int i=0;i<nWidth;i++){
      for(int j=0;j<nBmin;j++){
        if(signNumCP[i][j]>tempSign){
          //cout<<"i="<<i<<",j="<<j<<",signcp="<<signCP[i][j]<<",temp="<<tempSign<<endl;
          tempSign=signNumCP[i][j];
	  tempI=i;
	  tempJ=j;
	}
      }
    }
    
    //cout<<signCP[tempI][tempJ]<<endl;
    signNumCP[tempI][tempJ]=0;
    if(width[tempI]+bmin[tempJ]>150){
		//cout<<"!";
		continue;
	}
	signINum.push_back(tempI);
    signJNum.push_back(tempJ);
  }
  
  int numberBin=15;
  int numberBinFull=35;
  TH1F * th1Sign=new TH1F("Sign","Sign",numberBin,0,numberBin);
  TH1F * th1Eff=new TH1F("th1Eff","th1Eff",numberBin,0,numberBin);
  TH1F * th1Eff2=new TH1F("th1Eff2","th1Eff2",numberBin,0,numberBin);
  TH1F * th1SignFull=new TH1F("SignFull","SignFull",numberBinFull,0,numberBinFull);
  TH1F * th1EffFull=new TH1F("th1EffFull","th1EffFull",numberBinFull,0,numberBinFull);
  TH1F * th1Eff2Full=new TH1F("th1Eff2Full","th1Eff2Full",numberBinFull,0,numberBinFull);
  for(int i=0;i<numberBin;i++){
	  th1Sign->SetBinContent(i+1,sign[signI[i]][signJ[i]]);
	  th1Sign->SetBinError(i+1,signErr[signI[i]][signJ[i]]);
	  th1Eff->SetBinContent(i+1,eff[signI[i]][signJ[i]]);
	  th1Eff2->SetBinContent(i+1,eff2[signI[i]][signJ[i]]);
	  th1Eff->SetBinError(i+1,err[signI[i]][signJ[i]]);
	  th1Eff2->SetBinError(i+1,err2[signI[i]][signJ[i]]);
      th1Sign->GetXaxis()->SetBinLabel(i+1,Form("%d-%d",bmin[signJ[i]],bmin[signJ[i]]+width[signI[i]]));
    //cout<<"range=["<<bmin[signJ[i]]<<","<<bmin[signJ[i]]+width[signI[i]]<<"] and significant="<<sign[signI[i]][signJ[i]]<<endl;
  }
  
  //cout<<"signI.size="<<signI.size()<<endl;
  for(unsigned int i=0;i<signI.size();i++){
	  twikiSign[massP][signI[i]][signJ[i]]=i+1;
	  
	  th1SignFull->SetBinContent(i+1,sign[signI[i]][signJ[i]]);
	  th1SignFull->SetBinError(i+1,signErr[signI[i]][signJ[i]]);
	  th1EffFull->SetBinContent(i+1,eff[signI[i]][signJ[i]]);
	  th1Eff2Full->SetBinContent(i+1,eff2[signI[i]][signJ[i]]);
	  th1EffFull->SetBinError(i+1,err[signI[i]][signJ[i]]);
	  th1Eff2Full->SetBinError(i+1,err2[signI[i]][signJ[i]]);
      //th1SignFull->GetXaxis()->SetBinLabel(i+1,Form("%d-%d",bmin[signJ[i]],bmin[signJ[i]]+width[signI[i]]));
      th1SignFull->GetXaxis()->SetBinLabel(i+1,"");
	  twikiSignValue[massP][signI[i]][signJ[i]]=sign[signI[i]][signJ[i]];
	  twikiEffValue[massP][signI[i]][signJ[i]]=eff[signI[i]][signJ[i]];
    //cout<<"range=["<<bmin[signJ[i]]<<","<<bmin[signJ[i]]+width[signI[i]]<<"] and significant="<<sign[signI[i]][signJ[i]]<<endl;
  }
  gStyle->SetOptStat(0000000000);
  
  
  th1Sign->SetTitle(Form("largest 15 Sign.(eff),%sGev",masspoint[massP].data()));
  th1Sign->SetMinimum(0);
  th1Sign->SetMaximum(1);
  th1Sign->SetLineColor(4);
  th1Sign->Draw();
  th1Eff->SetLineColor(2);
  th1Eff2->SetLineColor(3);
  th1Eff->Draw("same");
  th1Eff2->Draw("same");
  TLegend* leg ;
  leg=new TLegend(x1NDC,y1NDC,x2NDC,y2NDC);
  setLeg(leg);
  leg->AddEntry(th1Sign,"significance");
  leg->AddEntry(th1Eff,"signal efficiency");
  leg->AddEntry(th1Eff2,"Bkg. efficiency");
  leg->Draw("same");
  
  if(massP==0)c1->Print(Form("pdfPR/%s.pdf(",output.data()));
  //else if(massP==12)c1->Print("pdf/signv2.pdf)");
  else c1->Print(Form("pdfPR/%s.pdf",output.data()));
  
  th1SignFull->SetTitle(Form("all windows(eff),%sGev",masspoint[massP].data()));
  th1SignFull->SetMinimum(0);
  th1SignFull->SetMaximum(1);
  th1SignFull->SetLineColor(4);
  th1SignFull->Draw("Hist");
  th1EffFull->SetLineColor(2);
  th1Eff2Full->SetLineColor(3);
  th1EffFull->Draw("same");
  th1Eff2Full->Draw("same");
  leg->Draw("same");
  
  c1->Print(Form("pdfPR/%s.pdf",output.data()));
  
  for(int i=0;i<numberBin;i++){
	  if(i==0){
		tgraphSigEff[massP]=eff[signINum[i]][signJNum[i]];
		tgraphSigEffError[massP]=err[signINum[i]][signJNum[i]];
		tgraphBkgEff[massP]=eff2[signINum[i]][signJNum[i]];
		tgraphBkgEffError[massP]=err2[signINum[i]][signJNum[i]];
		tgraphSgnEff[massP]=signNum[signINum[i]][signJNum[i]];
		tgraphSgnEffError[massP]=signNumErr[signINum[i]][signJNum[i]];
		tgraphBkgNum[massP]=eff2[signINum[i]][signJNum[i]]*th2->Integral();
		
	  
	  }
	  th1Sign->SetBinContent(i+1,signNum[signINum[i]][signJNum[i]]);
	  th1Sign->SetBinError(i+1,signNumErr[signINum[i]][signJNum[i]]);
	  th1Eff->SetBinContent(i+1,eff[signINum[i]][signJNum[i]]);
	  th1Eff2->SetBinContent(i+1,eff2[signINum[i]][signJNum[i]]);
	  th1Eff->SetBinError(i+1,err[signINum[i]][signJNum[i]]);
	  th1Eff2->SetBinError(i+1,err2[signINum[i]][signJNum[i]]);
      th1Eff->GetXaxis()->SetBinLabel(i+1,Form("%d-%d",bmin[signJNum[i]],bmin[signJNum[i]]+width[signINum[i]]));
    //cout<<"range=["<<bmin[signJ[i]]<<","<<bmin[signJ[i]]+width[signI[i]]<<"] and significant="<<sign[signI[i]][signJ[i]]<<endl;
  }
  
  
  for(unsigned int i=0;i<signINum.size();i++){
	  twikiSignNum[massP][signINum[i]][signJNum[i]]=i+1;
	  
	  th1SignFull->SetBinContent(i+1,signNum[signINum[i]][signJNum[i]]);
	  th1SignFull->SetBinError(i+1,signNumErr[signINum[i]][signJNum[i]]);
	  th1EffFull->SetBinContent(i+1,eff[signINum[i]][signJNum[i]]);
	  th1Eff2Full->SetBinContent(i+1,eff2[signINum[i]][signJNum[i]]);
	  th1EffFull->SetBinError(i+1,err[signINum[i]][signJNum[i]]);
	  th1Eff2Full->SetBinError(i+1,err2[signINum[i]][signJNum[i]]);
      //th1SignFull->GetXaxis()->SetBinLabel(i+1,Form("%d-%d",bmin[signJ[i]],bmin[signJ[i]]+width[signI[i]]));
      th1SignFull->GetXaxis()->SetBinLabel(i+1,"");
      th1EffFull->GetXaxis()->SetBinLabel(i+1,"");
      th1Eff2Full->GetXaxis()->SetBinLabel(i+1,"");
	  twikiSignNumValue[massP][signINum[i]][signJNum[i]]=signNum[signINum[i]][signJNum[i]];
	  th1EffFull->GetXaxis()->SetBinLabel(i+1,Form("%d-%d",bmin[signJNum[i]],bmin[signJNum[i]]+width[signINum[i]]));
	  //twikiEffNumValue[massP][signINum[i]][signJNum[i]]=effNum[signINum[i]][signJNum[i]];
    //cout<<"range=["<<bmin[signJ[i]]<<","<<bmin[signJ[i]]+width[signI[i]]<<"] and significant="<<sign[signI[i]][signJ[i]]<<endl;
  }
  
  
  
  th1Eff->SetMinimum(0);
  th1Eff->SetMaximum(1);
  th1EffFull->SetMinimum(0);
  th1EffFull->SetMaximum(1);
  th1Eff->SetYTitle("efficiency");
  th1EffFull->SetYTitle("efficiency");
  th1Eff->SetTitle(Form("largest 15 Sign.(Num),%sGev",masspoint[massP].data()));
  th1Eff->Draw();
  th1Eff2->Draw("same");
  
  Float_t rightmax = 2*th1Sign->GetBinContent(1);
  //cout<<rightmax<<endl;
  Float_t scale = gPad->GetUymax()/rightmax;
  //hint1->SetLineColor(kRed);
  th1Sign->Scale(scale);
  //hint1->Draw("same");
  //draw an axis on the right side
  c1->Update();
  th1Sign->Draw("same");
  //leg->Clear();
  //leg->AddEntry(th1Eff,"signal efficiency");
  //leg->AddEntry(th1Eff2,"Bkg. efficiency");
  //leg->SetY1(0.8335);
  leg->Draw("same");
  
  TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
         gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
   axis->SetTitle("significance");
   axis->SetTitleColor(4);
   axis->SetLineColor(4);
   axis->SetLabelColor(4);
   axis->Draw();
  
  c1->Print(Form("pdfPR/%s.pdf",output.data()));
  
  
  
  
  th1EffFull->SetTitle(Form("all windows(Num),%sGev",masspoint[massP].data()));
  th1EffFull->Draw();
  th1Eff2Full->Draw("same");
  c1->Update();
  rightmax = 1.1*th1SignFull->GetBinContent(1);
  scale = gPad->GetUymax()/rightmax;
  th1SignFull->Scale(scale);
  
  th1SignFull->Draw("Hist,same");
  leg->Draw("same");
  
  TGaxis *axis2 = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
         gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
   axis2->SetTitle("significance");
   axis2->SetTitleColor(4);
   axis2->SetLineColor(4);
   axis2->SetLabelColor(4);
   axis2->Draw();
  
 
  
  if(massP==12)c1->Print(Form("pdfPR/%s.pdf)",output.data()));
  else c1->Print(Form("pdfPR/%s.pdf",output.data()));
  
  }
  
  
  TLegend* leg ;
  leg=new TLegend(0.1441452,0.742447,0.300645,0.883966);
  setLeg(leg);
  c2->cd();
  TGraphErrors* tg1= new TGraphErrors(masspointConst,tgraphMass,tgraphSigEff,tgraphMassError,tgraphSigEffError);
  TGraphErrors* tg2= new TGraphErrors(masspointConst,tgraphMass,tgraphBkgEff,tgraphMassError,tgraphBkgEffError);
  TGraphErrors* tg3= new TGraphErrors(masspointConst,tgraphMass,tgraphSgnEff,tgraphMassError,tgraphSgnEffError);
  
  for(int i=0;i<masspointConst;i++){
	  //cout<<tgraphSigEff[i]<<","<<tgraphSigEff145[i];
	  cout<<"i="<<i<<","<<tgraphSigEff[i]<<","<<tgraphBkgNum[i]<<","<<tgraphSgnEff[i]<<endl;
	  cout<<"i="<<i<<","<<tgraphSigEff145[i]<<","<<tgraphBkgNum145[i]<<","<<tgraphSgnEff145[i]<<endl;
	  tgraphSigEff[i]=tgraphSigEff[i]/tgraphSigEff145[i];
	  //cout<<tgraphSigEff[i]<<endl;
	  tgraphSigEffError[i]=sqrt((tgraphSigEffError[i]/tgraphSigEff[i])*(tgraphSigEffError[i]/tgraphSigEff[i])+(tgraphSigEffError145[i]/tgraphSigEff145[i])*(tgraphSigEffError145[i]/tgraphSigEff145[i]));
	  tgraphSigEffError[i]*=tgraphSigEff[i];
	  tgraphBkgEff[i]=tgraphBkgEff[i]/tgraphBkgEff145[i];
	  tgraphBkgEffError[i]=sqrt((tgraphBkgEffError[i]/tgraphBkgEff[i])*(tgraphBkgEffError[i]/tgraphBkgEff[i])+(tgraphBkgEffError145[i]/tgraphBkgEff145[i])*(tgraphBkgEffError145[i]/tgraphBkgEff145[i]));
	  tgraphBkgEffError[i]*=tgraphBkgEff[i];
	  tgraphSgnEff[i]=tgraphSgnEff[i]/tgraphSgnEff145[i];
	  tgraphSgnEffError[i]=sqrt((tgraphSgnEffError[i]/tgraphSgnEff[i])*(tgraphSgnEffError[i]/tgraphSgnEff[i])+(tgraphSgnEffError145[i]/tgraphSgnEff145[i])*(tgraphSgnEffError145[i]/tgraphSgnEff145[i]));
	  tgraphSgnEffError[i]*=tgraphSgnEff[i];
	  
  }
  
  TGraphErrors* tg4= new TGraphErrors(masspointConst,tgraphMass,tgraphSigEff,tgraphMassError,tgraphSigEffError);
  TGraphErrors* tg5= new TGraphErrors(masspointConst,tgraphMass,tgraphBkgEff,tgraphMassError,tgraphBkgEffError);
  TGraphErrors* tg6= new TGraphErrors(masspointConst,tgraphMass,tgraphSgnEff,tgraphMassError,tgraphSgnEffError);
  
  
  tg1->SetMaximum(1);
  tg1->SetMinimum(0);
  tg1->SetTitle("");
  tg1->GetXaxis()->SetTitle("Mass of Z' [GeV]");
  tg1->Draw("APL*");
  tg2->SetLineColor(2);
  tg3->SetLineColor(3);
  tg1->SetMarkerStyle(20);
  tg2->SetMarkerStyle(21);
  tg3->SetMarkerStyle(22);
  leg->Clear();
  tg1->SetFillColor(kWhite);
  tg2->SetFillColor(kWhite);
  tg3->SetFillColor(kWhite);
  gStyle->SetFillColor(kWhite) ;
  leg->AddEntry(tg1,"signal eff.");
  leg->AddEntry(tg2,"bkg. eff.");
  leg->AddEntry(tg3,"significance");
  tg2->Draw("PL,same");
  tg3->Draw("PL,same");
  leg->Draw("same");
  c2->Print("png/1.png");
  
  
  tg4->GetXaxis()->SetTitle("Mass of Z' [GeV]");
  tg4->SetTitle("most significant window/105-145");
  tg4->SetMaximum(2.2);
  tg4->SetMinimum(0.3);
  tg4->SetMarkerStyle(20);
  tg5->SetMarkerStyle(21);
  tg6->SetMarkerStyle(22);
  tg4->Draw("APL");
  tg5->SetLineColor(2);
  tg6->SetLineColor(3);
  
  leg->Clear();
  tg4->SetFillColor(kWhite);
  tg5->SetFillColor(kWhite);
  tg6->SetFillColor(kWhite);
  gStyle->SetFillColor(kWhite) ;
  leg->AddEntry(tg4,"signal eff. ratio");
  leg->AddEntry(tg5,"bkg. eff. ratio");
  leg->AddEntry(tg6,"significance ratio");
  tg5->Draw("PL,same");
  tg6->Draw("PL,same");
  leg->Draw("same");
  c2->Print("png/2.png");
   
   
   
   
   
   
  ofstream myfile;
  myfile.open ("txt/twikiOPPRN.txt");
  myfile<<"|*windowRange*|";
  for (int massP=0;massP<4;massP++)myfile<<"*"<<masspoint[massP].data()<<"Num*|";
  //myfile<<"*";
  for (int massP=0;massP<4;massP++)myfile<<"*"<<masspoint[massP].data()<<"Eff*|";
  myfile<<"*avg.rank(eff)*|*avg.rank(Num)*|*avg.rank(total)*|"<<endl;
  for(int j=0;j<nBmin;j++){
	for(int i=0;i<nWidth;i++){
		if(width[i]+bmin[j]>150)continue;
		myfile<<"|"<<bmin[j]<<"to"<<bmin[j]+width[i]<<"|";
	    double temp1=0,temp2=0;
		for (int massP=0;massP<4;massP++){
			myfile<<twikiSignNum[massP][i][j]<<"("<<twikiSignNumValue[massP][i][j]<<")|";
			temp1+=twikiSign[massP][i][j];
			temp2+=twikiSignNum[massP][i][j];
		}
		for (int massP=0;massP<4;massP++){
			myfile<<twikiSign[massP][i][j]<<"("<<twikiSignValue[massP][i][j]<<")|";
		}
		myfile<<temp1/4<<"|"<<temp2/4<<"|"<<(temp1+temp2)/8<<"|"<<endl;
		//<<"to"<<bmin[i]+width[j]<<"|"<<endl;
		//cout<<j<<","<<i<<endl;
	}
  }
  myfile<<endl;
  
  
  
  myfile<<"|*rank(signalEff)*|";
  for (int massP=0;massP<4;massP++)myfile<<"*"<<masspoint[massP].data()<<"Eff*|";
  myfile<<endl;
  vector<int> signINum,signJNum;
 //for(int ij=0;ij<nWidth*nBmin;ij++){
    for(int ij=0;ij<15;ij++){
    myfile<<"|"<<ij+1<<"|";
	for (int massP=0;massP<4;massP++){
		double tempSign=1000;
    	int tempI=0,tempJ=0;
		//for(int i=0;i<nWidth;i++)for(int j=0;j<nBmin;j++)if(width[tempI]+bmin[tempJ]>150)twikiSign[massP][tempI][tempJ]=10000;
		
		for(int i=0;i<nWidth;i++){
			for(int j=0;j<nBmin;j++){
				if(twikiSign[massP][i][j]<tempSign){
					if(width[i]+bmin[j]>150)continue;
					//cout<<"i="<<i<<",j="<<j<<",signcp="<<signCP[i][j]<<",temp="<<tempSign<<endl;
					tempSign=twikiSign[massP][i][j];
					tempI=i;
					tempJ=j;
					//cout<<tempI<<","<<tempJ<<","<<tempSign<<endl;
				}
			}
		}
		//signINum.push_back(tempI);
		//signJNum.push_back(tempJ);
		//cout<<signCP[tempI][tempJ]<<endl;
	//cout<<"twikiSign[massP][tempI][tempJ]="<<twikiSign[massP][tempI][tempJ]<<"|"<<bmin[tempJ]<<"to"<<bmin[tempJ]+width[tempI]<<"|"<<endl;
		twikiSign[massP][tempI][tempJ]=10000;
		//if(width[tempI]+bmin[tempJ]>150)continue;
		myfile<<bmin[tempJ]<<"to"<<bmin[tempJ]+width[tempI]<<"("<<setprecision(3)<<twikiEffValue[massP][tempI][tempJ]<<")|";
		//cout<<twikiSign[massP][tempI][tempJ];
		
	}
	myfile<<endl;
  }
  myfile<<endl;
  
  myfile<<"|*rank(signalEff)*|";
  for (int massP=0;massP<4;massP++)myfile<<"*"<<masspoint[massP].data()<<"Num*|";
  myfile<<endl;
  //vector<int> signINum,signJNum;
 //for(int ij=0;ij<nWidth*nBmin;ij++){
    for(int ij=0;ij<15;ij++){
    myfile<<"|"<<ij+1<<"|";
	for (int massP=0;massP<4;massP++){
		double tempSign=10000;
    	int tempI=0,tempJ=0;
		//for(int i=0;i<nWidth;i++)for(int j=0;j<nBmin;j++)if(width[tempI]+bmin[tempJ]>150)twikiSignNum[massP][tempI][tempJ]=10000;
		
		for(int i=0;i<nWidth;i++){
			for(int j=0;j<nBmin;j++){
				if(twikiSignNum[massP][i][j]<tempSign){
					if(width[i]+bmin[j]>150)continue;
					//cout<<"i="<<i<<",j="<<j<<",signcp="<<signCP[i][j]<<",temp="<<tempSign<<endl;
					tempSign=twikiSignNum[massP][i][j];
					tempI=i;
					tempJ=j;
					//cout<<tempI<<","<<tempJ<<","<<tempSign<<endl;
				}
			}
		}
		//signINum.push_back(tempI);
		//signJNum.push_back(tempJ);
		//cout<<signCP[tempI][tempJ]<<endl;
		twikiSignNum[massP][tempI][tempJ]=10000;
		//if(width[tempI]+bmin[tempJ]>150)continue;
		myfile<<bmin[tempJ]<<"to"<<bmin[tempJ]+width[tempI]<<"("<<setprecision(3)<<twikiEffValue[massP][tempI][tempJ]<<")|";
		//cout<<twikiSign[massP][tempI][tempJ];
		
	}
	myfile<<endl;
  }
  myfile<<endl;
}

double  getEff(TH1F *h,int bmin,int bmax){
	if(h->Integral()==0){
		cout<<h->Integral();
		return 0;
	}
  double denom=h->Integral();//cout<<"de="<<denom;
  double nomin=h->Integral(bmin,bmax);//cout<<"no="<<nomin<<endl;

  
  return nomin/denom;
}

double  getErr(TH1F *h,int bmin,int bmax){
	if(h->Integral()==0)return 0;
  double temp=0,temp2=0;
  double denom=h->Integral();
  //double nomin=h->Integral(bmin,bmax);
  for(int i=1;i<h->GetNbinsX();i++){
		if(  i>bmin && i<bmax+1){
			temp+=h->GetBinError(i)*h->GetBinError(i);
		}
		else 
			temp2+=h->GetBinError(i)*h->GetBinError(i);
	}
  
  temp=temp/(h->Integral()*h->Integral())-2*h->Integral(bmin,bmax)*temp/(h->Integral()*h->Integral()*h->Integral())+h->Integral(bmin,bmax)*h->Integral(bmin,bmax)*(temp+temp2)/(h->Integral()*h->Integral()*h->Integral()*h->Integral());
  
  return sqrt(temp);
}

double  getErr2(TH1F *h,int bmin,int bmax){
  double temp=0;
  double denom=h->Integral();
  //double nomin=h->Integral(bmin,bmax);
  for(int i=bmin;i<bmax+1;i++){
	  temp+=h->GetBinError(i)*h->GetBinError(i);
  }
  temp=sqrt(temp);
  return temp;
}

double getSign(double d1,double d2){
  return d1/(1+sqrt(d2));
}