// example code to run Bulk Graviton->ZZ->ZlepZhad selections on electron-channel

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TString.h>
#include <map>
#include <TH1D.h>
#include <TFile.h>
#include "../../../untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TF1.h>

using namespace std;


float getPUPPIweight(float puppipt, float puppieta ){

	double barrleMean[13]={1.16648,1.15105,1.14657,1.14295,1.13818,1.13918,1.14078,1.1424,1.14941,1.15595,1.16176,1.16176,};
	double endcapMean[13]={1.32391,1.29096,1.27103,1.2598,1.24856,1.22939,1.21945,1.22035,1.226,1.22177,1.23195,1.23195,};
	double barrlex[13]={346.222,379.026,420.738,466.036,556.089,645.072,733.387,819.748,903.889,1119.54,1334.59,5000,};
	double endcapx[13]={328.682,354.991,391.634,430.037,499.404,564.743,623.403,674.7,725.021,848.722,953.155,5000,};
	TGraph* tg1=new TGraph(12,barrlex,barrleMean);
	TGraph* tg2=new TGraph(12,endcapx,endcapMean);

	float totalWeight=1;
	if( fabs(puppieta)  <= 1.3 )return tg1->Eval(puppipt);
	else return tg2->Eval(puppipt);
}

float getPUPPIweight_o(float puppipt, float puppieta ){

   TF1* puppisd_corrGEN      = new TF1("puppisd_corrGEN","[0]+[1]*pow(x*[2],-[3])");
  puppisd_corrGEN->SetParameters(1.00626, -1.06161, 0.07999,1.20454 );
  TF1* puppisd_corrRECO_cen = new TF1("puppisd_corrRECO_cen","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)");
  puppisd_corrRECO_cen->SetParameters( 1.05807,-5.91971e-05,2.296e-07, -1.98795e-10,6.67382e-14,  -7.80604e-18);
  TF1* puppisd_corrRECO_for = new TF1("puppisd_corrRECO_for","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)");
  puppisd_corrRECO_for->SetParameters( 1.26638,-0.000658496,  9.73779e-07,-5.93843e-10, 1.61619e-13, -1.6272e-17);
  float genCorr  = 1.;
  float recoCorr = 1.;
  float totalWeight = 1.;
  genCorr =  puppisd_corrGEN->Eval( puppipt );
  if( fabs(puppieta)  <= 1.3 ) recoCorr = puppisd_corrRECO_cen->Eval( puppipt );
  else if( fabs(puppieta) > 1.3 ) recoCorr = puppisd_corrRECO_for->Eval( puppipt );
  totalWeight = genCorr * recoCorr;
  return totalWeight;
}


void xAna_hh_massResolutionBase(std::string inputFile,TString outputFile,bool matchb=false, bool debug=false, bool cut=false){


  cout << "output file name = " << outputFile.Data() << endl;

  //get TTree from file ...
  TreeReader data(inputFile.data());
  
  Long64_t nTotal=0;
  Long64_t nPass[20]={0};

  const int nHistos=3;

  TH1F* h_massDiff = new TH1F("h_massDiff","",100,-0.5,0.5);
  //TH1F* h_mass     = new TH1F("h_mass","",100,0,200);
  TH1F* h_mass     = new TH1F("h_mass","",100,62.5,187.5);

  TH1F* h_AK8SD[nHistos];
  TH1F* h_AK8SDCorrThea[nHistos];
  TH1F* h_AK8SDHCorr[nHistos];
  TH1F* h_PR[nHistos];
  TH1F* h_PRCorr[nHistos];

  TH1F* h_diff_PR[nHistos];
  TH1F* h_diff_PRCorr[nHistos];
  TH1F* h_diff_AK8SD[nHistos];
  TH1F* h_diff_AK8SDCorrThea[nHistos];
  TH1F* h_diff_AK8SDHCorr[nHistos];
  
  std::string prefix[]={"leading","subleading","both"};
  for(int i=0; i<nHistos; i++)
    {

	h_AK8SD[i] = (TH1F*)h_mass->Clone(Form("h_AK8SD_%s",prefix[i].data()));
      h_AK8SD[i]->SetXTitle("AK8 Raw Puppi+Softdrop mass [GeV]");

	h_AK8SDCorrThea[i] = (TH1F*)h_mass->Clone(Form("h_AK8SDCorrThea_%s",prefix[i].data()));
      h_AK8SDCorrThea[i]->SetXTitle("AK8 Thea-corrected Puppi+Softdrop mass [GeV]");
	
	h_AK8SDHCorr[i] = (TH1F*)h_mass->Clone(Form("h_AK8SDHCorr_%s",prefix[i].data()));
      h_AK8SDHCorr[i]->SetXTitle("AK8 H-corrected Puppi+Softdrop mass [GeV]");

      h_PR[i] = (TH1F*)h_mass->Clone(Form("h_PR_%s",prefix[i].data()));
      h_PR[i]->SetXTitle("Raw CHS+Pruned mass [GeV]");

      h_PRCorr[i] = (TH1F*)h_mass->Clone(Form("h_PRCorr_%s",prefix[i].data()));
      h_PRCorr[i]->SetXTitle("L2L3-corrected CHS+Pruned mass [GeV]");

      // study the difference with respect to 125 GeV

     
	
	h_diff_AK8SD[i] = (TH1F*)h_massDiff->Clone(Form("h_diff_AK8SD_%s",prefix[i].data()));
      h_diff_AK8SD[i]->SetXTitle("AK8 Raw Puppi+Softdrop (m-125)/125");

	h_diff_AK8SDCorrThea[i] = (TH1F*)h_massDiff->Clone(Form("h_diff_AK8SDCorrThea_%s",prefix[i].data()));
      h_diff_AK8SDCorrThea[i]->SetXTitle("AK8 Thea-corrected Puppi+Softdrop (m-125)/125");
	
	h_diff_AK8SDHCorr[i] = (TH1F*)h_massDiff->Clone(Form("h_diff_AK8SDHCorr_%s",prefix[i].data()));
      h_diff_AK8SDHCorr[i]->SetXTitle("AK8 H-corrected Puppi+Softdrop (m-125)/125");

      h_diff_PR[i] = (TH1F*)h_massDiff->Clone(Form("h_diff_PR_%s",prefix[i].data()));
      h_diff_PR[i]->SetXTitle("Raw CHS+Pruned (m-125)/125");

      h_diff_PRCorr[i] = (TH1F*)h_massDiff->Clone(Form("h_diff_PRCorr_%s",prefix[i].data()));
      h_diff_PRCorr[i]->SetXTitle("L2L3-corrected CHS+Pruned (m-125)/125");

    }


  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){

  if((jEntry+1)%2)continue;
    if (jEntry % 1000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
    if(debug && jEntry>10)break;

    data.GetEntry(jEntry);
    nTotal++;

    Int_t nGenPar        = data.GetInt("nGenPar");
    Int_t* genParId      = data.GetPtrInt("genParId");
    Int_t* genParSt      = data.GetPtrInt("genParSt");
    Int_t* genMomParId   = data.GetPtrInt("genMomParId");
    Int_t* genDa1      = data.GetPtrInt("genDa1");
    Int_t* genDa2      = data.GetPtrInt("genDa2");

    int genHIndex[2]={-1,-1};
    int genbIndex[2][2]={{-1,-1},
			 {-1,-1}};		       

    for(int ig=0; ig < nGenPar; ig++){

      if(genParId[ig]!=25)continue;

      if(genHIndex[0]<0)
	{
	  genHIndex[0]=ig;
	  genbIndex[0][0]=genDa1[ig];
	  genbIndex[0][1]=genDa2[ig];
	}

      else if(genHIndex[1]<0)
	{
	  genHIndex[1]=ig;
	  genbIndex[1][0]=genDa1[ig];
	  genbIndex[1][1]=genDa2[ig];
	}

    }    

    if(genHIndex[0]<0 || genHIndex[1]<0)continue;
    if(genbIndex[0][0]<0 || genbIndex[0][1]<0)continue;
    if(genbIndex[1][0]<0 || genbIndex[1][1]<0)continue;

    nPass[0]++;

    if(genHIndex[0]==genHIndex[1])continue;
    nPass[1]++;

    TLorentzVector genH_l4[2];
    TLorentzVector genb_l4[2][2];
    TClonesArray* genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");

    for(int ih=0; ih<2; ih++)
      {
	genH_l4[ih] = *((TLorentzVector*)genParP4->At(genHIndex[ih]));
	for(int ib=0; ib<2; ib++)
	  {
	    genb_l4[ih][ib] = *((TLorentzVector*)genParP4->At(genbIndex[ih][ib]));
	  }
      }


   
    int nFATJet         = data.GetInt("FATnJet");
    const int nJets=nFATJet;
    TClonesArray* FATjetP4   = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    TClonesArray* AK8PuppijetP4 = (TClonesArray*) data.GetPtrTObject("AK8PuppijetP4");

    // check matching first

    bool findAMatch=false;
    const float dRMax=0.4;
    const float dRbMax=0.8;
    int matchedHJetIndex[2]={-1,-1};
		      
    for(int ij=0; ij<nJets; ij++)
      {
	TLorentzVector* thisJet = (TLorentzVector*)FATjetP4->At(ij);

	for(int jj=0; jj<nJets; jj++)
	  {

	    if(ij==jj)continue;
	    TLorentzVector* thatJet = (TLorentzVector*)FATjetP4->At(jj);
	    
	    if(thisJet->DeltaR(genH_l4[0])<dRMax && 
	       (!matchb || (matchb && 
			    thisJet->DeltaR(genb_l4[0][0])<dRbMax && 
			    thisJet->DeltaR(genb_l4[0][1])<dRbMax)) &&
	       thatJet->DeltaR(genH_l4[1])<dRMax &&
	       (!matchb || (matchb &&
			    thatJet->DeltaR(genb_l4[1][0])<dRbMax &&
			    thatJet->DeltaR(genb_l4[1][1])<dRbMax)))
	      {
		if(debug)
		  {
		    cout << "dRhb00= " <<  thisJet->DeltaR(genb_l4[0][0]) << endl;
		    cout << "dRhb01= " <<  thisJet->DeltaR(genb_l4[0][1]) << endl;
		    cout << "dRhb10= " <<  thatJet->DeltaR(genb_l4[1][0]) << endl;
		    cout << "dRhb11= " <<  thatJet->DeltaR(genb_l4[1][1]) << endl;
		  }
		if(ij<jj){
		  matchedHJetIndex[0]=ij;
		  matchedHJetIndex[1]=jj;
		}
		else
		  {
		    matchedHJetIndex[0]=jj;
		    matchedHJetIndex[1]=ij;
		  }
		findAMatch=true;
		break;
	      }

	    if(findAMatch)break;

	  }	

	if(findAMatch)break;

      }

    if(!findAMatch)continue;
    if(debug)
      cout << matchedHJetIndex[0] << "\t" << matchedHJetIndex[1] << endl;
    nPass[2]++;

    //0. has a good vertex
    Int_t nVtx        = data.GetInt("nVtx");
    if(nVtx<1 && cut)continue;
    nPass[3]++;



    Float_t*  fatjetTau1 = data.GetPtrFloat("FATjetTau1");
    Float_t*  fatjetTau2 = data.GetPtrFloat("FATjetTau2");
    Float_t*  fatjetCISVV2 = data.GetPtrFloat("FATjetCISVV2");
    Float_t*  fatjetPRmass = data.GetPtrFloat("FATjetPRmass");
    Float_t*  fatjetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    Float_t*  fatjetSDmass = data.GetPtrFloat("FATjetPuppiSDmass");
    Float_t*  fatjetSDmassL2L3Corr = data.GetPtrFloat("FATjetPuppiSDmassL2L3Corr"); 
    Float_t*  AK8PuppijetSDmass = data.GetPtrFloat("AK8PuppijetSDmass");

    vector<bool>    &passFatJetLooseID = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
    
    TLorentzVector recoH_l4[2];
    int nGoodJets=0;
    
    
    
    for(int i=0; i<2; i++)
      {
    	
	int ij = matchedHJetIndex[i];
     	TLorentzVector* thisJet = (TLorentzVector*)FATjetP4->At(ij);
	recoH_l4[i]= (*thisJet);
    	if(thisJet->Pt()<200)continue;
	if(fabs(thisJet->Eta())>2.4)continue;
	nGoodJets++;
      }

    if(nGoodJets<2)continue;
    nPass[5]++;
    
    

    float dEta=fabs(recoH_l4[0].Eta()-recoH_l4[1].Eta());
    if(dEta>1.3 && cut)continue;
    nPass[6]++;

    float M=(recoH_l4[0] + recoH_l4[1]).M();
    if(M<800 && cut)continue;
    nPass[7]++;


    int nHP=0;
    int nLP=0;
    for(int i=0; i<2; i++)
      {
    	
	int ij = matchedHJetIndex[i];

	float tau21_i = fatjetTau2[ij]/fatjetTau1[ij];
	bool isHP= (tau21_i < 0.6);
	if(isHP)nHP++;
      }

    if(nHP<2 && cut)continue;
    nPass[8]++;

	double SDMass[2];
	Int_t AK8PuppinJet        = data.GetInt("AK8PuppinJet");
	if(AK8PuppinJet<2)continue;
	bool matchThis=0,matchThat=0;
	//TClonesArray* AK8PuppijetP4 = (TClonesArray*) data.GetPtrTObject("AK8PuppijetP4");
	
	TLorentzVector* thisJet,* thatJet;
			thisJet= (TLorentzVector*)FATjetP4->At(0);
			thatJet= (TLorentzVector*)FATjetP4->At(1);
			for(int i=0;i<AK8PuppinJet;i++){
				TLorentzVector* thisAddJet ;
				thisAddJet= (TLorentzVector*)AK8PuppijetP4->At(i);
				if(!matchThis && thisAddJet->DeltaR(*thisJet)<0.8){
					matchThis=1;
					SDMass[0]=AK8PuppijetSDmass[i];
					continue;
				}
				if(!matchThat && thisAddJet->DeltaR(*thatJet)<0.8){
					matchThat=1;
					SDMass[1]=AK8PuppijetSDmass[i];
				}
				if(matchThis&& matchThat)break;
			}
    
    
    
    
    for(int i=0; i<2;i++)
      {
	
	float thea_corr = getPUPPIweight_o(thisJet->Pt(),thisJet->Eta());
	if(i==1)thea_corr = getPUPPIweight_o(thatJet->Pt(),thatJet->Eta());
	float H_corr= getPUPPIweight(thisJet->Pt(),thisJet->Eta());
	if(i==1)H_corr = getPUPPIweight(thatJet->Pt(),thatJet->Eta());
	
	h_PR[i]->Fill(fatjetPRmass[i]);
	h_PRCorr[i]->Fill(fatjetPRmassL2L3Corr[i]);
	h_AK8SD[i]->Fill(SDMass[i]);
	h_AK8SDCorrThea[i]->Fill(SDMass[i]*thea_corr);
	h_AK8SDHCorr[i]->Fill(SDMass[i]*H_corr);

	h_PR[2]->Fill(fatjetPRmass[i]);
	h_PRCorr[2]->Fill(fatjetPRmassL2L3Corr[i]);
	h_AK8SD[2]->Fill(SDMass[i]);
	h_AK8SDCorrThea[2]->Fill(SDMass[i]*thea_corr);
	h_AK8SDHCorr[2]->Fill(SDMass[i]*H_corr);

	h_diff_PR[i]->Fill((fatjetPRmass[i]-125)/125);
	h_diff_PRCorr[i]->Fill((fatjetPRmassL2L3Corr[i]-125)/125);
	h_diff_AK8SD[i]->Fill((SDMass[i]-125)/125);
	h_diff_AK8SDCorrThea[i]->Fill((SDMass[i]*thea_corr-125)/125);
	h_diff_AK8SDHCorr[i]->Fill((SDMass[i]*H_corr-125)/125);

	h_diff_PR[2]->Fill((fatjetPRmass[i]-125)/125);
	h_diff_PRCorr[2]->Fill((fatjetPRmassL2L3Corr[i]-125)/125);
	h_diff_AK8SD[2]->Fill((SDMass[i]-125)/125);
	h_diff_AK8SDCorrThea[2]->Fill((SDMass[i]*thea_corr-125)/125);
	h_diff_AK8SDHCorr[2]->Fill((SDMass[i]*H_corr-125)/125);
      }
    

  } // end of loop over entries

  std::cout << "nTotal    = " << nTotal << std::endl;
  for(int i=0;i<20;i++)
    if(nPass[i]>0)
      std::cout << "nPass[" << i << "]= " << nPass[i] << std::endl;

  TFile* outFile = new TFile(Form("output/%s",outputFile.Data()),"recreate");

  for(int i=0; i<nHistos; i++)
    {
      h_diff_PR[i]->Write();
      h_diff_PRCorr[i]->Write();
	h_diff_AK8SD[i]->Write();
	h_diff_AK8SDCorrThea[i]->Write();
	h_diff_AK8SDHCorr[i]->Write();
      
      h_PR[i]->Write();
      h_PRCorr[i]->Write();
	h_AK8SD[i]->Write();
	h_AK8SDCorrThea[i]->Write();
	h_AK8SDHCorr[i]->Write();
    }
  outFile->Close();

}

void xAna_hh_massResolution(){
	
	xAna_hh_massResolutionBase("/data7/chchen/AK8PuppijetGenSDmass/B/B1000.root","B1000.root");
	xAna_hh_massResolutionBase("/data7/chchen/AK8PuppijetGenSDmass/B/B1200.root","B1200.root");
	xAna_hh_massResolutionBase("/data7/chchen/AK8PuppijetGenSDmass/B/B1400.root","B1400.root");
	xAna_hh_massResolutionBase("/data7/chchen/AK8PuppijetGenSDmass/B/B1600.root","B1600.root");
	xAna_hh_massResolutionBase("/data7/chchen/AK8PuppijetGenSDmass/B/B1800.root","B1800.root");
	xAna_hh_massResolutionBase("/data7/chchen/AK8PuppijetGenSDmass/B/B2000.root","B2000.root");
	xAna_hh_massResolutionBase("/data7/chchen/AK8PuppijetGenSDmass/B/B2500.root","B2500.root");
	xAna_hh_massResolutionBase("/data7/chchen/AK8PuppijetGenSDmass/B/B3000.root","B3000.root");
	xAna_hh_massResolutionBase("/data7/chchen/AK8PuppijetGenSDmass/B/B4000.root","B4000.root");
	xAna_hh_massResolutionBase("/data7/chchen/AK8PuppijetGenSDmass/B/B4500.root","B4500.root");
	//xAna_hh_massResolutionBase("/data7/chchen/8_0_20/B/B1200.root",PU_up,PU_down);
}
/*
void xAna_hh_massResolution(){
	xAna_hh_massResolutionPt(200,300);
	xAna_hh_massResolutionPt(300,400);
	xAna_hh_massResolutionPt(400,500);
	xAna_hh_massResolutionPt(500,600);
	xAna_hh_massResolutionPt(600,700);
	xAna_hh_massResolutionPt(700,800);
	xAna_hh_massResolutionPt(800,900);
	xAna_hh_massResolutionPt(900,1000);
	xAna_hh_massResolutionPt(1250,1500);
	xAna_hh_massResolutionPt(1500,5000);
	
}
*/