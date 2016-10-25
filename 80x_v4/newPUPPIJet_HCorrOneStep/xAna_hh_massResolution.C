// example code to run Bulk Graviton->ZZ->ZlepZhad selections on electron-channel

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TString.h>
#include <map>
#include <TH1D.h>
#include <TFile.h>
#include "../untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TF1.h>

using namespace std;


float getPUPPIweight(float puppipt, float puppieta ){

  // reading the root file takes too long
  // TFile* file = TFile::Open( "weights/puppiCorr.root","READ");
  // TF1* puppisd_corrGEN      = (TF1*)file->Get("puppiJECcorr_gen");
  // TF1* puppisd_corrRECO_cen = (TF1*)file->Get("puppiJECcorr_reco_0eta1v3");
  // TF1* puppisd_corrRECO_for = (TF1*)file->Get("puppiJECcorr_reco_1v3eta2v5");

  TF1* puppisd_corrGEN      = new TF1("puppisd_corrGEN","[0]+[1]*pow(x*[2],-[3])");
  puppisd_corrGEN->SetParameters(
   				1.08383e+00,
   				1.80932e+01  ,
				1.13231e-02,
   				 4.05749e+00
   				 );
TF1* puppisd_corrGENEndcap      = new TF1("puppisd_corrGEN","[0]+[1]*pow(x*[2],-[3])");
  puppisd_corrGENEndcap->SetParameters(
   				1.06913e+00,
   				1.06485e+01,
   				1.50484e-02,
   				3.62054e+00
   				 );
  TF1* puppisd_corrRECO_cen = new TF1("puppisd_corrRECO_cen","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)");
  puppisd_corrRECO_cen->SetParameters(
   				      1.63082e+00,
   				    -1.93676e-03,
   				      2.59101e-06,
   				   -1.58447e-09,
   				     4.49982e-13,
   				-4.85023e-17
   				      );

  TF1* puppisd_corrRECO_for = new TF1("puppisd_corrRECO_for","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)");
  puppisd_corrRECO_for->SetParameters(
   				    2.07578e+00,
   				-3.92989e-03,
   				   6.64560e-06,
   				 -5.69331e-09,
   				   2.44060e-12,
   				   -4.13464e-16);
				   
TF1* oneBarel = new TF1("oneBarel","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)");
  oneBarel->SetParameters(/*
   				      2.32790e+00,
   				    -4.28162e-03,
   				     5.88872e-06,
   				   -3.69312e-09,
   				    1.07407e-12,
   				-1.18250e-16*/
				2.26084e+00,
				-4.03837e-03,
				5.54725e-06,
				-3.46742e-09,
				1.00377e-12,
				 -1.09930e-16
   				      );

  TF1* oneEndcap = new TF1("oneEndcap","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)");
  oneEndcap->SetParameters(
  /*
   				2.85414e+00,
   				-7.67771e-03,
   				 1.43612e-05,
   				-1.33354e-08,
   			 6.11380e-12,
   				-1.10337e-15*/
			2.78161e+00,
			-7.32847e-03,
			1.37222e-05,
			-1.27878e-08,
			5.89519e-12,
			-1.07097e-15
				);

  float genCorr  = 1.;
  float recoCorr = 1.;
  float totalWeight = 1.;
/*
  if( fabs(puppieta)  <= 1.3 )
  genCorr =  puppisd_corrGEN->Eval( puppipt );
else  genCorr =  puppisd_corrGENEndcap->Eval( puppipt );
  if( fabs(puppieta)  <= 1.3 ) recoCorr = puppisd_corrRECO_cen->Eval( puppipt );
  else
    if( fabs(puppieta) > 1.3 ) recoCorr = puppisd_corrRECO_for->Eval( puppipt );

  totalWeight = genCorr * recoCorr;
  */
   if( fabs(puppieta)  <= 1.3 )totalWeight =  oneBarel->Eval( puppipt );
else  totalWeight =  oneEndcap->Eval( puppipt );
  // file->Close();
  return totalWeight;
}


void xAna_hh_massResolutionBase(std::string inputFile,TString outputFile, bool matchb=false, bool debug=false, bool cut=false){


  cout << "output file name = " << outputFile.Data() << endl;


  //get TTree from file ...
  TreeReader data(inputFile.data());
  
  Long64_t nTotal=0;
  Long64_t nPass[20]={0};

  const int nHistos=3;

  TH1F* h_massDiff = new TH1F("h_massDiff","",100,-0.5,0.5);
  //TH1F* h_mass     = new TH1F("h_mass","",100,0,200);
  TH1F* h_mass     = new TH1F("h_mass","",100,62.5,187.5);

  TH1F* h_SD[nHistos];
  TH1F* h_SDCorr[nHistos];
  TH1F* h_SDCorrThea[nHistos];
  
  TH1F* h_AK8SD[nHistos];
  TH1F* h_AK8SDCorrThea[nHistos];
  
  TH1F* h_PR[nHistos];
  TH1F* h_PRCorr[nHistos];

  TH1F* h_diff_SD[nHistos];
  TH1F* h_diff_SDCorr[nHistos];
  TH1F* h_diff_SDCorrThea[nHistos];
  TH1F* h_diff_PR[nHistos];
  TH1F* h_diff_PRCorr[nHistos];
  
  TH1F* h_diff_AK8SD[nHistos];
  TH1F* h_diff_AK8SDCorrThea[nHistos];

  std::string prefix[]={"leading","subleading","both"};
  for(int i=0; i<nHistos; i++)
    {

      h_SD[i] = (TH1F*)h_mass->Clone(Form("h_SD_%s",prefix[i].data()));
      h_SD[i]->SetXTitle("Raw Puppi+Softdrop mass [GeV]");
	
	h_AK8SD[i] = (TH1F*)h_mass->Clone(Form("h_AK8SD_%s",prefix[i].data()));
      h_AK8SD[i]->SetXTitle("AK8 Raw Puppi+Softdrop mass [GeV]");

      h_SDCorr[i] = (TH1F*)h_mass->Clone(Form("h_SDCorr_%s",prefix[i].data()));
      h_SDCorr[i]->SetXTitle("L2L3-corrected Puppi+Softdrop mass [GeV]");

      h_SDCorrThea[i] = (TH1F*)h_mass->Clone(Form("h_SDCorrThea_%s",prefix[i].data()));
      h_SDCorrThea[i]->SetXTitle("Thea-corrected Puppi+Softdrop mass [GeV]");
	
	h_AK8SDCorrThea[i] = (TH1F*)h_mass->Clone(Form("h_AK8SDCorrThea_%s",prefix[i].data()));
      h_AK8SDCorrThea[i]->SetXTitle("AK8 Thea-corrected Puppi+Softdrop mass [GeV]");

      h_PR[i] = (TH1F*)h_mass->Clone(Form("h_PR_%s",prefix[i].data()));
      h_PR[i]->SetXTitle("Raw CHS+Pruned mass [GeV]");

      h_PRCorr[i] = (TH1F*)h_mass->Clone(Form("h_PRCorr_%s",prefix[i].data()));
      h_PRCorr[i]->SetXTitle("L2L3-corrected CHS+Pruned mass [GeV]");

      // study the difference with respect to 125 GeV

      h_diff_SD[i] = (TH1F*)h_massDiff->Clone(Form("h_diff_SD_%s",prefix[i].data()));
      h_diff_SD[i]->SetXTitle("Raw Puppi+Softdrop (m-125)/125");
	
	h_diff_AK8SD[i] = (TH1F*)h_massDiff->Clone(Form("h_diff_AK8SD_%s",prefix[i].data()));
      h_diff_AK8SD[i]->SetXTitle("AK8 Raw Puppi+Softdrop (m-125)/125");

      h_diff_SDCorr[i] = (TH1F*)h_massDiff->Clone(Form("h_diff_SDCorr_%s",prefix[i].data()));
      h_diff_SDCorr[i]->SetXTitle("L2L3-corrected Puppi+Softdrop (m-125)/125");

      h_diff_SDCorrThea[i] = (TH1F*)h_massDiff->Clone(Form("h_diff_SDCorrThea_%s",prefix[i].data()));
      h_diff_SDCorrThea[i]->SetXTitle("Thea-corrected Puppi+Softdrop (m-125)/125");
	
	h_diff_AK8SDCorrThea[i] = (TH1F*)h_massDiff->Clone(Form("h_diff_AK8SDCorrThea_%s",prefix[i].data()));
      h_diff_AK8SDCorrThea[i]->SetXTitle("AK8 Thea-corrected Puppi+Softdrop (m-125)/125");

      h_diff_PR[i] = (TH1F*)h_massDiff->Clone(Form("h_diff_PR_%s",prefix[i].data()));
      h_diff_PR[i]->SetXTitle("Raw CHS+Pruned (m-125)/125");

      h_diff_PRCorr[i] = (TH1F*)h_massDiff->Clone(Form("h_diff_PRCorr_%s",prefix[i].data()));
      h_diff_PRCorr[i]->SetXTitle("L2L3-corrected CHS+Pruned (m-125)/125");

    }


  for(Long64_t jEntry=25000; jEntry<data.GetEntriesFast() ;jEntry++){

    if (jEntry % 1000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());


    if(debug && jEntry>10)break;

    data.GetEntry(jEntry);
    nTotal++;

    //2. pass electron or muon trigger                                                                                                                                          
    std::string* trigName = data.GetPtrString("hlt_trigName");
    vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
    const Int_t nsize = data.GetPtrStringSize();

    bool passTrigger=false;
    for(int it=0; it< nsize; it++)
      {
	std::string thisTrig= trigName[it];
        bool results = trigResult[it];

        // std::cout << thisTrig << " : " << results << std::endl;                                                                                                              

        if( (thisTrig.find("HLT_PFHT900_v")!= std::string::npos && results==1)
            )
          {
            passTrigger=true;
            break;
          }


      }

    if(!passTrigger && cut)continue;
    nPass[4]++;



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


    if(debug){
      cout << genHIndex[0] << "\t" << genHIndex[1] << endl;
      genH_l4[0].Print();
      genH_l4[1].Print();

      cout << genbIndex[0][0] << "\t" << genbIndex[0][1] << "\t" 
	   << genbIndex[1][0] << "\t" << genbIndex[1][1] << endl;
      genH_l4[0].Print();
      genH_l4[1].Print();
      genb_l4[0][0].Print();
      genb_l4[0][1].Print();
      genb_l4[1][0].Print();
      genb_l4[1][1].Print();

    }
    int nFATJet         = data.GetInt("FATnJet");
    const int nJets=nFATJet;
    TClonesArray* fatjetP4   = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    TClonesArray* puppijetP4 = (TClonesArray*) data.GetPtrTObject("FATjetPuppiP4");
    TClonesArray* AK8PuppijetP4 = (TClonesArray*) data.GetPtrTObject("AK8PuppijetP4");

    // check matching first

    bool findAMatch=false;
    const float dRMax=0.4;
    const float dRbMax=0.8;
    int matchedHJetIndex[2]={-1,-1};
		      
    for(int ij=0; ij<nJets; ij++)
      {
	TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);

	for(int jj=0; jj<nJets; jj++)
	  {

	    if(ij==jj)continue;
	    TLorentzVector* thatJet = (TLorentzVector*)fatjetP4->At(jj);
	    
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

    
    bool findAK8Match=false;
    //const float dRMax=0.4;
    //const float dRbMax=0.8;
    int matchedHAK8JetIndex[2]={-1,-1};
		      int AK8nJet=data.GetInt("AK8PuppinJet");
    for(int ij=0; ij<AK8nJet; ij++)
      {
	TLorentzVector* thisJet = (TLorentzVector*)AK8PuppijetP4->At(ij);

	for(int jj=0; jj<AK8nJet; jj++)
	  {

	    if(ij==jj)continue;
	    TLorentzVector* thatJet = (TLorentzVector*)AK8PuppijetP4->At(jj);
	    
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
		  matchedHAK8JetIndex[0]=ij;
		  matchedHAK8JetIndex[1]=jj;
		}
		else
		  {
		    matchedHAK8JetIndex[0]=jj;
		    matchedHAK8JetIndex[1]=ij;
		  }
		findAK8Match=true;
		break;
	      }

	    if(findAK8Match)break;

	  }	

	if(findAK8Match)break;

      }

    if(!findAK8Match)continue;
   
    
    
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
     	TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
	recoH_l4[i]= (*thisJet);
    	if(thisJet->Pt()<200)continue;
	if(fabs(thisJet->Eta())>2.4)continue;
	nGoodJets++;
      }

    if(nGoodJets<2)continue;
    nPass[5]++;
    
    if(debug)
      {
	recoH_l4[0].Print();
	recoH_l4[1].Print();
      }

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

	long eventId=data.GetLong64("eventId");
	
	if(AK8nJet<2)continue;
    // now plot mass
    for(int i=0; i<2;i++)
      {
		
		//cout<< eventId<<endl;
	int jet=matchedHJetIndex[i];
	int AK8jet=matchedHAK8JetIndex[i];
	TLorentzVector* thisJet = (TLorentzVector*)puppijetP4->At(jet);
	float thea_corr = getPUPPIweight(thisJet->Pt(),thisJet->Eta());
	float thea_mass = fatjetSDmass[jet]*thea_corr;
	
	//if(jet>AK8nJet-1)break;
	if(thisJet->Pt()>99998)break;
	//cout<< eventId<<endl;
	
	TLorentzVector* thisAK8Jet = (TLorentzVector*)AK8PuppijetP4->At(AK8jet);
	thea_corr = getPUPPIweight(thisAK8Jet->Pt(),thisAK8Jet->Eta());

	if(debug)
	  cout << thisJet->Pt() << "\t" << thisJet->Eta() << "\t" << thea_corr << endl;
	
	h_SD[i]->Fill(fatjetSDmass[jet]);
	h_SDCorr[i]->Fill(fatjetSDmassL2L3Corr[jet]);
	h_SDCorrThea[i]->Fill(thea_mass);
	h_PR[i]->Fill(fatjetPRmass[jet]);
	h_PRCorr[i]->Fill(fatjetPRmassL2L3Corr[jet]);
	
	h_AK8SD[i]->Fill(AK8PuppijetSDmass[AK8jet]);
	h_AK8SDCorrThea[i]->Fill(AK8PuppijetSDmass[AK8jet]*thea_corr);

	h_SD[2]->Fill(fatjetSDmass[jet]);
	h_SDCorr[2]->Fill(fatjetSDmassL2L3Corr[jet]);
	h_SDCorrThea[2]->Fill(thea_mass);
	h_PR[2]->Fill(fatjetPRmass[jet]);
	h_PRCorr[2]->Fill(fatjetPRmassL2L3Corr[jet]);
	
	h_AK8SD[2]->Fill(AK8PuppijetSDmass[AK8jet]);
	h_AK8SDCorrThea[2]->Fill(AK8PuppijetSDmass[AK8jet]*thea_corr);


	h_diff_SD[i]->Fill((fatjetSDmass[jet]-125)/125);
	h_diff_SDCorr[i]->Fill((fatjetSDmassL2L3Corr[jet]-125)/125);
	h_diff_SDCorrThea[i]->Fill((thea_mass-125)/125);
	h_diff_PR[i]->Fill((fatjetPRmass[jet]-125)/125);
	h_diff_PRCorr[i]->Fill((fatjetPRmassL2L3Corr[jet]-125)/125);
	
	h_diff_AK8SD[i]->Fill((AK8PuppijetSDmass[AK8jet]-125)/125);
	h_diff_AK8SDCorrThea[i]->Fill((AK8PuppijetSDmass[AK8jet]*thea_corr-125)/125);

	h_diff_SD[2]->Fill((fatjetSDmass[jet]-125)/125);
	h_diff_SDCorr[2]->Fill((fatjetSDmassL2L3Corr[jet]-125)/125);
	h_diff_SDCorrThea[2]->Fill((thea_mass-125)/125);
	h_diff_PR[2]->Fill((fatjetPRmass[jet]-125)/125);
	h_diff_PRCorr[2]->Fill((fatjetPRmassL2L3Corr[jet]-125)/125);
	
	h_diff_AK8SD[2]->Fill((AK8PuppijetSDmass[AK8jet]-125)/125);
	h_diff_AK8SDCorrThea[2]->Fill((AK8PuppijetSDmass[AK8jet]*thea_corr-125)/125);

      }
    

  } // end of loop over entries

  std::cout << "nTotal    = " << nTotal << std::endl;
  for(int i=0;i<20;i++)
    if(nPass[i]>0)
      std::cout << "nPass[" << i << "]= " << nPass[i] << std::endl;

  TFile* outFile = new TFile(outputFile.Data(),"recreate");

  for(int i=0; i<nHistos; i++)
    {
      h_diff_SD[i]->Write();
      h_diff_SDCorr[i]->Write();
      h_diff_SDCorrThea[i]->Write();
      h_diff_PR[i]->Write();
      h_diff_PRCorr[i]->Write();

      h_SD[i]->Write();
      h_SDCorr[i]->Write();
      h_SDCorrThea[i]->Write();
      h_PR[i]->Write();
      h_PRCorr[i]->Write();
	
	h_AK8SD[i]->Write();
	h_AK8SDCorrThea[i]->Write();
	
	h_diff_AK8SD[i]->Write();
	h_diff_AK8SDCorrThea[i]->Write();
    }

  outFile->Close();


}

void xAna_hh_massResolution(){
	
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
