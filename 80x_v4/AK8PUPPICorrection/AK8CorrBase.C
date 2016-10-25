
void AK8CorrBase(int wMs,int wM, string st,string st2,string option=""){	

	//1=signal ,0=QCD ,2=data
	int nameRoot=1;
	if((st2.find("QCD")!= std::string::npos)||
	(st2.find("bGen")!= std::string::npos)||
	(st2.find("bEnriched")!= std::string::npos))nameRoot=0;
	if(st2.find("data")!= std::string::npos)nameRoot=2;
	cout<<"nameRoot = "<<nameRoot<<endl;
	
	//option-----------------------------------------------------------
	
	int JESOption=0;
	if(option.find("JESUp")!= std::string::npos)JESOption=1;
	if(option.find("JESDown")!= std::string::npos)JESOption=2;
	cout<<"JESOption = "<<JESOption<<endl;
	
	
	TFile *f;
	TTree *tree;
	int nPass[20]={0};
	int total=0;
	
	TH1D* th1;
	th1=new TH1D("mass","mass",150,0,150);
	
	TH1D* th3;
	th3=new TH1D("mass","mass",1100,200,2400);
	
	double ptBins[14]={200,300,400,500,600,700,800,900,1000,1250,1500,1750,2000,2500};
	
	TH1D* th2[6][14];
	
	for(int i=0;i<14;i++){
		th2[0][i]=(TH1D*)th1->Clone(Form("genBarelMass%.0f",ptBins[i]));
		th2[1][i]=(TH1D*)th1->Clone(Form("genEndcapMass%.0f",ptBins[i]));
		th2[2][i]=(TH1D*)th1->Clone(Form("recoBarelMass%.0f",ptBins[i]));
		th2[3][i]=(TH1D*)th1->Clone(Form("recoEndcapMass%.0f",ptBins[i]));
		th2[4][i]=(TH1D*)th3->Clone(Form("ptBarel%.0f",ptBins[i]));
		th2[5][i]=(TH1D*)th3->Clone(Form("ptEndcap%.0f",ptBins[i]));
	}

	for(int i=0;i<14;i++){
		th2[0][i]->Sumw2();
		th2[1][i]->Sumw2();
		th2[2][i]->Sumw2();
		th2[3][i]->Sumw2();
		th2[4][i]->Sumw2();
		th2[5][i]->Sumw2();
	}
	
	//---------------------------------
	
	for (int w=wMs;w<wM;w++){
		if(w%20==0)cout<<w<<endl;
		
		if (nameRoot!=1)f = TFile::Open(Form("%s%d.root",st.data(),w));
		else f = TFile::Open(st.data());
		if (!f || !f->IsOpen())continue;
		
		TDirectory * dir;
		if (nameRoot!=1)dir = (TDirectory*)f->Get(Form("%s%d.root:/tree",st.data(),w));
		else dir = (TDirectory*)f->Get(Form("%s:/tree",st.data()));
		
		dir->GetObject("treeMaker",tree);
		TreeReader data(tree);
		total+=data.GetEntriesFast();
		for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
			if(jEntry%2)continue;
			data.GetEntry(jEntry);
			
			
			
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


			
			
			TClonesArray* AK8PuppijetP4 = (TClonesArray*) data.GetPtrTObject("AK8PuppijetP4");

			// check matching first    
			bool findAK8Match=false;
			const float dRMax=0.4;
			const float dRbMax=0.8;
			int matchedHAK8JetIndex[2]={-1,-1};
		      int AK8nJet=data.GetInt("AK8PuppinJet");
			if(AK8nJet<2)continue;
			bool matchb=1;
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
			
			
				
			
			
			Float_t*  AK8PuppijetGenSDmass = data.GetPtrFloat("AK8PuppijetGenSDmass");
			Float_t*  AK8PuppijetSDmass = data.GetPtrFloat("AK8PuppijetSDmass");
			int* AK8PuppinSubSDJet=data.GetPtrInt("AK8PuppinSubSDJet");
			
			 TLorentzVector *leadjet = (TLorentzVector*)AK8PuppijetP4->At(matchedHAK8JetIndex[0]);
			TLorentzVector *subljet = (TLorentzVector*)AK8PuppijetP4->At(matchedHAK8JetIndex[1]);
			
			 int n_neutrinos_lead=0;
    TLorentzVector neutrinos_lead_p4;

    int n_neutrinos_subl=0;
    TLorentzVector neutrinos_subl_p4;    

    for(int ig=0; ig < nGenPar; ig++){
      
      int pid =abs(genParId[ig]);
      int mom_pid = abs(genMomParId[ig]);
      if(genParSt[ig]!=1)continue;
      
      if(pid!=12 && pid!=14 && pid!=16)continue;
    
      // look for daughters of charm and b hadrons
      // bhadrons can decay to l nu + charm hadrons 
      // charm hadrons could further decay to l nu + other hadrons
      
      if(mom_pid < 400 || mom_pid > 600)continue;
      
      TLorentzVector* thisGen = (TLorentzVector*)genParP4->At(ig);
      
      if(thisGen->DeltaR(*leadjet)<0.8) { 
	n_neutrinos_lead++;
	neutrinos_lead_p4 += *thisGen; 
      }
      else if(thisGen->DeltaR(*subljet)<0.8) { 
	n_neutrinos_subl++;
	neutrinos_subl_p4 += *thisGen;
      }
      else continue;
      
    }
    
    if(n_neutrinos_lead > 0||  n_neutrinos_subl > 0)continue;
			
			for(int i=0; i<2;i++){
		
				
				int AK8jet=matchedHAK8JetIndex[i];
				
				TLorentzVector* thisAK8Jet = (TLorentzVector*)AK8PuppijetP4->At(AK8jet);
				
				if(thisAK8Jet->Pt()<200)continue;
				if(fabs(thisAK8Jet->Eta())>2.4)continue;
				if(AK8PuppinSubSDJet[AK8jet]!=2)continue;
				
				th1->Fill(AK8PuppijetGenSDmass[AK8jet]);
				
				for(int j=0;j<13;j++){
					if(thisAK8Jet->Pt()>ptBins[j] && thisAK8Jet->Pt()<ptBins[j+1]){
						
						if(fabs(thisAK8Jet->Eta())<1.3){
							th2[0][j]->Fill(AK8PuppijetGenSDmass[AK8jet]);
							th2[2][j]->Fill(AK8PuppijetSDmass[AK8jet]);
							th2[4][j]->Fill(thisAK8Jet->Pt());
						}
						else{
							th2[1][j]->Fill(AK8PuppijetGenSDmass[AK8jet]);
							th2[3][j]->Fill(AK8PuppijetSDmass[AK8jet]);
							th2[5][j]->Fill(thisAK8Jet->Pt());
						}
							
					}
				}
				
				if(thisAK8Jet->Pt()>ptBins[13]){
					if(fabs(thisAK8Jet->Eta())<1.3){
							th2[0][13]->Fill(AK8PuppijetGenSDmass[AK8jet]);
							th2[2][13]->Fill(AK8PuppijetSDmass[AK8jet]);
							th2[4][13]->Fill(thisAK8Jet->Pt());
						}
						else{
							th2[1][13]->Fill(AK8PuppijetGenSDmass[AK8jet]);
							th2[3][13]->Fill(AK8PuppijetSDmass[AK8jet]);
							th2[5][13]->Fill(thisAK8Jet->Pt());
						}
				}
				
			}
			
		}
	}	
	cout<<"entries="<<total<<endl;	
	
	
	TFile* outFile ;
	outFile= new TFile(Form("corr2/%s.root",st2.data()),"recreate");
	th1->Write();
	for(int i=0;i<14;i++){
		th2[0][i]->Write();
		th2[1][i]->Write();
		th2[2][i]->Write();
		th2[3][i]->Write();
		th2[4][i]->Write();
		th2[5][i]->Write();
	}
	outFile->Close();
}
