
void AK8CorrBase(int wMs,int wM, string st,string st2,string option=""){	

	//1=signal ,0=QCD ,2=data
	int nameRoot=1;
	if((st2.find("QCD")!= std::string::npos)||
	(st2.find("bGen")!= std::string::npos)||
	(st2.find("bEnriched")!= std::string::npos))nameRoot=0;
	if(st2.find("data")!= std::string::npos)nameRoot=2;
	
	bool fixGen=0;
	if(st2.find("B1000")!= std::string::npos)fixGen=1;
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
	th3=new TH1D("mass","mass",1500,200,3200);
	
	//double ptBins[14]={200,300,400,500,600,700,800,900,1000,1250,1500,1750,2000,2500};
	
	TH1D* th2[4];
	
	th2[0]=(TH1D*)th1->Clone("recoBarrelMass");
	th2[1]=(TH1D*)th1->Clone("recoEndcapMass");
	th2[2]=(TH1D*)th3->Clone("ptBarrel");
	th2[3]=(TH1D*)th3->Clone("ptEndcap");
	
	for(int i=0;i<4;i++){
		//th2[0][i]->Sumw2();
		//th2[1][i]->Sumw2();
		th2[i]->Sumw2();
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


			
			
			TClonesArray* FATjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");

			// check matching first    
			bool findMatch=false;
			const float dRMax=0.4;
			const float dRbMax=0.8;
			int matchedHAK8JetIndex[2]={-1,-1};
		      int FATnJet=data.GetInt("FATnJet");
			if(FATnJet<2)continue;
			bool matchb=1;
			for(int ij=0; ij<FATnJet; ij++)
			{
				TLorentzVector* thisJet = (TLorentzVector*)FATjetP4->At(ij);

				for(int jj=0; jj<FATnJet; jj++)
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
					
					if(ij<jj){
					matchedHAK8JetIndex[0]=ij;
					matchedHAK8JetIndex[1]=jj;
					}
					else
					{
					matchedHAK8JetIndex[0]=jj;
					matchedHAK8JetIndex[1]=ij;
					}
					findMatch=true;
					break;
					}

				if(findMatch)break;
	
				}	

				if(findMatch)break;

			}
	
			if(!findMatch)continue;
			
			Float_t*  AK8PuppijetSDmass = data.GetPtrFloat("AK8PuppijetSDmass");
			int* AK8PuppinSubSDJet=data.GetPtrInt("AK8PuppinSubSDJet");
			
			TLorentzVector* thisJet,* thatJet;
			thisJet= (TLorentzVector*)FATjetP4->At(0);
			thatJet= (TLorentzVector*)FATjetP4->At(1);
			if(thisJet->Pt()<200||thatJet->Pt()<200)continue;
			if(fabs(thisJet->Eta())>2.4||fabs(thatJet->Eta())>2.4)continue;
			
			double SDMass[2];
			
			Int_t AK8PuppinJet        = data.GetInt("AK8PuppinJet");
			bool matchThis=0,matchThat=0;
			TClonesArray* AK8PuppijetP4 = (TClonesArray*) data.GetPtrTObject("AK8PuppijetP4");
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
			
			if(fabs(thisJet->Eta())<1.3){
				th2[0]->Fill(SDMass[0]);
				th2[2]->Fill(thisJet->Pt());
			}
			else {
				th2[1]->Fill(SDMass[0]);
				th2[3]->Fill(thisJet->Pt());
			}
			
			if(fabs(thatJet->Eta())<1.3){
				th2[0]->Fill(SDMass[1]);
				th2[2]->Fill(thatJet->Pt());
			}
			else {
				th2[1]->Fill(SDMass[1]);
				th2[3]->Fill(thatJet->Pt());
			}
			

		}
	}	
	cout<<"entries="<<total<<endl;	
	TFile* outFile ;
	outFile= new TFile(Form("corr/%s.root",st2.data()),"recreate");
	th1->Write();
	for(int i=0;i<4;i++){
		th2[i]->Write();
	}
	outFile->Close();
}
