void HHbbbbMassAnalyzerBase(int wMs,int wM, string st,string st2,int option=0){	
	//0=signal ,1=QCD ,2=data-----------------------------------------------------------
	int nameRoot=1;
	if(st2.find("QCD")!= std::string::npos)nameRoot=0;
	if(st2.find("data")!= std::string::npos)nameRoot=2;
	
	bool isWW=0;
	if(st.find("WW")!= std::string::npos)isWW=1;
	
	
	//tuple tree and cutflow variables------------------------------------------------------------------------------------
	TFile *f2;
	if(isWW) f2=TFile::Open("TProfileWW.root");
	else f2=TFile::Open("TProfile.root");
	
	TProfile* tp1,*tp2;
	if (option==1){
		tp1=new TProfile("JECPR","JECPR",10,200,2200,0,0.2);
		tp2=new TProfile("JECSDPuppi","JECSDPuppi",10,200,2200,0,0.2);
	}
	else  {
		tp1= (TProfile *)f2->FindObjectAny("JECPR");
		tp2= (TProfile *)f2->FindObjectAny("JECSDPuppi");
	}
	TFile *f;
	TTree *tree;
	int nPass[20]={0},total=0;
	double nPassB[6]={0};
	
	TH1D* th5[20];
	th5[0]=new TH1D("FATjetPRmass","FATjetPRmass",90,0,180);
	th5[1]=new TH1D("FATjetPRmassL2L3Corr","FATjetPRmassL2L3Corr",90,0,180);
	th5[2]=new TH1D("FATjetPuppiSDmass","FATjetPuppiSDmass",90,0,180);
	th5[3]=new TH1D("FATjetPuppiSDmassL2L3Corr","FATjetPuppiSDmassL2L3Corr",90,0,180);
	
	th5[4]=new TH1D("FATjetPRmassS","FATjetPRmass",90,0,180);
	th5[5]=new TH1D("FATjetPRmassL2L3CorrS","FATjetPRmassL2L3Corr",90,0,180);
	th5[6]=new TH1D("FATjetPuppiSDmassS","FATjetPuppiSDmass",90,0,180);
	th5[7]=new TH1D("FATjetPuppiSDmassL2L3CorrS","FATjetPuppiSDmassL2L3Corr",90,0,180);
	
	th5[8]=new TH1D("PV","PV",30,0,60);
	th5[9]=new TH1D("PVFATjetPRmass","PV",30,0,60);
	th5[10]=new TH1D("PVFATjetPRmassL2L3Corr","PV",30,0,60);
	th5[11]=new TH1D("PVFATjetPuppiSDmass","PV",30,0,60);
	th5[12]=new TH1D("PVFATjetPuppiSDmassL2L3Corr","PV",30,0,60);
	th5[13]=new TH1D("PVFATjetPRmassTau21","PV",30,0,60);
	th5[14]=new TH1D("PVFATjetPuppiSDmassPuppiTau21","PV",30,0,60);
	th5[15]=new TH1D("PVFATjetPRmassL2L3CorrTau21","PV",30,0,60);
	th5[16]=new TH1D("PVFATjetPuppiSDmassL2L3CorrPuppiTau21","PV",30,0,60);
	
	th5[17]=new TH1D("tau21","tau21",20,0,1);
	th5[18]=new TH1D("tau21Puppi","tau21Puppi",20,0,1);
	
	//NCUtuple loop----------------------------------------------------------------------------------------
	for (int w=wMs;w<wM;w++){
		if(w%20==0)cout<<w<<endl;
		//Get ntuple----------------------------------------------------------------------------------------
		cout<<st.data()<<endl;
		f = TFile::Open(st.data());
		if (!f || !f->IsOpen())continue;
		TDirectory * dir;
		dir = (TDirectory*)f->Get(Form("%s:/tree",st.data()));
		dir->GetObject("treeMaker",tree);
		TreeReader data(tree);
		total+=data.GetEntriesFast();
		
		for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){//event loop----------------------------------------------------------------------------------------
			data.GetEntry(jEntry);
			
			Int_t nGenPar        = data.GetInt("nGenPar");
			Int_t* genParId      = data.GetPtrInt("genParId");
			Int_t* genParSt      = data.GetPtrInt("genParSt");
			Int_t* genMomParId   = data.GetPtrInt("genMomParId");

			int genHIndex[2]={-1,-1};
			const int bosonID= isWW?  24:25;
			
				for(int ig=0; ig < nGenPar; ig++){

      if(abs(genParId[ig])!=bosonID)continue;

      

      if(genHIndex[0]<0)
	genHIndex[0]=ig;

      else if(genHIndex[1]<0)
	genHIndex[1]=ig;

    }    
    

    if(genHIndex[0]<0 || genHIndex[1]<0)continue;
    

    if(genHIndex[0]==genHIndex[1])continue;
	
	 bool hasLepton=false;
    Int_t* genDa1      = data.GetPtrInt("genDa1");
    Int_t* genDa2      = data.GetPtrInt("genDa2");

    for(int ig=0; ig < nGenPar; ig++){

      if(abs(genParId[ig])!=bosonID)continue;
      int da1=genDa1[ig];
      int da2=genDa2[ig];

      if(da1<0 || da2<0)continue;
      int da1pdg = genParId[da1];
      int da2pdg = genParId[da2];

      if(abs(da1pdg)>10 && abs(da1pdg)<17)
       	hasLepton=true;
      if(abs(da2pdg)>10 && abs(da2pdg)<17)
       	hasLepton=true;

      if(hasLepton)break;

    }

    if(hasLepton)continue;
	
	
	nPass[0]++;
	
	TLorentzVector genH_l4[2];
    TClonesArray* genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");

    
    for(int ih=0; ih<2; ih++)
      genH_l4[ih] = *((TLorentzVector*)genParP4->At(genHIndex[ih]));
	
			
			Int_t nVtx        = data.GetInt("nVtx");
			//0. has a good vertex
			if(nVtx<1)continue;nPass[0]++;
			//1.trigger
			std::string* trigName = data.GetPtrString("hlt_trigName");
		 	vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
			bool passTrigger=false;
			for(int it=0; it< data.GetPtrStringSize(); it++){
				std::string thisTrig= trigName[it];
				bool results = trigResult[it];
				//if(trigResult[it]==1)cout<<it<<","<<thisTrig<<","<<trigResult[it]<<endl;
				if( ((thisTrig.find("HLT_PFHT800")!= std::string::npos||
						thisTrig.find("HLT_PFHT650")!= std::string::npos||
						thisTrig.find("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v")!= std::string::npos||
						thisTrig.find("HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v")!= std::string::npos||
						thisTrig.find("HLT_AK8PFJet360_TrimMass30_v")!= std::string::npos||
						thisTrig.find("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v")!= std::string::npos
						) && results==1)){
					passTrigger=true;
					break;
				}
			}
			//if(!passTrigger)continue;
			nPass[1]++;
		
			int nFATJet         = data.GetInt("FATnJet");
			TClonesArray* fatjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
			vector<bool>    &FATjetPassIDTight = *((vector<bool>*) data.GetPtr("FATjetPassIDTight"));
			//2.nJets
			if(nFATJet<2)continue;nPass[2]++;
			TLorentzVector* thisJet ,* thatJet;
			thisJet= (TLorentzVector*)fatjetP4->At(0);
			thatJet = (TLorentzVector*)fatjetP4->At(1);
			//3. Pt 
			if(thisJet->Pt()<300||thatJet->Pt()<300)continue;
			nPass[3]++;
			//4tightId-----------------------------------------
			if(FATjetPassIDTight[0]==0||FATjetPassIDTight[1]==0)continue;
			Float_t*  FATjetCEmEF = data.GetPtrFloat("FATjetCEmEF");
			Float_t*  FATjetMuEF = data.GetPtrFloat("FATjetMuEF");
			if(FATjetMuEF[0]>0.8||FATjetMuEF[1]>0.8)continue;
			if(FATjetCEmEF[0]>0.9||FATjetCEmEF[1]>0.9)continue;
			nPass[4]++;
			//5. Eta-----------------------------------------
			if(fabs(thisJet->Eta())>2.4||fabs(thatJet->Eta())>2.4)continue;
			nPass[5]++;
			//6. DEta-----------------------------------------
			float dEta = fabs(thisJet->Eta()-thatJet->Eta());
			//if(dEta>1.3)continue;
			nPass[6]++;
			//7. Mjj-----------------------------------------
			float mjj = (*thisJet+*thatJet).M();
			if(mjj<1000)continue;
			nPass[7]++;
			//8. fatjetPRmassL2L3Corr-----------------------------------------
			Float_t*  fatjetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
			//if(fatjetPRmassL2L3Corr[0]<massDown||fatjetPRmassL2L3Corr[0]>massUp)continue;
			//if(fatjetPRmassL2L3Corr[1]<massDown||fatjetPRmassL2L3Corr[1]>massUp)continue;
			nPass[8]++;
			//9.-----------------------------------------
			
			
			vector<float>   *subjetSDCSV =  data.GetPtrVectorFloat("FATsubjetSDCSV");
			if(subjetSDCSV[0][0]<0.46 && isWW==0)continue;
			if(subjetSDCSV[0][1]<0.46 && isWW==0)continue;
			if(subjetSDCSV[1][0]<0.46 && isWW==0)continue;
			if(subjetSDCSV[1][1]<0.46 && isWW==0)continue;
			nPass[10]++;
			
			
			Float_t*  FATjetPRmass = data.GetPtrFloat("FATjetPRmass");
			Float_t*  FATjetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
			Float_t*  FATjetPuppiSDmass = data.GetPtrFloat("FATjetPuppiSDmass");
			Float_t*  FATjetPuppiSDmassL2L3Corr = data.GetPtrFloat("FATjetPuppiSDmassL2L3Corr");
			
			th5[8]->Fill(nVtx);
			tp1->Fill(thisJet->Pt(),((FATjetPRmassL2L3Corr[0]-FATjetPRmass[0])/FATjetPRmass[0]));
			tp2->Fill(thisJet->Pt(),((FATjetPuppiSDmassL2L3Corr[0]-FATjetPuppiSDmass[0])/FATjetPuppiSDmass[0]));
			
			Float_t*  fatjetTau1 = data.GetPtrFloat("FATjetTau1");
			Float_t*  fatjetTau2 = data.GetPtrFloat("FATjetTau2");
			
			Float_t*  FATjetPuppiTau1 = data.GetPtrFloat("FATjetPuppiTau1");
			Float_t*  FATjetPuppiTau2 = data.GetPtrFloat("FATjetPuppiTau2");
			
			double tau21_1=(fatjetTau2[0]/fatjetTau1[0]),tau21_2=(fatjetTau2[1]/fatjetTau1[1]);
			
			
			double massUp=135,massDown=105;
			if(isWW)massUp=105;
			if(isWW)massDown=65;
			
			if(massDown<FATjetPRmass[0] && FATjetPRmass[0]<massUp)th5[9]->Fill(nVtx);
			if(massDown<FATjetPRmassL2L3Corr[0] && FATjetPRmassL2L3Corr[0]<massUp)th5[10]->Fill(nVtx);
			if(massDown<FATjetPuppiSDmass[0] && FATjetPuppiSDmass[0]<massUp)th5[11]->Fill(nVtx);
			if(massDown<FATjetPuppiSDmassL2L3Corr[0] && FATjetPuppiSDmassL2L3Corr[0]<massUp)th5[12]->Fill(nVtx);
			
		
			double puppitau21_1=(FATjetPuppiTau2[0]/FATjetPuppiTau1[0]);
			
			th5[17]->Fill(tau21_1);
			th5[18]->Fill(puppitau21_1);
			if( puppitau21_1<0.4&&massDown<FATjetPuppiSDmass[0] && FATjetPuppiSDmass[0]<massUp)th5[14]->Fill(nVtx);
			if( puppitau21_1<0.4&&massDown<FATjetPuppiSDmassL2L3Corr[0] && FATjetPuppiSDmassL2L3Corr[0]<massUp)th5[16]->Fill(nVtx);
			
			if(tau21_1<0.4){
				if(massDown<FATjetPRmass[0] && FATjetPRmass[0]<massUp)th5[13]->Fill(nVtx);
				if(massDown<FATjetPRmassL2L3Corr[0] && FATjetPRmassL2L3Corr[0]<massUp)th5[15]->Fill(nVtx);
				
			}
			
			
			
			if(dEta>1.3)continue;
			
			bool thisJetMatchZero=thisJet->DeltaR(genH_l4[0])<0.4 ?1:0;
			
			if((thisJet->DeltaR(genH_l4[1])<0.4 )||(thisJet->DeltaR(genH_l4[0])<0.4 )){
				th5[0]->Fill(FATjetPRmass[0]);
				th5[1]->Fill(FATjetPRmassL2L3Corr[0]);
				th5[2]->Fill(FATjetPuppiSDmass[0]);
				th5[3]->Fill(FATjetPuppiSDmassL2L3Corr[0]);
				
			}
			
			
			if((thatJet->DeltaR(genH_l4[1])<0.4 )||(thatJet->DeltaR(genH_l4[0])<0.4 && !thisJetMatchZero )){
			th5[4]->Fill(FATjetPRmass[1]);
			th5[5]->Fill(FATjetPRmassL2L3Corr[1]);
			th5[6]->Fill(FATjetPuppiSDmass[1]);
			th5[7]->Fill(FATjetPuppiSDmassL2L3Corr[1]);
			
			}
		}//end event loop----------------------------------------------------------------------------------------
	}	//end ntuple loop----------------------------------------------------------------------------------------
	
	cout<<"total="<<total<<endl;
	for(int i=0;i<16;i++)cout<<"nPass["<<i<<"]="<<nPass[i]<<endl;
	
	TFile* outFile = new TFile(Form("output/%s.root",st2.data()),"recreate");
	for(int i=0;i<19;i++)th5[i]->Write();
	outFile->Close();
	
	TFile* outFile2;
if(isWW)outFile2= new TFile("TProfileWW.root","recreate");
else 	outFile2= new TFile("TProfile.root","recreate");
	tp1->Write();
	tp2->Write();
	outFile2->Close();
	
	if(option==2){
		TFile* outFile3;
if(isWW)outFile3= new TFile("output/TProfileWW.root","recreate");
else 	outFile3= new TFile("output/TProfile.root","recreate");
		tp1->Write();
		tp2->Write();
		outFile3->Close();
	}
}