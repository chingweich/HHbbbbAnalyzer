
void HHbbbbBtagAdditionThinJet(int wMs,int wM, string st,string st2,double Xsec,int nameRoot=0){	
	TFile *f;
	TTree *tree;
	int nPass[20]={0};int total=0;
	int dataPassingcsc=0;
	
	TH1D* th1[4];
	
	th1[0]=new TH1D("addThinJet_3b","addThinJet_3b",10,-0.5,9.5);
	th1[1]=new TH1D("addThinJet_4b","addThinJet_4b",10,-0.5,9.5);
	
	th1[2]=new TH1D("addThinJet_3b_withoutDr","addThinJet_3b_withoutDr",10,-0.5,9.5);
	th1[3]=new TH1D("addThinJet_4b_withoutDr","addThinJet_4b_withoutDr",10,-0.5,9.5);
	
	std::vector<TString> eventlist;                                                                                                                                            
	if(nameRoot==2){
		std::vector<TString> eventlist;                                                                                                                                            
		std::ifstream fin("somfilename.txt");
		eventlist.clear();
		std::string line;
		if (fin.is_open())  {
			while (fin.good())  {
				getline(fin,line);
				eventlist.push_back(line);
		}
		fin.close();
		}
	}
		
	for (int w=wMs;w<wM;w++){
		
		if (nameRoot!=1)f = TFile::Open(Form("%s%d.root",st.data(),w));
		else f = TFile::Open(st.data());
		if (!f || !f->IsOpen())continue;
		
		TDirectory * dir;
		
		if (nameRoot!=1)dir = (TDirectory*)f->Get(Form("%s%d.root:/tree",st.data(),w));
		else dir = (TDirectory*)f->Get(Form("%s:/tree",st.data()));
		
		if(w%20==0)cout<<w<<endl;
		
		dir->GetObject("treeMaker",tree);
		TreeReader data(tree);
		total+=data.GetEntriesFast();
		for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
			data.GetEntry(jEntry);
			
			
			//Int_t event        = data.GetInt("eventId");
			//0. has a good vertex
			Int_t nVtx        = data.GetInt("nVtx");
			if(nVtx<1)continue;
			nPass[0]++;
	
			std::string* trigName = data.GetPtrString("hlt_trigName");
		 	vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
		 	const Int_t nsize = data.GetPtrStringSize();

			bool passTrigger=false;
			for(int it=0; it< nsize; it++){
				std::string thisTrig= trigName[it];
				bool results = trigResult[it];
				//cout<<it<<"="<<thisTrig<<endl;
				// std::cout << thisTrig << " : " << results << std::endl;
				if( ((thisTrig.find("HLT_PFHT800")!= std::string::npos||
						thisTrig.find("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v")!= std::string::npos||
						thisTrig.find("HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v")!= std::string::npos||
						thisTrig.find("HLT_AK8PFJet360_TrimMass30_v")!= std::string::npos||
						thisTrig.find("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v")!= std::string::npos
						) && results==1)){
					passTrigger=true;
					break;
				}
			}
			if(!passTrigger && (nameRoot!=0))continue;
			nPass[1]++;
		
			int nFATJet         = data.GetInt("FATnJet");
			const int nJets=nFATJet;
			TClonesArray* fatjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
			Float_t*  fatjetTau1 = data.GetPtrFloat("FATjetTau1");
			Float_t*  fatjetTau2 = data.GetPtrFloat("FATjetTau2");
			//Float_t*  fatjetTau4 = data.GetPtrFloat("FATjetTau4");
			Float_t*  fatjetCISVV2 = data.GetPtrFloat("FATjetCISVV2");
			Float_t*  fatjetPRmass = data.GetPtrFloat("FATjetPRmass");
			Float_t*  fatjetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
			Float_t*  fatjetSDmass = data.GetPtrFloat("FATjetSDmass");
			Int_t*   nSubSoftDropJet = data.GetPtrInt("FATnSubSDJet");
			vector<float>   *subjetSDCSV =  data.GetPtrVectorFloat("FATsubjetSDCSV");
			//  vector<float>   *subjetSDCSV =  data.GetPtrVectorFloat("FATsubjetSDCSV", nFATJet);
			vector<float>   *subjetSDPx  =  data.GetPtrVectorFloat("FATsubjetSDPx", nFATJet);
			vector<float>   *subjetSDPy  =  data.GetPtrVectorFloat("FATsubjetSDPy", nFATJet);
			vector<float>   *subjetSDPz  =  data.GetPtrVectorFloat("FATsubjetSDPz", nFATJet);
			vector<float>   *subjetSDE   =  data.GetPtrVectorFloat("FATsubjetSDE", nFATJet);
			vector<Int_t>   *FATsubjetSDHadronFlavor =  data.GetPtrVectorInt("FATsubjetSDHadronFlavor");
			//vector<bool>    &passFatJetLooseID = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
			vector<bool>    &FATjetPassIDTight = *((vector<bool>*) data.GetPtr("FATjetPassIDTight"));
			if(nJets<2)continue;
			//-----------------------------------------
			TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(0);
			TLorentzVector* thatJet = (TLorentzVector*)fatjetP4->At(1);
			//2. Pt and tightId-----------------------------------------
			if(thisJet->Pt()<200)continue;
			if(thatJet->Pt()<200)continue;
			if(FATjetPassIDTight[0]==0)continue;
			if(FATjetPassIDTight[1]==0)continue;
			nPass[2]++;
			//3. Eta-----------------------------------------
			if(fabs(thisJet->Eta())>2.4)continue;
			if(fabs(thatJet->Eta())>2.4)continue;
			nPass[3]++;
			//4. DEta-----------------------------------------
			float dEta = fabs(thisJet->Eta()-thatJet->Eta());
			if(dEta>1.3)continue;
			nPass[4]++;
			//5. Mjj-----------------------------------------
			float mjj = (*thisJet+*thatJet).M();
			if(mjj<1000)continue;
			nPass[5]++;
			//6. fatjetPRmassL2L3Corr-----------------------------------------
			if(fatjetPRmassL2L3Corr[0]<105||fatjetPRmassL2L3Corr[0]>135)continue;
			if(fatjetPRmassL2L3Corr[1]<105||fatjetPRmassL2L3Corr[1]>135)continue;
			nPass[6]++;
			//7.-----------------------------------------
			double tau21_1=(fatjetTau2[0]/fatjetTau1[0]),
		           tau21_2=(fatjetTau2[1]/fatjetTau1[1]);
			if(tau21_1>0.75||tau21_2>0.75)continue;
			if(tau21_1>0.6 &&tau21_2>0.6) continue;
			nPass[7]++;
			
			bool isHPHP=0;
			if(tau21_1<0.6 && tau21_2<0.6 )isHPHP=1;
			
			if(nameRoot==2){
				//data:csc2015
				TString thisEvent;                                                   
				Int_t  runId = data.GetInt("runId");
				Int_t  lumiSection = data.GetInt("lumiSection");
				Int_t  eventId = data.GetInt("eventId");
				thisEvent.Form("%d:%d:%d",runId,lumiSection,eventId);                                                                                                                 
				if ( std::find(eventlist.begin(), eventlist.end(), thisEvent) != eventlist.end()){
					std::cout<<" match found to skip event "<<std::endl;
					continue;
				}                                                                              
				dataPassingcsc++;
			}
			
			
			//8.btag
			int nbtag=0;
			if(subjetSDCSV[0][0]>0.605)nbtag++;
			if(subjetSDCSV[0][1]>0.605)nbtag++;
			if(subjetSDCSV[1][0]>0.605)nbtag++;
			if(subjetSDCSV[1][1]>0.605)nbtag++;
		
			int THINnJet         = data.GetInt("THINnJet");
			TClonesArray* THINgenjetP4 = (TClonesArray*) data.GetPtrTObject("THINgenjetP4");
			vector<bool>    &THINjetPassIDTight = *((vector<bool>*) data.GetPtr("THINjetPassIDTight"));
			int AddTHINnJet=0;
			int AddTHINnJet2=0;
			for(int i=0;i<THINnJet;i++){
				if(THINjetPassIDTight[i]==0)continue;
				TLorentzVector* thisThinJet = (TLorentzVector*)THINgenjetP4->At(i);
				AddTHINnJet2++;
				if(thisThinJet->DeltaR(*thisJet)<1)continue;
				if(thisThinJet->DeltaR(*thatJet)<1)continue;
				if(thisThinJet->Pt()<30)continue;
				if(fabs(thisThinJet->Eta())<3)continue;
				AddTHINnJet++;
				//cout<<w<<",etr="<<jEntry<<","<<i<<","<<"1="<<thisThinJet->DeltaR(*thisJet)<<",2="<<thisThinJet->DeltaR(*thatJet)<<endl;
			}
			
			if(nbtag==3)th1[0]->Fill(AddTHINnJet);
			if(nbtag==4)th1[1]->Fill(AddTHINnJet);
			
			if(nbtag==3)th1[2]->Fill(AddTHINnJet2);
			if(nbtag==4)th1[3]->Fill(AddTHINnJet2);
		}
	}	
	cout<<"entries="<<total<<endl;	
	if(nameRoot==2)cout<<"dataPassingcsc="<<dataPassingcsc<<endl;
	for(int i=0;i<14;i++)cout<<"nPass["<<i<<"]="<<nPass[i]<<endl;
	
	TFile* outFile = new TFile(Form("root_files_addThinJet/%s_3.root",st2.data()),"recreate");
	if(nameRoot==0){
		for(int i=0;i<4;i++){
			th1[i]->Sumw2();
			th1[i]->Scale(2245.87*Xsec/total);
		}
	}
	for (int i=0;i<4;i++)th1[i]->Write();
	outFile->Close();
	
}
