void HHbbbbBtagEffBaseO(int wMs,int wM, string st,string st2,double Xsec,int nameRoot=0){	
	TFile *f1;
	if(nameRoot==2)f1=TFile::Open("root_files_btaggedEff/data.root");
	else f1=TFile::Open(Form("root_files_btaggedEff/%s.root",st2.data()));
	TH1D* th4[5];
	
	
	th4[0]=new TH1D("weight_0b","weight_0b",40,0,2);
	th4[1]=new TH1D("weight_1b","weight_1b",40,0,2);
	th4[2]=new TH1D("weight_2b","weight_2b",40,0,2);
	th4[3]=new TH1D("weight_3b","weight_3b",40,0,2);
	th4[4]=new TH1D("weight_4b","weight_4b",40,0,2);
	
	
	BTagCalibration calib("CSVv2L", "CSVv2.csv");
	BTagCalibrationReader LF(&calib,               // calibration instance
                             BTagEntry::OP_LOOSE,  // operating point
                             "comb",               // measurement type
                             "central");           // systematics type
	
	BTagCalibrationReader HF(&calib,               // calibration instance
                             BTagEntry::OP_LOOSE,  // operating point
                             "mujets",               // measurement type
                             "central");           // systematics type
	
	
	TFile *f;
	TTree *tree;
	int nPass[20]={0};int total=0;
	int dataPassingcsc=0;
	double nPassB[5]={0};
	
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
			if(!passTrigger )continue;
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
			
			
			//8.btag

			TLorentzVector* thisSub1=new TLorentzVector(0,0,0,0);
			TLorentzVector* thisSub2=new TLorentzVector(0,0,0,0);
			TLorentzVector* thatSub1=new TLorentzVector(0,0,0,0);
			TLorentzVector* thatSub2=new TLorentzVector(0,0,0,0);
			thisSub1->SetPxPyPzE(subjetSDPx[0][0],subjetSDPy[0][0],subjetSDPz[0][0],subjetSDE[0][0]);
			thisSub2->SetPxPyPzE(subjetSDPx[0][1],subjetSDPy[0][1],subjetSDPz[0][1],subjetSDE[0][1]);
			thatSub1->SetPxPyPzE(subjetSDPx[1][0],subjetSDPy[1][0],subjetSDPz[1][0],subjetSDE[1][0]);
			thatSub2->SetPxPyPzE(subjetSDPx[1][1],subjetSDPy[1][1],subjetSDPz[1][1],subjetSDE[1][1]);
			
			
			double sf1=1,sf2=1,sf3=1,sf4=1;
			
			
			float MaxBJetPt = 670., MaxLJetPt = 1000.;
			
			double pt1=thisSub1->Pt(),pt2=thisSub2->Pt(),pt3=thatSub1->Pt(),pt4=thatSub2->Pt();
			if(FATsubjetSDHadronFlavor[0][0]!=0 && pt1>MaxBJetPt )pt1=MaxBJetPt-0.1;
			if(FATsubjetSDHadronFlavor[0][0]==0 && pt1>MaxLJetPt )pt1=MaxLJetPt-0.1;
			
			if(FATsubjetSDHadronFlavor[0][1]!=0 && pt2>MaxBJetPt )pt2=MaxBJetPt-0.1;
			if(FATsubjetSDHadronFlavor[0][1]==0 && pt2>MaxLJetPt )pt2=MaxLJetPt-0.1;
			
			if(FATsubjetSDHadronFlavor[1][0]!=0 && pt3>MaxBJetPt )pt3=MaxBJetPt-0.1;
			if(FATsubjetSDHadronFlavor[1][0]==0 && pt3>MaxLJetPt )pt3=MaxLJetPt-0.1;
			
			if(FATsubjetSDHadronFlavor[1][1]!=0 && pt4>MaxBJetPt )pt4=MaxBJetPt-0.1;
			if(FATsubjetSDHadronFlavor[1][1]==0 && pt4>MaxLJetPt )pt4=MaxLJetPt-0.1;
			
			//cout<<pt1<<","<<pt2<<","<<pt3<<","<<pt4<<endl;
			//cout<<thisSub1->Eta()<<","<<thisSub2->Eta()<<","<<thatSub1->Eta()<<","<<thatSub2->Eta()<<endl;
			
			if(FATsubjetSDHadronFlavor[0][0]==5)sf1=HF.eval(BTagEntry::FLAV_B,thisSub1->Eta(),pt1); 
			else if(FATsubjetSDHadronFlavor[0][0]==4)sf1=HF.eval(BTagEntry::FLAV_C,thisSub1->Eta(),pt1); 
			else sf1=LF.eval(BTagEntry::FLAV_UDSG,thisSub1->Eta(),pt1); 
			
			if(FATsubjetSDHadronFlavor[0][1]==5)sf2=HF.eval(BTagEntry::FLAV_B,thisSub2->Eta(),pt2); 
			else if(FATsubjetSDHadronFlavor[0][1]==4)sf2=HF.eval(BTagEntry::FLAV_C,thisSub2->Eta(),pt2); 
			else sf2=LF.eval(BTagEntry::FLAV_UDSG,thisSub2->Eta(),pt2); 
			
			if(FATsubjetSDHadronFlavor[1][0]==5)sf3=HF.eval(BTagEntry::FLAV_B,thatSub1->Eta(),pt3); 
			else if(FATsubjetSDHadronFlavor[1][0]==4)sf3=HF.eval(BTagEntry::FLAV_C,thatSub1->Eta(),pt3); 
			else sf3=LF.eval(BTagEntry::FLAV_UDSG,thatSub1->Eta(),pt3); 
			
			if(FATsubjetSDHadronFlavor[1][1]==5)sf4=HF.eval(BTagEntry::FLAV_B,thatSub2->Eta(),pt4); 
			else if(FATsubjetSDHadronFlavor[1][1]==4)sf4=HF.eval(BTagEntry::FLAV_C,thatSub2->Eta(),pt4); 
			else sf4=LF.eval(BTagEntry::FLAV_UDSG,thatSub2->Eta(),pt4); 
			
			int nbtag=0;
			if(subjetSDCSV[0][0]>0.605)nbtag++;
			if(subjetSDCSV[0][1]>0.605)nbtag++;
			if(subjetSDCSV[1][0]>0.605)nbtag++;
			if(subjetSDCSV[1][1]>0.605)nbtag++;
			
			if(nbtag==0){
				nPassB[0]+=(1-sf1)*(1-sf2)*(1-sf3)*(1-sf4);
				th4[0]->Fill((1-sf1)*(1-sf2)*(1-sf3)*(1-sf4));
			}
			if(nbtag==1){
				double temp=(sf1)*(1-sf2)*(1-sf3)*(1-sf4)
				+(1-sf1)*(sf2)*(1-sf3)*(1-sf4)
				+(1-sf1)*(1-sf2)*(sf3)*(1-sf4)
				+(1-sf1)*(1-sf2)*(1-sf3)*(sf4);
				nPassB[1]+=temp;
				th4[1]->Fill(temp);
			}
			if(nbtag==2){
				double temp=(sf1)*(sf2)*(1-sf3)*(1-sf4)
				+(sf1)*(1-sf2)*(sf3)*(1-sf4)
				+(sf1)*(1-sf2)*(1-sf3)*(sf4)
				+(1-sf1)*(sf2)*(sf3)*(1-sf4)
				+(1-sf1)*(sf2)*(1-sf3)*(sf4)
				+(1-sf1)*(1-sf2)*(sf3)*(sf4);
				nPassB[2]+=temp;
				th4[2]->Fill(temp);
			}
			if(nbtag==3){
				double temp=(sf1)*(sf2)*(sf3)*(1-sf4)
				+(sf1)*(sf2)*(1-sf3)*(sf4)
				+(sf1)*(1-sf2)*(sf3)*(sf4)
				+(1-sf1)*(sf2)*(sf3)*(sf4);
				nPassB[3]+=temp;
				th4[3]->Fill(temp);
			}
			if(nbtag==4){
				double temp=(sf1)*(sf2)*(sf3)*(sf4);
				nPassB[4]+=temp;
				th4[4]->Fill(temp);
			}
			
		}
	}	
	cout<<"entries="<<total<<endl;	
	
	for(int i=0;i<5;i++)cout<<"nPass["<<i<<"]="<<nPassB[i]<<endl;
	
	TH1D * th2o=new TH1D("Nbtagjet","Nbtagjet",5,-0.5,4.5);
	th2o->SetBinContent(1,nPassB[0]);
	th2o->SetBinContent(2,nPassB[1]);
	th2o->SetBinContent(3,nPassB[2]);
	if(nameRoot!=2)th2o->SetBinContent(4,nPassB[3]);
	if(nameRoot!=2)th2o->SetBinContent(5,nPassB[4]);
	
	TH1D * th2s=new TH1D("NbtagjetS","NbtagjetS",5,-0.5,4.5);
	th2s->SetBinContent(1,nPassB[0]);
	th2s->SetBinContent(2,nPassB[1]);
	th2s->SetBinContent(3,nPassB[2]);
	if(nameRoot!=2)th2s->SetBinContent(4,nPassB[3]);
	if(nameRoot!=2)th2s->SetBinContent(5,nPassB[4]);
	
	TFile* outFile = new TFile(Form("root_files_btaggedEff/%s_O.root",st2.data()),"recreate");
	th2o->Write();
	th2s->Scale(2245.87*Xsec/total);
	th2s->Write();
	for(int i=0;i<5;i++)th4[i]->Write();
	outFile->Close();
	
}