
void HH4bBtagMakeEffBase_80(int wMs,int wM, string st,string st2,string option=""){	

	//1=signal ,0=QCD ,2=data
	int nameRoot=1;
	if(st2.find("QCD")!= std::string::npos)nameRoot=0;
	if(st2.find("bGen")!= std::string::npos)nameRoot=0;
	if(st2.find("bEnriched")!= std::string::npos)nameRoot=0;
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
	
	double xBinsForHeavyFlavor[5]={30,140,180,240,3000};
	double xBinsForLightFlavor[11]={20,100,200,300,400,500,600,700,800,900,3000};
	
	TH2D* th1[6];
	th1[0]=new TH2D("effD_b","effD_b",4,xBinsForHeavyFlavor,30,-3,3);
	th1[1]=new TH2D("effN_b","effN_b",4,xBinsForHeavyFlavor,30,-3,3);
	th1[2]=new TH2D("effD_c","effD_c",4,xBinsForHeavyFlavor,30,-3,3);
	th1[3]=new TH2D("effN_c","effN_c",4,xBinsForHeavyFlavor,30,-3,3);
	th1[4]=new TH2D("effD_l","effD_l",4,xBinsForLightFlavor,30,-3,3);
	th1[5]=new TH2D("effN_l","effN_l",4,xBinsForLightFlavor,30,-3,3);
	
	TH1D* th2[6];
	th2[0]=new TH1D("effD_b_1d","effD_b_1d",4,xBinsForHeavyFlavor);
	th2[1]=new TH1D("effN_b_1d","effN_b_1d",4,xBinsForHeavyFlavor);
	th2[2]=new TH1D("effD_c_1d","effD_c_1d",4,xBinsForHeavyFlavor);
	th2[3]=new TH1D("effN_c_1d","effN_c_1d",4,xBinsForHeavyFlavor);
	th2[4]=new TH1D("effD_l_1d","effD_l_1d",10,xBinsForLightFlavor);
	th2[5]=new TH1D("effN_l_1d","effN_l_1d",10,xBinsForLightFlavor);

	for(int i=0;i<6;i++){
			th1[i]->Sumw2();
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
			data.GetEntry(jEntry);
			Int_t nVtx        = data.GetInt("nVtx");
			//0. has a good vertex
			if(nVtx<1)continue;
			nPass[0]++;
			
			//1.trigger
			std::string* trigName = data.GetPtrString("hlt_trigName");
		 	vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
			bool passTrigger=false;
			for(int it=0; it< data.GetPtrStringSize();it++){
				std::string thisTrig= trigName[it];
				bool results = trigResult[it];
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
			if(!passTrigger && nameRoot==2)continue;
			nPass[1]++;

			const int nFATJet=data.GetInt("FATnJet");
			//2.nJets
			if(nFATJet<2)continue;nPass[2]++;
			TClonesArray* fatjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
			float*  FATjetCorrUncUp = data.GetPtrFloat("FATjetCorrUncUp"); 
			float*  FATjetCorrUncDown = data.GetPtrFloat("FATjetCorrUncDown"); 
			TLorentzVector* thisJet ,* thatJet;
			if(JESOption==0){
				thisJet=(TLorentzVector*)fatjetP4->At(0);
				thatJet=(TLorentzVector*)fatjetP4->At(1);
			}
			else if (JESOption==1){
				TLorentzVector thisJetScaleUp= (*((TLorentzVector*)fatjetP4->At(0)))*(1+FATjetCorrUncUp[0]);
				TLorentzVector thatJetScaleUp= (*((TLorentzVector*)fatjetP4->At(1)))*(1+FATjetCorrUncUp[1]);
				thisJet= &thisJetScaleUp;
				thatJet= &thatJetScaleUp;
			}
			else if (JESOption==2){
				TLorentzVector thisJetScaleDown= (*((TLorentzVector*)fatjetP4->At(0)))*(1-FATjetCorrUncDown[0]);
				TLorentzVector thatJetScaleDown= (*((TLorentzVector*)fatjetP4->At(1)))*(1-FATjetCorrUncDown[1]);
				thisJet= &thisJetScaleDown;
				thatJet= &thatJetScaleDown;
			}
			
			//3. Pt 
			if(thisJet->Pt()<200)continue;
			if(thatJet->Pt()<200)continue;
			nPass[3]++;
			//4tightId-----------------------------------------
			vector<bool>    &FATjetPassIDTight = *((vector<bool>*) data.GetPtr("FATjetPassIDTight"));
			if(FATjetPassIDTight[0]==0)continue;
			if(FATjetPassIDTight[1]==0)continue;
			Float_t*  FATjetCEmEF = data.GetPtrFloat("FATjetCEmEF");
			Float_t*  FATjetMuEF = data.GetPtrFloat("FATjetMuEF");
			if(FATjetMuEF[0]>0.8)continue;
			if(FATjetCEmEF[0]>0.9)continue;
			if(FATjetMuEF[1]>0.8)continue;
			if(FATjetCEmEF[1]>0.9)continue;
			nPass[4]++;
			//5. Eta-----------------------------------------
			if(fabs(thisJet->Eta())>2.4)continue;
			if(fabs(thatJet->Eta())>2.4)continue;
			nPass[5]++;
			//6. DEta-----------------------------------------
			float dEta = fabs(thisJet->Eta()-thatJet->Eta());
			if(dEta>1.3)continue;
			nPass[6]++;
			//7. Mjj-----------------------------------------
			float mjjRed = (*thisJet+*thatJet).M()+250-thisJet->M()-thatJet->M();
			if(mjjRed<1000)continue;
			nPass[7]++;
			//8. fatjetPRmassL2L3Corr-----------------------------------------
			Float_t*  fatjetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
			if(fatjetPRmassL2L3Corr[0]<105||fatjetPRmassL2L3Corr[0]>135)continue;
			if(fatjetPRmassL2L3Corr[1]<105||fatjetPRmassL2L3Corr[1]>135)continue;
			nPass[8]++;
			//9.-----------------------------------------
			Float_t*  fatjetTau1 = data.GetPtrFloat("FATjetTau1");
			Float_t*  fatjetTau2 = data.GetPtrFloat("FATjetTau2");
			double tau21_1=(fatjetTau2[0]/fatjetTau1[0]),tau21_2=(fatjetTau2[1]/fatjetTau1[1]);
			//if(tau21_1>0.6 || tau21_2>0.6) continue;
			nPass[9]++;
			
			//10.btag
			vector<float>   *subjetSDCSV =  data.GetPtrVectorFloat("FATsubjetSDCSV");
			vector<float>   *subjetSDPx  =  data.GetPtrVectorFloat("FATsubjetSDPx");
			vector<float>   *subjetSDPy  =  data.GetPtrVectorFloat("FATsubjetSDPy");
			vector<float>   *subjetSDPz  =  data.GetPtrVectorFloat("FATsubjetSDPz");
			vector<float>   *subjetSDE   =  data.GetPtrVectorFloat("FATsubjetSDE");
			vector<Int_t>   *FATsubjetSDHadronFlavor =  data.GetPtrVectorInt("FATsubjetSDHadronFlavor");
			TLorentzVector* thisSub1=new TLorentzVector(0,0,0,0);
			TLorentzVector* thisSub2=new TLorentzVector(0,0,0,0);
			TLorentzVector* thatSub1=new TLorentzVector(0,0,0,0);
			TLorentzVector* thatSub2=new TLorentzVector(0,0,0,0);
			thisSub1->SetPxPyPzE(subjetSDPx[0][0],subjetSDPy[0][0],subjetSDPz[0][0],subjetSDE[0][0]);
			thisSub2->SetPxPyPzE(subjetSDPx[0][1],subjetSDPy[0][1],subjetSDPz[0][1],subjetSDE[0][1]);
			thatSub1->SetPxPyPzE(subjetSDPx[1][0],subjetSDPy[1][0],subjetSDPz[1][0],subjetSDE[1][0]);
			thatSub2->SetPxPyPzE(subjetSDPx[1][1],subjetSDPy[1][1],subjetSDPz[1][1],subjetSDE[1][1]);
			
			if(FATsubjetSDHadronFlavor[0][0]==5)th1[0]->Fill(thisSub1->Pt(),thisSub1->Eta());
			else if(FATsubjetSDHadronFlavor[0][0]==4)th1[2]->Fill(thisSub1->Pt(),thisSub1->Eta());
			else th1[4]->Fill(thisSub1->Pt(),thisSub1->Eta());
			
			if(FATsubjetSDHadronFlavor[0][1]==5)th1[0]->Fill(thisSub2->Pt(),thisSub2->Eta());
			else if(FATsubjetSDHadronFlavor[0][1]==4)th1[2]->Fill(thisSub2->Pt(),thisSub2->Eta());
			else th1[4]->Fill(thisSub2->Pt(),thisSub2->Eta());
			
			if(FATsubjetSDHadronFlavor[1][0]==5)th1[0]->Fill(thatSub1->Pt(),thatSub1->Eta());
			else if(FATsubjetSDHadronFlavor[1][0]==4)th1[2]->Fill(thatSub1->Pt(),thatSub1->Eta());
			else th1[4]->Fill(thatSub1->Pt(),thatSub1->Eta());
			
			if(FATsubjetSDHadronFlavor[1][1]==5)th1[0]->Fill(thatSub2->Pt(),thatSub2->Eta());
			else if(FATsubjetSDHadronFlavor[1][1]==4)th1[2]->Fill(thatSub2->Pt(),thatSub2->Eta());
			else th1[4]->Fill(thatSub2->Pt(),thatSub2->Eta());
		
			if(FATsubjetSDHadronFlavor[0][0]==5)th2[0]->Fill(thisSub1->Pt());
			else if(FATsubjetSDHadronFlavor[0][0]==4)th2[2]->Fill(thisSub1->Pt());
			else th2[4]->Fill(thisSub1->Pt());
				
			if(FATsubjetSDHadronFlavor[0][1]==5)th2[0]->Fill(thisSub2->Pt());
			else if(FATsubjetSDHadronFlavor[0][1]==4)th2[2]->Fill(thisSub2->Pt());
			else th2[4]->Fill(thisSub2->Pt());
			
			if(FATsubjetSDHadronFlavor[1][0]==5)th2[0]->Fill(thatSub1->Pt());
			else if(FATsubjetSDHadronFlavor[1][0]==4)th2[2]->Fill(thatSub1->Pt());
			else th2[4]->Fill(thatSub1->Pt());
			
			if(FATsubjetSDHadronFlavor[1][1]==5)th2[0]->Fill(thatSub2->Pt());
			else if(FATsubjetSDHadronFlavor[1][1]==4)th2[2]->Fill(thatSub2->Pt());
			else th2[4]->Fill(thatSub2->Pt());
			
			if(subjetSDCSV[0][0]>0.46){
				if(FATsubjetSDHadronFlavor[0][0]==5)th1[1]->Fill(thisSub1->Pt(),thisSub1->Eta());
				else if(FATsubjetSDHadronFlavor[0][0]==4)th1[3]->Fill(thisSub1->Pt(),thisSub1->Eta());
				else th1[5]->Fill(thisSub1->Pt(),thisSub1->Eta());
				
				if(FATsubjetSDHadronFlavor[0][0]==5)th2[1]->Fill(thisSub1->Pt());
				else if(FATsubjetSDHadronFlavor[0][0]==4)th2[3]->Fill(thisSub1->Pt());
				else th2[5]->Fill(thisSub1->Pt());
			}
			if(subjetSDCSV[0][1]>0.46){
				if(FATsubjetSDHadronFlavor[0][1]==5)th1[1]->Fill(thisSub2->Pt(),thisSub2->Eta());
				else if(FATsubjetSDHadronFlavor[0][1]==4)th1[3]->Fill(thisSub2->Pt(),thisSub2->Eta());
				else th1[5]->Fill(thisSub2->Pt(),thisSub2->Eta());
				
				if(FATsubjetSDHadronFlavor[0][1]==5)th2[1]->Fill(thisSub2->Pt());
				else if(FATsubjetSDHadronFlavor[0][1]==4)th2[3]->Fill(thisSub2->Pt());
				else th2[5]->Fill(thisSub2->Pt());
			}
			if(subjetSDCSV[1][0]>0.46){
				if(FATsubjetSDHadronFlavor[1][0]==5)th1[1]->Fill(thatSub1->Pt(),thatSub1->Eta());
				else if(FATsubjetSDHadronFlavor[1][0]==4)th1[3]->Fill(thatSub1->Pt(),thatSub1->Eta());
				else th1[5]->Fill(thatSub1->Pt(),thatSub1->Eta());
				
				if(FATsubjetSDHadronFlavor[1][0]==5)th2[1]->Fill(thatSub1->Pt());
				else if(FATsubjetSDHadronFlavor[1][0]==4)th2[3]->Fill(thatSub1->Pt());
				else th2[5]->Fill(thatSub1->Pt());
			}
			if(subjetSDCSV[1][1]>0.46){
				if(FATsubjetSDHadronFlavor[1][1]==5)th1[1]->Fill(thatSub2->Pt(),thatSub2->Eta());
				else if(FATsubjetSDHadronFlavor[1][1]==4)th1[3]->Fill(thatSub2->Pt(),thatSub2->Eta());
				else th1[5]->Fill(thatSub2->Pt(),thatSub2->Eta());
				
				if(FATsubjetSDHadronFlavor[1][1]==5)th2[1]->Fill(thatSub2->Pt());
				else if(FATsubjetSDHadronFlavor[1][1]==4)th2[3]->Fill(thatSub2->Pt());
				else th2[5]->Fill(thatSub2->Pt());
			}
		}
	}	
	cout<<"entries="<<total<<endl;	
	for(int i=0;i<9;i++)cout<<"nPass["<<i<<"]="<<nPass[i]<<endl;
	
	TFile* outFile ;
	if(JESOption==0)outFile= new TFile(Form("btagEffSource/%s.root",st2.data()),"recreate");
	else if(JESOption==1)outFile= new TFile(Form("btagEffSource/%s_JESUp.root",st2.data()),"recreate");
	else if(JESOption==2)outFile= new TFile(Form("btagEffSource/%s_JESDown.root",st2.data()),"recreate");
	
	for (int i=0;i<6;i++)th1[i]->Write();
	for (int i=0;i<6;i++)th2[i]->Write();
	outFile->Close();
}
