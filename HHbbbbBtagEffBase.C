void HHbbbbBtagEffBase(int wMs,int wM, string st,string st2,double Xsec,int nameRoot=0){	
	TFile *f1;
	if(nameRoot==2)f1=TFile::Open("root_files_btaggedEff/data.root");
	else f1=TFile::Open(Form("root_files_btaggedEff/%s.root",st2.data()));
	TH2D* th1[6];
	th1[0]=(TH2D*)f1->FindObjectAny("effD_b");
	th1[1]=(TH2D*)f1->FindObjectAny("effN_b");
	th1[2]=(TH2D*)f1->FindObjectAny("effD_c");
	th1[3]=(TH2D*)f1->FindObjectAny("effN_c");
	th1[4]=(TH2D*)f1->FindObjectAny("effD_l");
	th1[5]=(TH2D*)f1->FindObjectAny("effN_l");
	
	th1[1]->Divide(th1[0]);
	th1[3]->Divide(th1[2]);
	th1[5]->Divide(th1[4]);
	for(int i=0;i<6;i++)th1[i]->Sumw2();
	
	TFile *f2;
	f2=TFile::Open("root_files_btaggedEff/data.root");
	TH2D* th2[6];
	th2[0]=(TH2D*)f1->FindObjectAny("effD_b");
	th2[1]=(TH2D*)f1->FindObjectAny("effN_b");
	th2[2]=(TH2D*)f1->FindObjectAny("effD_c");
	th2[3]=(TH2D*)f1->FindObjectAny("effN_c");
	th2[4]=(TH2D*)f1->FindObjectAny("effD_l");
	th2[5]=(TH2D*)f1->FindObjectAny("effN_l");
	
	th2[1]->Divide(th2[0]);
	th2[3]->Divide(th2[2]);
	th2[5]->Divide(th2[4]);
	for(int i=0;i<6;i++)th2[i]->Sumw2();
	
	TH2D* th3[6];
	
	th3[0]=new TH2D("zeroSF_b","zeroSF_b",200,0,2000,60,-3,3);
	th3[1]=new TH2D("zeroSF_c","zeroSF_c",200,0,2000,60,-3,3);
	th3[2]=new TH2D("zeroSF_l","zeroSF_l",200,0,2000,60,-3,3);
	
	th3[3]=new TH2D("SF_vs_Pt_b","SF_vs_Pt_b",120,0,1200,40,0.8,1.2);
	th3[4]=new TH2D("SF_vs_Pt_c","SF_vs_Pt_c",120,0,1200,40,0.8,1.2);
	th3[5]=new TH2D("SF_vs_Pt_l","SF_vs_Pt_l",120,0,1200,40,0.8,1.2);

	for(int i=0;i<6;i++)th3[i]->Sumw2();
	TH1D* th4[14];
	
	th4[0]=new TH1D("SF_jet0_sub0_pass","SF_jet0_sub0_pass",40,0.8,1.2);
	th4[1]=new TH1D("SF_jet0_sub1_pass","SF_jet0_sub1_pass",40,0.8,1.2);
	th4[2]=new TH1D("SF_jet1_sub0_pass","SF_jet1_sub0_pass",40,0.8,1.2);
	th4[3]=new TH1D("SF_jet1_sub1_pass","SF_jet1_sub1_pass",40,0.8,1.2);
	th4[4]=new TH1D("SF_jet0_sub0_fail","SF_jet0_sub0_fail",40,0.8,1.2);
	th4[5]=new TH1D("SF_jet0_sub1_fail","SF_jet0_sub1_fail",40,0.8,1.2);
	th4[6]=new TH1D("SF_jet1_sub0_fail","SF_jet1_sub0_fail",40,0.8,1.2);
	th4[7]=new TH1D("SF_jet1_sub1_fail","SF_jet1_sub1_fail",40,0.8,1.2);
	th4[8]=new TH1D("weight","weight",40,0,2);
	th4[9]=new TH1D("weight_0b","weight_0b",40,0,2);
	th4[10]=new TH1D("weight_1b","weight_1b",40,0,2);
	th4[11]=new TH1D("weight_2b","weight_2b",40,0,2);
	th4[12]=new TH1D("weight_2b_ll","weight_2b_ll",40,0,2);
	th4[13]=new TH1D("weight_2b_onel","weight_2b_onel",40,0,2);
	for(int i=0;i<14;i++)th4[i]->Sumw2();
	
	TH1D * th5[51];
	
	th5[0]=new TH1D("Pt_j0_sj0_0b","Pt_j0_sj0_0b",200,0,2000);
	th5[1]=new TH1D("Pt_j0_sj1_0b","Pt_j0_sj1_0b",200,0,2000);
	th5[2]=new TH1D("Pt_j1_sj0_0b","Pt_j1_sj0_0b",200,0,2000);
	th5[3]=new TH1D("Pt_j1_sj1_0b","Pt_j1_sj1_0b",200,0,2000);
	
	th5[4]=new TH1D("Pt_j0_sj0_1b","Pt_j0_sj0_1b",200,0,2000);
	th5[5]=new TH1D("Pt_j0_sj1_1b","Pt_j0_sj1_1b",200,0,2000);
	th5[6]=new TH1D("Pt_j1_sj0_1b","Pt_j1_sj0_1b",200,0,2000);
	th5[7]=new TH1D("Pt_j1_sj1_1b","Pt_j1_sj1_1b",200,0,2000);
	
	th5[8]=new TH1D("Pt_j0_sj0_2b","Pt_j0_sj0_2b",200,0,2000);
	th5[9]=new TH1D("Pt_j0_sj1_2b","Pt_j0_sj1_2b",200,0,2000);
	th5[10]=new TH1D("Pt_j1_sj0_2b","Pt_j1_sj0_2b",200,0,2000);
	th5[11]=new TH1D("Pt_j1_sj1_2b","Pt_j1_sj1_2b",200,0,2000);
	
	th5[12]=new TH1D("deltaR0_0b","deltaR0_0b",20,0,1);
	th5[13]=new TH1D("deltaR1_0b","deltaR1_0b",20,0,1);
	th5[14]=new TH1D("deltaR0_1b","deltaR0_1b",20,0,1);
	th5[15]=new TH1D("deltaR1_1b","deltaR1_1b",20,0,1);
	th5[16]=new TH1D("deltaR0_2b","deltaR0_2b",20,0,1);
	th5[17]=new TH1D("deltaR1_2b","deltaR1_2b",20,0,1);
	
	th5[18]=new TH1D("Eta_j0_sj0_0b","Eta_j0_sj0_0b",60,-3,3);
	th5[19]=new TH1D("Eta_j0_sj1_0b","Eta_j0_sj1_0b",60,-3,3);
	th5[20]=new TH1D("Eta_j1_sj0_0b","Eta_j1_sj0_0b",60,-3,3);
	th5[21]=new TH1D("Eta_j1_sj1_0b","Eta_j1_sj1_0b",60,-3,3);
	
	th5[22]=new TH1D("Eta_j0_sj0_1b","Eta_j0_sj0_1b",60,-3,3);
	th5[23]=new TH1D("Eta_j0_sj1_1b","Eta_j0_sj1_1b",60,-3,3);
	th5[24]=new TH1D("Eta_j1_sj0_1b","Eta_j1_sj0_1b",60,-3,3);
	th5[25]=new TH1D("Eta_j1_sj1_1b","Eta_j1_sj1_1b",60,-3,3);
	
	th5[26]=new TH1D("Eta_j0_sj0_2b","Eta_j0_sj0_2b",60,-3,3);
	th5[27]=new TH1D("Eta_j0_sj1_2b","Eta_j0_sj1_2b",60,-3,3);
	th5[28]=new TH1D("Eta_j1_sj0_2b","Eta_j1_sj0_2b",60,-3,3);
	th5[29]=new TH1D("Eta_j1_sj1_2b","Eta_j1_sj1_2b",60,-3,3);
	
	th5[30]=new TH1D("Pt_j0_0b","Pt_j0_0b",200,0,2000);
	th5[31]=new TH1D("Pt_j1_0b","Pt_j1_0b",200,0,2000);
	th5[32]=new TH1D("Pt_j0_1b","Pt_j0_1b",200,0,2000);
	th5[33]=new TH1D("Pt_j1_1b","Pt_j1_1b",200,0,2000);
	th5[34]=new TH1D("Pt_j0_2b","Pt_j0_2b",200,0,2000);
	th5[35]=new TH1D("Pt_j1_2b","Pt_j1_2b",200,0,2000);
	
	th5[36]=new TH1D("Eta_j0_0b","Eta_j0_0b",60,-3,3);
	th5[37]=new TH1D("Eta_j1_0b","Eta_j1_0b",60,-3,3);
	th5[38]=new TH1D("Eta_j0_1b","Eta_j0_1b",60,-3,3);
	th5[39]=new TH1D("Eta_j1_1b","Eta_j1_1b",60,-3,3);
	th5[40]=new TH1D("Eta_j0_2b","Eta_j0_2b",60,-3,3);
	th5[41]=new TH1D("Eta_j1_2b","Eta_j1_2b",60,-3,3);
	
	th5[42]=new TH1D("totalMass_0b","totalMass_0b",75,1000,4000);
	th5[43]=new TH1D("totalMass_1b","totalMass_1b",75,1000,4000);
	th5[44]=new TH1D("totalMass_2b","totalMass_2b",75,1000,4000);
	
	th5[45]=new TH1D("prMass_j0_0b","prMass_j0_0b",15,90,150);
	th5[46]=new TH1D("prMass_j1_0b","prMass_j1_0b",15,90,150);
	th5[47]=new TH1D("prMass_j0_1b","prMass_j0_1b",15,90,150);
	th5[48]=new TH1D("prMass_j1_1b","prMass_j1_1b",15,90,150);
	th5[49]=new TH1D("prMass_j0_2b","prMass_j0_2b",15,90,150);
	th5[50]=new TH1D("prMass_j1_2b","prMass_j1_2b",15,90,150);
	
	TH1D * th6[51];
	
	th6[0]=new TH1D("Pt_j0_sj0_0bs","Pt_j0_sj0_0b",200,0,2000);
	th6[1]=new TH1D("Pt_j0_sj1_0bs","Pt_j0_sj1_0b",200,0,2000);
	th6[2]=new TH1D("Pt_j1_sj0_0bs","Pt_j1_sj0_0b",200,0,2000);
	th6[3]=new TH1D("Pt_j1_sj1_0bs","Pt_j1_sj1_0b",200,0,2000);
	
	th6[4]=new TH1D("Pt_j0_sj0_1bs","Pt_j0_sj0_1b",200,0,2000);
	th6[5]=new TH1D("Pt_j0_sj1_1bs","Pt_j0_sj1_1b",200,0,2000);
	th6[6]=new TH1D("Pt_j1_sj0_1bs","Pt_j1_sj0_1b",200,0,2000);
	th6[7]=new TH1D("Pt_j1_sj1_1bs","Pt_j1_sj1_1b",200,0,2000);
	
	th6[8]=new TH1D("Pt_j0_sj0_2bs","Pt_j0_sj0_2b",200,0,2000);
	th6[9]=new TH1D("Pt_j0_sj1_2bs","Pt_j0_sj1_2b",200,0,2000);
	th6[10]=new TH1D("Pt_j1_sj0_2bs","Pt_j1_sj0_2b",200,0,2000);
	th6[11]=new TH1D("Pt_j1_sj1_2bs","Pt_j1_sj1_2b",200,0,2000);
	
	th6[12]=new TH1D("deltaR0_0bs","deltaR0_0b",20,0,1);
	th6[13]=new TH1D("deltaR1_0bs","deltaR1_0b",20,0,1);
	th6[14]=new TH1D("deltaR0_1bs","deltaR0_1b",20,0,1);
	th6[15]=new TH1D("deltaR1_1bs","deltaR1_1b",20,0,1);
	th6[16]=new TH1D("deltaR0_2bs","deltaR0_2b",20,0,1);
	th6[17]=new TH1D("deltaR1_2bs","deltaR1_2b",20,0,1);
	
	th6[18]=new TH1D("Eta_j0_sj0_0bs","Eta_j0_sj0_0b",60,-3,3);
	th6[19]=new TH1D("Eta_j0_sj1_0bs","Eta_j0_sj1_0b",60,-3,3);
	th6[20]=new TH1D("Eta_j1_sj0_0bs","Eta_j1_sj0_0b",60,-3,3);
	th6[21]=new TH1D("Eta_j1_sj1_0bs","Eta_j1_sj1_0b",60,-3,3);
	
	th6[22]=new TH1D("Eta_j0_sj0_1bs","Eta_j0_sj0_1b",60,-3,3);
	th6[23]=new TH1D("Eta_j0_sj1_1bs","Eta_j0_sj1_1b",60,-3,3);
	th6[24]=new TH1D("Eta_j1_sj0_1bs","Eta_j1_sj0_1b",60,-3,3);
	th6[25]=new TH1D("Eta_j1_sj1_1bs","Eta_j1_sj1_1b",60,-3,3);
	
	th6[26]=new TH1D("Eta_j0_sj0_2bs","Eta_j0_sj0_2b",60,-3,3);
	th6[27]=new TH1D("Eta_j0_sj1_2bs","Eta_j0_sj1_2b",60,-3,3);
	th6[28]=new TH1D("Eta_j1_sj0_2bs","Eta_j1_sj0_2b",60,-3,3);
	th6[29]=new TH1D("Eta_j1_sj1_2bs","Eta_j1_sj1_2b",60,-3,3);
	
	th6[30]=new TH1D("Pt_j0_0bs","Pt_j0_0b",200,0,2000);
	th6[31]=new TH1D("Pt_j1_0bs","Pt_j1_0b",200,0,2000);
	th6[32]=new TH1D("Pt_j0_1bs","Pt_j0_1b",200,0,2000);
	th6[33]=new TH1D("Pt_j1_1bs","Pt_j1_1b",200,0,2000);
	th6[34]=new TH1D("Pt_j0_2bs","Pt_j0_2b",200,0,2000);
	th6[35]=new TH1D("Pt_j1_2bs","Pt_j1_2b",200,0,2000);
	
	th6[36]=new TH1D("Eta_j0_0bs","Eta_j0_0b",60,-3,3);
	th6[37]=new TH1D("Eta_j1_0bs","Eta_j1_0b",60,-3,3);
	th6[38]=new TH1D("Eta_j0_1bs","Eta_j0_1b",60,-3,3);
	th6[39]=new TH1D("Eta_j1_1bs","Eta_j1_1b",60,-3,3);
	th6[40]=new TH1D("Eta_j0_2bs","Eta_j0_2b",60,-3,3);
	th6[41]=new TH1D("Eta_j1_2bs","Eta_j1_2b",60,-3,3);
	
	th6[42]=new TH1D("totalMass_0bs","totalMass_0b",75,1000,4000);
	th6[43]=new TH1D("totalMass_1bs","totalMass_1b",75,1000,4000);
	th6[44]=new TH1D("totalMass_2bs","totalMass_2b",75,1000,4000);
	
	th6[45]=new TH1D("prMass_j0_0bs","prMass_j0_0b",15,90,150);
	th6[46]=new TH1D("prMass_j1_0bs","prMass_j1_0b",15,90,150);
	th6[47]=new TH1D("prMass_j0_1bs","prMass_j0_1b",15,90,150);
	th6[48]=new TH1D("prMass_j1_1bs","prMass_j1_1b",15,90,150);
	th6[49]=new TH1D("prMass_j0_2bs","prMass_j0_2b",15,90,150);
	th6[50]=new TH1D("prMass_j1_2bs","prMass_j1_2b",15,90,150);
	
	for(int i=0;i<51;i++)th5[i]->Sumw2();
	for(int i=0;i<51;i++)th6[i]->Sumw2();
	
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
			
			if(FATsubjetSDHadronFlavor[0][0]==5){
				sf1=HF.eval(BTagEntry::FLAV_B,thisSub1->Eta(),pt1); 
				th3[3]->Fill(pt1,sf1);
			}
			else if(FATsubjetSDHadronFlavor[0][0]==4){
				sf1=HF.eval(BTagEntry::FLAV_C,thisSub1->Eta(),pt1); 
				th3[4]->Fill(pt1,sf1);
			}
			else {
				sf1=LF.eval(BTagEntry::FLAV_UDSG,thisSub1->Eta(),pt1); 
				th3[5]->Fill(pt1,sf1);
			}
			
			if(FATsubjetSDHadronFlavor[0][1]==5){
				sf2=HF.eval(BTagEntry::FLAV_B,thisSub2->Eta(),pt2); 
				th3[3]->Fill(pt2,sf2);
			}
			else if(FATsubjetSDHadronFlavor[0][1]==4){
				sf2=HF.eval(BTagEntry::FLAV_C,thisSub2->Eta(),pt2); 
				th3[4]->Fill(pt2,sf2);
			}
			else{
				sf2=LF.eval(BTagEntry::FLAV_UDSG,thisSub2->Eta(),pt2); 
				th3[5]->Fill(pt2,sf2);
			} 
			
			if(FATsubjetSDHadronFlavor[1][0]==5){
				sf3=HF.eval(BTagEntry::FLAV_B,thatSub1->Eta(),pt3); 
				th3[3]->Fill(pt3,sf3);
			}
			else if(FATsubjetSDHadronFlavor[1][0]==4){
				sf3=HF.eval(BTagEntry::FLAV_C,thatSub1->Eta(),pt3); 
				th3[4]->Fill(pt3,sf3);
			}
			else {
				sf3=LF.eval(BTagEntry::FLAV_UDSG,thatSub1->Eta(),pt3); 
				th3[5]->Fill(pt3,sf3);
			}
			
			if(FATsubjetSDHadronFlavor[1][1]==5){
				sf4=HF.eval(BTagEntry::FLAV_B,thatSub2->Eta(),pt4); 
				th3[3]->Fill(pt4,sf4);
			}
			else if(FATsubjetSDHadronFlavor[1][1]==4){
				sf4=HF.eval(BTagEntry::FLAV_C,thatSub2->Eta(),pt4); 
				th3[4]->Fill(pt4,sf4);
			}
			else {
				sf4=LF.eval(BTagEntry::FLAV_UDSG,thatSub2->Eta(),pt4); 
				th3[5]->Fill(pt4,sf4);
			}
			
			if(sf1==0 && FATsubjetSDHadronFlavor[0][0]==5) cout<<"pt1="<<pt1<<endl;
			
			if(sf1==0 && FATsubjetSDHadronFlavor[0][0]==5 ) th3[0]->Fill(pt1,thisSub1->Eta());
			if(sf1==0 && FATsubjetSDHadronFlavor[0][0]==4 ) th3[1]->Fill(pt1,thisSub1->Eta());
			if(sf1==0 && FATsubjetSDHadronFlavor[0][0]!=5 && FATsubjetSDHadronFlavor[0][0]!=4 ) th3[2]->Fill(pt1,thisSub1->Eta());
			
			if(sf2==0 && FATsubjetSDHadronFlavor[0][1]==5 ) th3[0]->Fill(pt2,thisSub2->Eta());
			if(sf2==0 && FATsubjetSDHadronFlavor[0][1]==4 ) th3[1]->Fill(pt2,thisSub2->Eta());
			if(sf2==0 && FATsubjetSDHadronFlavor[0][1]!=5 && FATsubjetSDHadronFlavor[0][1]!=4 ) th3[2]->Fill(pt2,thisSub2->Eta());
			
			if(sf3==0 && FATsubjetSDHadronFlavor[1][0]==5 ) th3[0]->Fill(pt3,thatSub1->Eta());
			if(sf3==0 && FATsubjetSDHadronFlavor[1][0]==4 ) th3[1]->Fill(pt3,thatSub1->Eta());
			if(sf3==0 && FATsubjetSDHadronFlavor[1][0]!=5 && FATsubjetSDHadronFlavor[1][0]!=4 ) th3[2]->Fill(pt3,thatSub1->Eta());
			
			if(sf4==0 && FATsubjetSDHadronFlavor[1][1]==5 ) th3[0]->Fill(pt4,thatSub2->Eta());
			if(sf4==0 && FATsubjetSDHadronFlavor[1][1]==4 ) th3[1]->Fill(pt4,thatSub2->Eta());
			if(sf4==0 && FATsubjetSDHadronFlavor[1][1]!=5 && FATsubjetSDHadronFlavor[1][1]!=4 ) th3[2]->Fill(pt4,thatSub2->Eta());
						
			double eff1,eff2,eff3,eff4;

			if(FATsubjetSDHadronFlavor[0][0]==5)eff1=th1[1]->GetBinContent(ceil(thisSub1->Pt()/10),ceil(thisSub1->Eta()/0.1)+30);
			else if(FATsubjetSDHadronFlavor[0][0]==4)eff1=th1[3]->GetBinContent(ceil(thisSub1->Pt()/10),ceil(thisSub1->Eta()/0.1)+30);
			else eff1=th1[5]->GetBinContent(ceil(thisSub1->Pt()/10),ceil(thisSub1->Eta()/0.1)+30);
			
			if(FATsubjetSDHadronFlavor[0][1]==5)eff2=th1[1]->GetBinContent(ceil(thisSub2->Pt()/10),ceil(thisSub2->Eta()/0.1)+30);
			else if(FATsubjetSDHadronFlavor[0][1]==4)eff2=th1[3]->GetBinContent(ceil(thisSub2->Pt()/10),ceil(thisSub2->Eta()/0.1)+30);
			else eff2=th1[5]->GetBinContent(ceil(thisSub2->Pt()/10),ceil(thisSub2->Eta()/0.1)+30);
			
			if(FATsubjetSDHadronFlavor[1][0]==5)eff3=th1[1]->GetBinContent(ceil(thatSub1->Pt()/10),ceil(thatSub1->Eta()/0.1)+30);
			else if(FATsubjetSDHadronFlavor[1][0]==4)eff3=th1[3]->GetBinContent(ceil(thatSub1->Pt()/10),ceil(thatSub1->Eta()/0.1)+30);
			else eff3=th1[5]->GetBinContent(ceil(thatSub1->Pt()/10),ceil(thatSub1->Eta()/0.1)+30);
			
			if(FATsubjetSDHadronFlavor[1][1]==5)eff4=th1[1]->GetBinContent(ceil(thatSub2->Pt()/10),ceil(thatSub2->Eta()/0.1)+30);
			else if(FATsubjetSDHadronFlavor[1][1]==4)eff4=th1[3]->GetBinContent(ceil(thatSub2->Pt()/10),ceil(thatSub2->Eta()/0.1)+30);
			else eff4=th1[5]->GetBinContent(ceil(thatSub2->Pt()/10),ceil(thatSub2->Eta()/0.1)+30);
			
			int nbtag=0;
			if(subjetSDCSV[0][0]>0.605)nbtag++;
			if(subjetSDCSV[0][1]>0.605)nbtag++;
			if(subjetSDCSV[1][0]>0.605)nbtag++;
			if(subjetSDCSV[1][1]>0.605)nbtag++;
			
			double scaleFactor=1;
			if(subjetSDCSV[0][0]>=0.605)scaleFactor*=sf1;
			else scaleFactor*=((1-eff1*sf1)/(1-eff1));
			//if(nbtag==2)cout<<jEntry<<",SF="<<scaleFactor<<",CSV="<<subjetSDCSV[0][0]<<",eff="<<eff1<<",SD="<<FATsubjetSDHadronFlavor[0][0]<<",sf="<<sf1<<endl;
			if(subjetSDCSV[0][1]>=0.605)scaleFactor*=sf2;
			else scaleFactor*=((1-eff2*sf2)/(1-eff2));
			//if(nbtag==2)cout<<jEntry<<",SF="<<scaleFactor<<",CSV="<<subjetSDCSV[0][1]<<",eff="<<eff2<<",SD="<<FATsubjetSDHadronFlavor[0][1]<<",sf="<<sf2<<endl;
			if(subjetSDCSV[1][0]>=0.605)scaleFactor*=sf3;
			else scaleFactor*=((1-eff3*sf3)/(1-eff3));
			//if(nbtag==2)cout<<jEntry<<",SF="<<scaleFactor<<",CSV="<<subjetSDCSV[1][0]<<",eff="<<eff3<<",SD="<<FATsubjetSDHadronFlavor[1][0]<<",sf="<<sf3<<endl;
			if(subjetSDCSV[1][1]>=0.605)scaleFactor*=sf4;
			else scaleFactor*=((1-eff4*sf4)/(1-eff4));
			//if(nbtag==2)cout<<jEntry<<",SF="<<scaleFactor<<",CSV="<<subjetSDCSV[1][1]<<",eff="<<eff4<<",SD="<<FATsubjetSDHadronFlavor[1][1]<<",sf="<<sf4<<endl;
			
			if(scaleFactor>1.3 && (nbtag==4||nbtag==3) ){
				cout<<jEntry<<",SF="<<scaleFactor<<",CSV="<<subjetSDCSV[0][0]<<",eff="<<eff1<<",SD="<<FATsubjetSDHadronFlavor[0][0]<<",sf="<<sf1<<endl;
				cout<<jEntry<<",SF="<<scaleFactor<<",CSV="<<subjetSDCSV[0][1]<<",eff="<<eff2<<",SD="<<FATsubjetSDHadronFlavor[0][1]<<",sf="<<sf2<<endl;
				cout<<jEntry<<",SF="<<scaleFactor<<",CSV="<<subjetSDCSV[1][0]<<",eff="<<eff3<<",SD="<<FATsubjetSDHadronFlavor[1][0]<<",sf="<<sf3<<endl;
				cout<<jEntry<<",SF="<<scaleFactor<<",CSV="<<subjetSDCSV[1][1]<<",eff="<<eff4<<",SD="<<FATsubjetSDHadronFlavor[1][1]<<",sf="<<sf4<<endl;
				
				
			}
			
			int nbtag2=0;
			if(subjetSDCSV[0][0]>0.605 &&FATsubjetSDHadronFlavor[0][0]==0)nbtag2++;
			if(subjetSDCSV[0][1]>0.605 &&FATsubjetSDHadronFlavor[0][1]==0)nbtag2++;
			if(subjetSDCSV[1][0]>0.605 &&FATsubjetSDHadronFlavor[1][0]==0)nbtag2++;
			if(subjetSDCSV[1][1]>0.605 &&FATsubjetSDHadronFlavor[1][1]==0)nbtag2++;
			
			
			
			if(subjetSDCSV[0][0]>0.605)th4[0]->Fill(sf1);
			else th4[4]->Fill(sf1);
			if(subjetSDCSV[0][1]>0.605)th4[1]->Fill(sf1);
			else th4[5]->Fill(sf1);
			if(subjetSDCSV[1][0]>0.605)th4[2]->Fill(sf1);
			else th4[6]->Fill(sf1);
			if(subjetSDCSV[1][1]>0.605)th4[3]->Fill(sf1);
			else th4[7]->Fill(sf1);
			th4[8]->Fill(scaleFactor);
			
			if(nbtag2==2 && nbtag==2)th4[12]->Fill(scaleFactor);
			if(nbtag2==1 && nbtag==2)th4[13]->Fill(scaleFactor);
			
			
			double dr1=thisSub1->DeltaR(*thisSub2),dr2=thatSub1->DeltaR(*thatSub2);
			
			if(nbtag==0){
				th4[9]->Fill(scaleFactor);
				nPassB[0]+=scaleFactor;
				th5[0]->Fill(thisSub1->Pt());
				th5[1]->Fill(thisSub2->Pt());
				th5[2]->Fill(thatSub1->Pt());
				th5[3]->Fill(thatSub2->Pt());
				th5[18]->Fill(thisSub1->Eta());
				th5[19]->Fill(thisSub2->Eta());
				th5[20]->Fill(thatSub1->Eta());
				th5[21]->Fill(thatSub2->Eta());
				th5[12]->Fill(dr1);
				th5[13]->Fill(dr2);
				th5[30]->Fill((*thisSub1+*thisSub2).Pt());
				th5[31]->Fill((*thatSub1+*thatSub2).Pt());
				th5[36]->Fill((*thisSub1+*thisSub2).Eta());
				th5[37]->Fill((*thatSub1+*thatSub2).Eta());
				th5[42]->Fill((*thisJet+*thatJet).M());
				th5[45]->Fill(fatjetPRmassL2L3Corr[0]);
				th5[46]->Fill(fatjetPRmassL2L3Corr[1]);
				
				th6[0]->Fill(thisSub1->Pt(),scaleFactor);
				th6[1]->Fill(thisSub2->Pt(),scaleFactor);
				th6[2]->Fill(thatSub1->Pt(),scaleFactor);
				th6[3]->Fill(thatSub2->Pt(),scaleFactor);
				th6[18]->Fill(thisSub1->Eta(),scaleFactor);
				th6[19]->Fill(thisSub2->Eta(),scaleFactor);
				th6[20]->Fill(thatSub1->Eta(),scaleFactor);
				th6[21]->Fill(thatSub2->Eta(),scaleFactor);
				th6[12]->Fill(dr1,scaleFactor);
				th6[13]->Fill(dr2,scaleFactor);
				th6[30]->Fill((*thisSub1+*thisSub2).Pt(),scaleFactor);
				th6[31]->Fill((*thatSub1+*thatSub2).Pt(),scaleFactor);
				th6[36]->Fill((*thisSub1+*thisSub2).Eta(),scaleFactor);
				th6[37]->Fill((*thatSub1+*thatSub2).Eta(),scaleFactor);
				th6[42]->Fill((*thisJet+*thatJet).M(),scaleFactor);
				th6[45]->Fill(fatjetPRmassL2L3Corr[0],scaleFactor);
				th6[46]->Fill(fatjetPRmassL2L3Corr[1],scaleFactor);
			}
			if(nbtag==1){
				th4[10]->Fill(scaleFactor);
				nPassB[1]+=scaleFactor;
				th5[4]->Fill(thisSub1->Pt());
				th5[5]->Fill(thisSub2->Pt());
				th5[6]->Fill(thatSub1->Pt());
				th5[7]->Fill(thatSub2->Pt());
				th5[22]->Fill(thisSub1->Eta());
				th5[23]->Fill(thisSub2->Eta());
				th5[24]->Fill(thatSub1->Eta());
				th5[25]->Fill(thatSub2->Eta());
				th5[14]->Fill(dr1);
				th5[15]->Fill(dr2);
				th5[32]->Fill((*thisSub1+*thisSub2).Pt());
				th5[33]->Fill((*thatSub1+*thatSub2).Pt());
				th5[38]->Fill((*thisSub1+*thisSub2).Eta());
				th5[39]->Fill((*thatSub1+*thatSub2).Eta());
				th5[43]->Fill((*thisJet+*thatJet).M());
				th5[47]->Fill(fatjetPRmassL2L3Corr[0]);
				th5[48]->Fill(fatjetPRmassL2L3Corr[1]);
				
				th6[4]->Fill(thisSub1->Pt(),scaleFactor);
				th6[5]->Fill(thisSub2->Pt(),scaleFactor);
				th6[6]->Fill(thatSub1->Pt(),scaleFactor);
				th6[7]->Fill(thatSub2->Pt(),scaleFactor);
				th6[22]->Fill(thisSub1->Eta(),scaleFactor);
				th6[23]->Fill(thisSub2->Eta(),scaleFactor);
				th6[24]->Fill(thatSub1->Eta(),scaleFactor);
				th6[25]->Fill(thatSub2->Eta(),scaleFactor);
				th6[14]->Fill(dr1,scaleFactor);
				th6[15]->Fill(dr2,scaleFactor);
				th6[32]->Fill((*thisSub1+*thisSub2).Pt(),scaleFactor);
				th6[33]->Fill((*thatSub1+*thatSub2).Pt(),scaleFactor);
				th6[38]->Fill((*thisSub1+*thisSub2).Eta(),scaleFactor);
				th6[39]->Fill((*thatSub1+*thatSub2).Eta(),scaleFactor);
				th6[43]->Fill((*thisJet+*thatJet).M(),scaleFactor);
				th6[47]->Fill(fatjetPRmassL2L3Corr[0],scaleFactor);
				th6[48]->Fill(fatjetPRmassL2L3Corr[1],scaleFactor);
			}
			if(nbtag==2){
				th4[11]->Fill(scaleFactor);
				nPassB[2]+=scaleFactor;
				th5[8]->Fill(thisSub1->Pt());
				th5[9]->Fill(thisSub2->Pt());
				th5[10]->Fill(thatSub1->Pt());
				th5[11]->Fill(thatSub2->Pt());
				th5[26]->Fill(thisSub1->Eta());
				th5[27]->Fill(thisSub2->Eta());
				th5[28]->Fill(thatSub1->Eta());
				th5[29]->Fill(thatSub2->Eta());
				th5[16]->Fill(dr1);
				th5[17]->Fill(dr2);
				th5[34]->Fill((*thisSub1+*thisSub2).Pt());
				th5[35]->Fill((*thatSub1+*thatSub2).Pt());
				th5[40]->Fill((*thisSub1+*thisSub2).Eta());
				th5[41]->Fill((*thatSub1+*thatSub2).Eta());
				th5[44]->Fill((*thisJet+*thatJet).M());
				th5[49]->Fill(fatjetPRmassL2L3Corr[0]);
				th5[50]->Fill(fatjetPRmassL2L3Corr[1]);
				
				th6[8]->Fill(thisSub1->Pt(),scaleFactor);
				th6[9]->Fill(thisSub2->Pt(),scaleFactor);
				th6[10]->Fill(thatSub1->Pt(),scaleFactor);
				th6[11]->Fill(thatSub2->Pt(),scaleFactor);
				th6[26]->Fill(thisSub1->Eta(),scaleFactor);
				th6[27]->Fill(thisSub2->Eta(),scaleFactor);
				th6[28]->Fill(thatSub1->Eta(),scaleFactor);
				th6[29]->Fill(thatSub2->Eta(),scaleFactor);
				th6[16]->Fill(dr1,scaleFactor);
				th6[17]->Fill(dr2,scaleFactor);
				th6[34]->Fill((*thisSub1+*thisSub2).Pt(),scaleFactor);
				th6[35]->Fill((*thatSub1+*thatSub2).Pt(),scaleFactor);
				th6[40]->Fill((*thisSub1+*thisSub2).Eta(),scaleFactor);
				th6[41]->Fill((*thatSub1+*thatSub2).Eta(),scaleFactor);
				th6[44]->Fill((*thisJet+*thatJet).M(),scaleFactor);
				th6[49]->Fill(fatjetPRmassL2L3Corr[0],scaleFactor);
				th6[50]->Fill(fatjetPRmassL2L3Corr[1],scaleFactor);
			}
			if(nbtag==3)nPassB[3]+=scaleFactor;
			if(nbtag==4)nPassB[4]+=scaleFactor;
			
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
	
	TFile* outFile = new TFile(Form("root_files_btaggedEff/%s_2.root",st2.data()),"recreate");
	th2o->Write();
	th2s->Scale(2245.87*Xsec/total);
	th2s->Write();
	for(int i=0;i<6;i++)th3[i]->Write();
	for(int i=0;i<14;i++)th4[i]->Write();
	if(nameRoot==0){
		for(int i=0;i<51;i++){
			//th5[i]->Sumw2();
			th5[i]->Scale(2245.87*Xsec/total);
			//th6[i]->Sumw2();
			th6[i]->Scale(2245.87*Xsec/total);
		}
	}
	for(int i=0;i<51;i++)th5[i]->Write();
	for(int i=0;i<51;i++)th6[i]->Write();
	outFile->Close();
	
}