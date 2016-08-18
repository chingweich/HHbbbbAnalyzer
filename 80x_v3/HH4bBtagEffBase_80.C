
void HH4bBtagEffBase_80(int wMs,int wM, string st,string st2,string option=""){	
	//1=signal ,0=QCD ,2=data-----------------------------------------------------------
	int nameRoot=1;
	if(st2.find("QCD")!= std::string::npos)nameRoot=0;
	if(st2.find("data")!= std::string::npos)nameRoot=2;
	//---------------------Thea correction
	TFile *f3;
	f3=TFile::Open("puppiCorr.root");
	TF1 * puppisd_corrGEN,* puppisd_corrRECO_cen,* puppisd_corrRECO_for;
	puppisd_corrGEN=(TF1 *) f3->FindObjectAny("puppiJECcorr_gen");
	puppisd_corrRECO_cen=(TF1 *) f3->FindObjectAny("puppiJECcorr_reco_0eta1v3");
	puppisd_corrRECO_for=(TF1 *) f3->FindObjectAny("puppiJECcorr_reco_1v3eta2v5");
	TF1* tf1[3];
	tf1[0]=puppisd_corrGEN;
	tf1[1]=puppisd_corrRECO_cen;
	tf1[2]=puppisd_corrRECO_for;
	//option-----------------------------------------------------------
	int JESOption=0;
	if(option.find("JESUp")!= std::string::npos)JESOption=1;
	if(option.find("JESDown")!= std::string::npos)JESOption=2;
	if(option.find("BtagUp")!= std::string::npos)JESOption=3;
	if(option.find("BtagDown")!= std::string::npos)JESOption=4;
	if(option.find("tau21Up")!= std::string::npos)JESOption=5;
	if(option.find("tau21Down")!= std::string::npos)JESOption=6;
	cout<<"JESOption = "<<JESOption<<endl;
	bool printHighPtSubjet=0;
	bool isFast=1;
	//tuple tree and cutflow variables------------------------------------------------------------------------------------
	TFile *f;
	TTree *tree;
	int nPass[20]={0},total=0,dataPassingcsc=0;
	double nPassB[6]={0};
	double fixScaleNum[2]={0};
	//using for Btag Eff -----------------------------------------------------------------------------
	string btagsystematicsType="central";
	if(JESOption==3)btagsystematicsType="up";
	else if(JESOption==4)btagsystematicsType="down";
	BTagCalibration calib("CSVv2L", "subjet_CSVv2_ichep.csv");
	BTagCalibrationReader LF(BTagEntry::OP_LOOSE,"central", {"up", "down"});  
	BTagCalibrationReader HFC(BTagEntry::OP_LOOSE, "central", {"up", "down"});      // other sys types
	BTagCalibrationReader HF(BTagEntry::OP_LOOSE,"central",{"up", "down"});      // other sys types
	LF.load(calib, BTagEntry::FLAV_UDSG,    // btag flavour
            "incl");               // measurement type
	HFC.load(calib, BTagEntry::FLAV_C,    // btag flavour
            "lt");               // measurement type
	HF.load(calib, BTagEntry::FLAV_B,    // btag flavour
            "lt");               // measurement type
	
	TFile *f1;
	if(nameRoot==2)f1=TFile::Open("btagEffSource/data.root");
	else if (nameRoot!=2 && (JESOption==0||JESOption==3||JESOption==4||JESOption==5||JESOption==6))f1=TFile::Open(Form("btagEffSource/%s.root",st2.data()));
	else if (nameRoot!=2 && JESOption==1)f1=TFile::Open(Form("btagEffSource/%s_JESUp.root",st2.data()));
	else if (nameRoot!=2 && JESOption==2)f1=TFile::Open(Form("btagEffSource/%s_JESDown.root",st2.data()));
	TH1D* th1[6];
	string btaggingEff[6]={"effD_b","effN_b","effD_c","effN_c","effD_l","effN_l"};
	for(int i=0;i<6;i++){
		th1[i]=(TH1D*)f1->FindObjectAny(Form("%s_1d",btaggingEff[i].data()));
		if(i==1||i==3||i==5)th1[i]->Divide(th1[i-1]);
	}
	//check for zero btagging SF----------------------------------------------------------------------------------------
	TH2D* th3[6];
	th3[0]=new TH2D("zeroSF_b","zeroSF_b",200,0,2000,60,-3,3);
	th3[1]=new TH2D("zeroSF_c","zeroSF_c",200,0,2000,60,-3,3);
	th3[2]=new TH2D("zeroSF_l","zeroSF_l",200,0,2000,60,-3,3);
	th3[3]=new TH2D("SF_vs_Pt_b","SF_vs_Pt_b",120,0,1200,40,0.8,1.2);
	th3[4]=new TH2D("SF_vs_Pt_c","SF_vs_Pt_c",120,0,1200,40,0.8,1.2);
	th3[5]=new TH2D("SF_vs_Pt_l","SF_vs_Pt_l",120,0,1200,40,0.8,1.2);
	for(int i=0;i<6;i++)th3[i]->Sumw2();
	//check for high btagging SF----------------------------------------------------------------------------------------
	TH1D* th4[14];
	string SF_jet_sub[8]={"SF_jet0_sub0_pass","SF_jet0_sub1_pass","SF_jet1_sub0_pass","SF_jet1_sub1_pass","SF_jet0_sub0_fail","SF_jet0_sub1_fail","SF_jet1_sub0_fail","SF_jet1_sub1_fail"};
	for(int i=0;i<8;i++)th4[i]=new TH1D(Form("%s",SF_jet_sub[i].data()),Form("%s",SF_jet_sub[i].data()),40,0.8,1.2);
	string weightName[6]={"weight","weight_0b","weight_1b","weight_2b","weight_2b_ll","weight_2b_onel"};
	for(int i=0;i<6;i++)th4[i+8]=new TH1D(Form("%s",weightName[i].data()),Form("%s",weightName[i].data()),40,0,2);
	for(int i=0;i<14;i++)th4[i]->Sumw2();
	//saving variables----------------------------------------------------------------------------------------
	TH1D * th5[300],* th6[300];
	for(int i=0;i<2;i++){
		for(int j=0;j<2;j++){
			for(int k=0;k<5;k++){
				th5[(i*2+j)*5+k]=new TH1D(Form("Pt_j%d_sj%d_%db",i,j,k),Form("Pt_j%d_sj%d_%db",i,j,k),200,0,2000);
				th5[(i*2+j)*5+k+20]=new TH1D(Form("Eta_j%d_sj%d_%db",i,j,k),Form("Eta_j%d_sj%d_%db",i,j,k),60,-3,3);
				th5[(i*2+j)*5+k+85]=new TH1D(Form("subCSV_j%d_sj%d_%db",i,j,k),Form("subCSV_j%d_sj%d_%db",i,j,k),20,0,1);
				th5[(i*2+j)*5+k+194]=new TH1D(Form("subCSVCut_j%d_sj%d_%db",i,j,k),Form("subCSVCut_j%d_sj%d_%db",i,j,k),20,0,1);
			}
		}
		for(int k=0;k<5;k++){
			th5[i*5+k+40]=new TH1D(Form("deltaR_j%d_%db",i,k),Form("deltaR_j%d_%db",i,k),20,0,1);
			th5[i*5+k+50]=new TH1D(Form("Pt_j%d_%db",i,k),Form("Pt_j%d_%db",i,k),200,0,2000);
			th5[i*5+k+60]=new TH1D(Form("Eta_j%d_%db",i,k),Form("Eta_j%d_%db",i,k),60,-3,3);
			th5[i*5+k+70]=new TH1D(Form("prMassL2L3_j%d_%db",i,k),Form("prMassL2L3_j%d_%db",i,k),15,90,150);
			th5[i*5+k+105]=new TH1D(Form("tau21_j%d_%db",i,k),Form("tau21_j%d_%db",i,k),25,0,1);
			th5[i*5+k+120]=new TH1D(Form("PuppiSDMassL2L3_j%d_%db",i,k),Form("PuppiSDMassL2L3_j%d_%db",i,k),15,90,150);
			th5[i*5+k+130]=new TH1D(Form("puppiTau21_j%d_%db",i,k),Form("puppiTau21_j%d_%db",i,k),25,0,1);
			th5[i*5+k+140]=new TH1D(Form("prMass_j%d_%db",i,k),Form("prMass_j%d_%db",i,k),15,90,150);
			th5[i*5+k+150]=new TH1D(Form("PuppiSDMass_j%d_%db",i,k),Form("PuppiSDMass_j%d_%db",i,k),15,90,150);
			th5[i*5+k+170]=new TH1D(Form("doubleSV_j%d_%db",i,k),Form("doubleSV_j%d_%db",i,k),40,-1,1);
			th5[i*5+k+184]=new TH1D(Form("FatSV_j%d_%db",i,k),Form("FatSV_j%d_%db",i,k),20,0,1);
			th5[i*5+k+250]=new TH1D(Form("PuppiSDMassThea_j%d_%db",i,k),Form("PuppiSDMassThea_j%d_%db",i,k),15,90,150);
		}
	}
	for(int k=0;k<5;k++){
		th5[k+80]=new TH1D(Form("totalMass_%db",k),Form("totalMass_%db",k),200,1000,5000);
		th5[k+115]=new TH1D(Form("deltaEta_%db",k),Form("deltaEta_%db",k),40,0,2);
		th5[k+160]=new TH1D(Form("logPt_%db",k),Form("logPt_%db",k),70,0,7);
		th5[k+165]=new TH1D(Form("totalMassRed_%db",k),Form("totalMassRed_%db",k),200,1000,5000);
	}
	th5[180]= new TH1D("h_nvtx","h_nvtx",60,0,60);
	th5[181]= new TH1D("h_ntrue","h_ntrue",60,0,60);
	th5[182]= new TH1D("h_nvtx_cut","h_nvtx_cut",60,0,60);
	th5[183]= new TH1D("h_ntrue_cut","h_ntrue_cut",60,0,60);
	
	string prMass_no[4]={"prMass","prMassL2L3","PuppiSDMass","PuppiSDMassL2L3"};
	string tau21_no[2]={"tau21","puppiTau21"};
	for (int i=0;i<4;i++){
		th5[214+i*2]=new TH1D(Form("%s_j0_noPr_noTau21",prMass_no[i].data()),Form("%s_j0_noPr_noTau21",prMass_no[i].data()),200,0,200);
		th5[215+i*2]=new TH1D(Form("%s_j1_noPr_noTau21",prMass_no[i].data()),Form("%s_j1_noPr_noTau21",prMass_no[i].data()),200,0,200);
		th5[222+i*2]=new TH1D(Form("%s_j0_noPr",prMass_no[i].data()),Form("%s_j0_noPr",prMass_no[i].data()),200,0,200);
		th5[223+i*2]=new TH1D(Form("%s_j1_noPr",prMass_no[i].data()),Form("%s_j1_noPr",prMass_no[i].data()),200,0,200);
		th5[230+i*2]=new TH1D(Form("%s_j0_noTau21",prMass_no[i].data()),Form("%s_j0_noTau21",prMass_no[i].data()),15,90,150);
		th5[231+i*2]=new TH1D(Form("%s_j1_noTau21",prMass_no[i].data()),Form("%s_j1_noTau21",prMass_no[i].data()),15,90,150);
	}
	for (int i=0;i<2;i++){
		th5[238+i*2]=new TH1D(Form("%s_j0_noPr_noTau21",tau21_no[i].data()),Form("%s_j0_noPr_noTau21",tau21_no[i].data()),25,0,1);
		th5[239+i*2]=new TH1D(Form("%s_j1_noPr_noTau21",tau21_no[i].data()),Form("%s_j1_noPr_noTau21",tau21_no[i].data()),25,0,1);
		th5[242+i*2]=new TH1D(Form("%s_j0_noPr",tau21_no[i].data()),Form("%s_j0_noPr",tau21_no[i].data()),25,0,1);
		th5[243+i*2]=new TH1D(Form("%s_j1_noPr",tau21_no[i].data()),Form("%s_j1_noPr",tau21_no[i].data()),25,0,1);
		th5[246+i*2]=new TH1D(Form("%s_j0_noTau21",tau21_no[i].data()),Form("%s_j0_noTau21",tau21_no[i].data()),25,0,1);
		th5[247+i*2]=new TH1D(Form("%s_j1_noTau21",tau21_no[i].data()),Form("%s_j1_noTau21",tau21_no[i].data()),25,0,1);
	}
	string flavor[4]={"bb","b","cc","udscg"};
	for (int i=0;i<4;i++){
		th5[260+i]=new TH1D(Form("Pt_j0_0b_%s",flavor[i].data()),Form("Pt_j0_0b_%s",flavor[i].data()),200,0,2000);
		th5[264+i]=new TH1D(Form("Pt_j0_1b_%s",flavor[i].data()),Form("Pt_j0_1b_%s",flavor[i].data()),200,0,2000);
		th5[268+i]=new TH1D(Form("Pt_j0_2b_%s",flavor[i].data()),Form("Pt_j0_2b_%s",flavor[i].data()),200,0,2000);
		th5[272+i]=new TH1D(Form("Pt_j0_DSV_%s",flavor[i].data()),Form("Pt_j0_DSV_%s",flavor[i].data()),200,0,2000);
	}
	
	for(int i=0;i<276;i++){
		th6[i]=(TH1D* )th5[i]->Clone(Form("%ss",th5[i]->GetTitle()));
		th5[i]->Sumw2();
		th6[i]->Sumw2();
	}
	//pileup uncertainty----------------------------------------------------------------------------------------
	TH1D* th7[14];
	th7[0]=new TH1D("totalMass","totalMass",200,1000,5000);
	th7[1]=new TH1D("totalMass_pileup_up","totalMass_pileup_up",200,1000,5000);
	th7[2]=new TH1D("totalMass_pileup_down","totalMass_pileup_down",200,1000,5000);
	th7[3]=new TH1D("pileupEff","pileupEff",15,0.5,15.5);
	double totalPileup[3]={0},passPileup[3]={0};
	standalone_LumiReWeighting LumiWeights_central(0),LumiWeights_up(1),LumiWeights_down(-1);
	//PDF uncertainty
	th7[4]=new TH1D("PDFEff","PDFEff",101,0.5,101.5);
	double passPDF[101]={0},totalPDF[101]={0};
	//QCD uncertainty
	for(int i=5;i<14;i++)th7[i]=new TH1D(Form("uns_QCD_%d",i-5),Form("uns_QCD_%d",i-5),200,1000,5000);
	
	
	
	//NCUtuple loop----------------------------------------------------------------------------------------
	for (int w=wMs;w<wM;w++){
		if(w%20==0)cout<<w<<endl;
		//Get ntuple----------------------------------------------------------------------------------------
		if (nameRoot!=1)f = TFile::Open(Form("%s%d.root",st.data(),w));
		else f = TFile::Open(st.data());
		if (!f || !f->IsOpen())continue;
		TDirectory * dir;
		if (nameRoot!=1)dir = (TDirectory*)f->Get(Form("%s%d.root:/tree",st.data(),w));
		else dir = (TDirectory*)f->Get(Form("%s:/tree",st.data()));
		dir->GetObject("treeMaker",tree);
		TreeReader data(tree);
		total+=data.GetEntriesFast();
		
		TFile* fileFast=TFile::Open(Form("fast/%s/%d.root",st2.data(),w));
		TTree* treeFast=(TTree* )fileFast->Get("hh4bFast");
		TreeReader dataFast(treeFast);
		//if(isFast)dataFast=TreeReader(treeFast);
		
		for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){//event loop----------------------------------------------------------------------------------------
			data.GetEntry(jEntry);
			if(isFast)dataFast.GetEntry(jEntry);
			
			Int_t nVtx        = data.GetInt("nVtx");
			Float_t ntrue= data.GetFloat("pu_nTrueInt");
			//Float_t*  pdfscaleSysWeights= data.GetPtrFloat("pdfscaleSysWeights");
			th5[180]->Fill(nVtx);
			th5[181]->Fill(ntrue);
			double PU_weight[3]={1,1,1};
			
			if(nameRoot!=2){
				if(ntrue<51){
					PU_weight[0] = LumiWeights_central.weight(ntrue);
					PU_weight[1]= LumiWeights_up.weight(ntrue);
					PU_weight[2] = LumiWeights_down.weight(ntrue);
				}
				else {
					PU_weight[0] = LumiWeights_central.weight(50);
					PU_weight[1] = LumiWeights_up.weight(50);
					PU_weight[2]= LumiWeights_down.weight(50);
				}
			}
			
			th6[180]->Fill(nVtx,PU_weight[0]);
			th6[181]->Fill(ntrue,PU_weight[0]);
			fixScaleNum[0]+=PU_weight[0];
			for(int i=0;i<3;i++)totalPileup[i]+=PU_weight[i];
			for(int i=0;i<101;i++)totalPDF[i]+=PU_weight[0];
			
			Int_t nPassFast;  
			if(isFast)nPassFast= dataFast.GetInt("nPassB");
			//cout<<nPassFast<<endl;
			if(isFast &&nPassFast<8)continue;
			//0. has a good vertex
			if(nVtx<1)continue;nPass[0]++;
			//1.trigger
			std::string* trigName = data.GetPtrString("hlt_trigName");
		 	vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
			bool passTrigger=false;
			for(int it=0; it< data.GetPtrStringSize(); it++){
				std::string thisTrig= trigName[it];
				bool results = trigResult[it];
				if( ((thisTrig.find("HLT_PFHT800")!= std::string::npos||
						//thisTrig.find("HLT_PFHT650")!= std::string::npos||
						thisTrig.find("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v")!= std::string::npos||
						thisTrig.find("HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v")!= std::string::npos||
						thisTrig.find("HLT_AK8PFJet360_TrimMass30_v")!= std::string::npos||
						thisTrig.find("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v")!= std::string::npos
						) && results==1)){
					passTrigger=true;
					break;
				}
			}
			if(!passTrigger && nameRoot==2)continue;nPass[1]++;
			
			int nFATJet         = data.GetInt("FATnJet");
			TClonesArray* fatjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
			float*  FATjetCorrUncUp = data.GetPtrFloat("FATjetCorrUncUp"); 
			float*  FATjetCorrUncDown = data.GetPtrFloat("FATjetCorrUncDown"); 
			vector<bool>    &FATjetPassIDTight = *((vector<bool>*) data.GetPtr("FATjetPassIDTight"));
			//2.nJets
			if(nFATJet<2)continue;nPass[2]++;
			TLorentzVector* thisJet ,* thatJet;
			if (JESOption==1){
				TLorentzVector test0= (*((TLorentzVector*)fatjetP4->At(0)))*(1+FATjetCorrUncUp[0] );
				TLorentzVector test1= (*((TLorentzVector*)fatjetP4->At(1)))*(1+FATjetCorrUncUp[1] );
				thisJet= &test0;
				thatJet= &test1;
			}
			else if (JESOption==2){
				TLorentzVector test0= (*((TLorentzVector*)fatjetP4->At(0)))*(1-FATjetCorrUncDown[0] );
				TLorentzVector test1= (*((TLorentzVector*)fatjetP4->At(1)))*(1-FATjetCorrUncDown[1] );
				thisJet= &test0;
				thatJet= &test1;
			}
			else{
				thisJet= (TLorentzVector*)fatjetP4->At(0);
				thatJet = (TLorentzVector*)fatjetP4->At(1);
			}
			//3. Pt 
			if(thisJet->Pt()<200||thatJet->Pt()<200)continue;
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
			if(dEta>1.3)continue;
			nPass[6]++;
			//7. Mjj-----------------------------------------
			float mjj = (*thisJet+*thatJet).M();
			float mjjRed = (*thisJet+*thatJet).M()+250-thisJet->M()-thatJet->M();
			if(mjjRed<1000)continue;
			nPass[7]++;
			//8. fatjetPRmassL2L3Corr-----------------------------------------
			Float_t*  fatjetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
			Float_t*  FATjetPuppiSDmassL2L3Corr = data.GetPtrFloat("FATjetPuppiSDmassL2L3Corr");
			Float_t*  FATjetPRmass = data.GetPtrFloat("FATjetPRmass");
			Float_t*  FATjetPuppiSDmass = data.GetPtrFloat("FATjetPuppiSDmass");
			
			
			
			vector<float>   *subjetSDCSV =  data.GetPtrVectorFloat("FATsubjetSDCSV");
			vector<float>   *subjetSDPx  =  data.GetPtrVectorFloat("FATsubjetSDPx", nFATJet);
			vector<float>   *subjetSDPy  =  data.GetPtrVectorFloat("FATsubjetSDPy", nFATJet);
			vector<float>   *subjetSDPz  =  data.GetPtrVectorFloat("FATsubjetSDPz", nFATJet);
			vector<float>   *subjetSDE   =  data.GetPtrVectorFloat("FATsubjetSDE", nFATJet);
			vector<Int_t>   *FATsubjetSDHadronFlavor =  data.GetPtrVectorInt("FATsubjetSDHadronFlavor");
			int nbtag=0,nbtag2=0;
			float MaxBJetPt = 670., MaxLJetPt = 1000.;
			double sf[2][2],eta[2],pt[2],dr[2],subjetPt[2][2],subjetEta[2][2],eff[2][2],btaggingscaleFactor=1;
			pt[0]=thisJet->Pt();
			pt[1]=thatJet->Pt();
			TLorentzVector* subjetP4[2][2];
			for(int i=0;i<2;i++){
				for(int j=0;j<2;j++){
					sf[i][j]=1;
					subjetP4[i][j]=new TLorentzVector(0,0,0,0);
					subjetP4[i][j]->SetPxPyPzE(subjetSDPx[i][j],subjetSDPy[i][j],subjetSDPz[i][j],subjetSDE[i][j]);
					subjetPt[i][j]=subjetP4[i][j]->Pt();
				}
			}
			for(int i=0;i<2;i++){
				for(int j=0;j<2;j++){
					
					if(subjetPt[i][j]>2000 && printHighPtSubjet){
						Float_t  pfMetCorrPt = data.GetFloat("pfMetCorrPt");
						Long64_t  runId=data.GetLong64("runId");
						Long64_t  lumiSection=data.GetLong64("lumiSection");
						Long64_t  eventId=data.GetLong64("eventId");
						TLorentzVector* tempTL;
						if(i==0)tempTL=thisJet;
						else tempTL=thatJet;
						cout<<"run="<<runId<<",lumi="<<lumiSection<<",event="<<eventId<<",met="<<pfMetCorrPt<<endl<<
						"     ,FatPt="<<pt[i]<<",fatEta="<<tempTL->Eta()<<",fatPhi="<<tempTL->Phi()<<endl<<
						"     ,leadingSubPt="<<subjetPt[i][0]<<",Eta="<<subjetP4[i][0]->Eta()<<",phi="<<subjetP4[i][0]->Phi()<<",flavor="<<FATsubjetSDHadronFlavor[i][0]<<",csv="<<subjetSDCSV[i][0]<<endl<<
						"     ,subleadingSubPt="<<subjetPt[i][1]<<",Eta="<<subjetP4[i][1]->Eta()<<",phi="<<subjetP4[i][1]->Phi()<<",flavor="<<FATsubjetSDHadronFlavor[i][1]<<",csv="<<subjetSDCSV[i][1]<<endl;
					}
						
					
					
					subjetEta[i][j]=subjetP4[i][j]->Eta();
					//get btagging eff------------------------------------------------------------
					int getPtBin=1;
					if(subjetPt[i][j]<140)getPtBin=1;
					else if (140<=subjetPt[i][j] && subjetPt[i][j]<180)getPtBin=2;
					else if (180<=subjetPt[i][j] && subjetPt[i][j]<240)getPtBin=3;
					else getPtBin=4;
					//if(FATsubjetSDHadronFlavor[i][j]==5)eff[i][j]=th1[1]->GetBinContent(getPtBin,ceil(subjetEta[i][j]/0.2)+15);
					//else if(FATsubjetSDHadronFlavor[i][j]==4)eff[i][j]=th1[3]->GetBinContent(getPtBin,ceil(subjetEta[i][j]/0.2)+15);
					//else eff[i][j]=th1[5]->GetBinContent(getPtBin,ceil(subjetEta[i][j]/0.2)+15);
					if(FATsubjetSDHadronFlavor[i][j]==5)eff[i][j]=th1[1]->GetBinContent(getPtBin);
					else if(FATsubjetSDHadronFlavor[i][j]==4)eff[i][j]=th1[3]->GetBinContent(getPtBin);
					else {
						int temp=0;
						if(subjetPt[i][j]>=900)temp=10;
						else temp=ceil(subjetPt[i][j]/100);
						
						bool checkBinContentIfZero=0;
						while(checkBinContentIfZero==0){
							if(th1[4]->GetBinContent(temp)==0){
								temp--;
							}
							else checkBinContentIfZero=1;
						}
						eff[i][j]=th1[5]->GetBinContent(temp);
					}
					//if(eff[i][j]==0)cout<<FATsubjetSDHadronFlavor[i][j]<<",="<<subjetPt[i][j]<<",="<<eff[i][j]<<"!!"<<ceil(subjetPt[i][j]/100)<<endl;
					//check maxPt-------------------------------------------------------------
					//if(FATsubjetSDHadronFlavor[i][j]!=0 && subjetPt[i][j]>MaxBJetPt )subjetPt[i][j]=MaxBJetPt-0.1;
					//if(FATsubjetSDHadronFlavor[i][j]==0 && subjetPt[i][j]>MaxLJetPt )subjetPt[i][j]=MaxLJetPt-0.1;
					//Get SF from csv------------------------------------------------------------
					if(FATsubjetSDHadronFlavor[i][j]==5){
						sf[i][j]=HF.eval_auto_bounds("central",BTagEntry::FLAV_B, subjetEta[i][j],subjetPt[i][j]); 
						//sf[i][j]=HF.eval(BTagEntry::FLAV_B,subjetEta[i][j],subjetPt[i][j]); 
						th3[3]->Fill(subjetPt[i][j],sf[i][j]);
					}
					else if(FATsubjetSDHadronFlavor[i][j]==4){
						sf[i][j]=HFC.eval_auto_bounds("central",BTagEntry::FLAV_C, subjetEta[i][j],subjetPt[i][j]); 
						//sf[i][j]=HF.eval(BTagEntry::FLAV_C,subjetEta[i][j],subjetPt[i][j]); 
						th3[4]->Fill(subjetPt[i][j],sf[i][j]);
					}
					else {
						sf[i][j]=LF.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, subjetEta[i][j],subjetPt[i][j]); 
						//sf[i][j]=LF.eval(BTagEntry::FLAV_UDSG,subjetEta[i][j],subjetPt[i][j]); 
						th3[5]->Fill(subjetPt[i][j],sf[i][j]);
					}
					//check zero ------------------------------------------------------------
					if(sf[i][j]==0 && FATsubjetSDHadronFlavor[i][j]==5 ) th3[0]->Fill(subjetPt[i][j],subjetEta[i][j]);
					if(sf[i][j]==0 && FATsubjetSDHadronFlavor[i][j]==4 ) th3[1]->Fill(subjetPt[i][j],subjetEta[i][j]);
					if(sf[i][j]==0 && FATsubjetSDHadronFlavor[i][j]!=5 && FATsubjetSDHadronFlavor[i][j]!=4 ) th3[2]->Fill(subjetPt[i][j],subjetEta[i][j]);
					//conut nbtag---------------------------------------------------------
					if(subjetSDCSV[i][j]>0.46)nbtag++;
					//get tot. btagging SF
					if(subjetSDCSV[i][j]>=0.46)btaggingscaleFactor*=sf[i][j];
					else btaggingscaleFactor*=((1-eff[i][j]*sf[i][j])/(1-eff[i][j]));
					//##############check light jet SF##########
					if(subjetSDCSV[i][j]>0.46 &&FATsubjetSDHadronFlavor[i][j]==0)nbtag2++;
					if(subjetSDCSV[i][j]>0.46)th4[i*2+j]->Fill(sf[i][j]);
					else th4[i*2+j+4]->Fill(sf[i][j]);
					//subjetPt[i][j]=subjetP4[i][j]->Pt();
				}
				dr[i]=subjetP4[i][0]->DeltaR(*subjetP4[i][1]);
			}
			
			if(btaggingscaleFactor>3){
				for(int i=0;i<2;i++){
					for(int j=0;j<2;j++){
						//cout<<eff[i][j]<<",="<<sf[i][j]<<endl;
					}
				}
				//cout<<jEntry<<endl;
			}
			
			
			double scaleFactor=PU_weight[0]*btaggingscaleFactor;
			Float_t*  fatjetTau1 = data.GetPtrFloat("FATjetTau1");
			Float_t*  fatjetTau2 = data.GetPtrFloat("FATjetTau2");
			Float_t*  FATjetPuppiTau1 = data.GetPtrFloat("FATjetPuppiTau1");
			Float_t*  FATjetPuppiTau2 = data.GetPtrFloat("FATjetPuppiTau2");
			double tau21[2];
			tau21[0]=(fatjetTau2[0]/fatjetTau1[0]),tau21[1]=(fatjetTau2[1]/fatjetTau1[1]);
			double puppiTau21[2];
			puppiTau21[0]=(FATjetPuppiTau2[0]/FATjetPuppiTau1[0]),puppiTau21[1]=(FATjetPuppiTau2[1]/FATjetPuppiTau1[1]);
			
			for (int i=0;i<2;i++){
				th5[214+i]->Fill(FATjetPRmass[i]);
				th5[216+i]->Fill(fatjetPRmassL2L3Corr[i]);
				th5[218+i]->Fill(FATjetPuppiSDmass[i]);
				th5[220+i]->Fill(FATjetPuppiSDmassL2L3Corr[i]);
				th5[238+i]->Fill(tau21[i]);
				th5[240+i]->Fill(puppiTau21[i]);
				
				th6[214+i]->Fill(FATjetPRmass[i],scaleFactor);
				th6[216+i]->Fill(fatjetPRmassL2L3Corr[i],scaleFactor);
				th6[218+i]->Fill(FATjetPuppiSDmass[i],scaleFactor);
				th6[220+i]->Fill(FATjetPuppiSDmassL2L3Corr[i],scaleFactor);
				th6[238+i]->Fill(tau21[i],scaleFactor);
				th6[240+i]->Fill(puppiTau21[i],scaleFactor);
			}
			
			if(!(tau21[0]>0.6 || tau21[1]>0.6)){
				for (int i=0;i<2;i++){
					th5[222+i]->Fill(FATjetPRmass[i]);
					th5[224+i]->Fill(fatjetPRmassL2L3Corr[i]);
					th5[226+i]->Fill(FATjetPuppiSDmass[i]);
					th5[228+i]->Fill(FATjetPuppiSDmassL2L3Corr[i]);
					th5[242+i]->Fill(tau21[i]);
					th5[244+i]->Fill(puppiTau21[i]);
					
					th6[222+i]->Fill(FATjetPRmass[i],scaleFactor);
					th6[224+i]->Fill(fatjetPRmassL2L3Corr[i],scaleFactor);
					th6[226+i]->Fill(FATjetPuppiSDmass[i],scaleFactor);
					th6[228+i]->Fill(FATjetPuppiSDmassL2L3Corr[i],scaleFactor);
					th6[242+i]->Fill(tau21[i],scaleFactor);
					th6[244+i]->Fill(puppiTau21[i],scaleFactor);
				}
			}
			
			
			if(fatjetPRmassL2L3Corr[0]<105||fatjetPRmassL2L3Corr[0]>135)continue;
			if(fatjetPRmassL2L3Corr[1]<105||fatjetPRmassL2L3Corr[1]>135)continue;
			nPass[8]++;
			for (int i=0;i<2;i++){
				th5[230+i]->Fill(FATjetPRmass[i]);
				th5[232+i]->Fill(fatjetPRmassL2L3Corr[i]);
				th5[234+i]->Fill(FATjetPuppiSDmass[i]);
				th5[236+i]->Fill(FATjetPuppiSDmassL2L3Corr[i]);
				th5[246+i]->Fill(tau21[i]);
				th5[248+i]->Fill(puppiTau21[i]);
				
				th6[230+i]->Fill(FATjetPRmass[i],scaleFactor);
				th6[232+i]->Fill(fatjetPRmassL2L3Corr[i],scaleFactor);
				th6[234+i]->Fill(FATjetPuppiSDmass[i],scaleFactor);
				th6[236+i]->Fill(FATjetPuppiSDmassL2L3Corr[i],scaleFactor);
				th6[246+i]->Fill(tau21[i],scaleFactor);
				th6[248+i]->Fill(puppiTau21[i],scaleFactor);
			}
			
			//9.-----------------------------------------
			if(tau21[0]>0.6 || tau21[1]>0.6) continue;
			nPass[9]++;
			double tau21_SF=1.031*0.881;
			if(JESOption==5)tau21_SF=(1.031+0.126)*(0.881+0.49);
			if(JESOption==6)tau21_SF=(1.031-0.126)*(0.881-0.49);
			if(tau21[0]<0.6 && tau21[1]<0.6 ){
				tau21_SF=1.031*1.031;
				if(JESOption==5)tau21_SF=(1.031+0.126)*(1.031+0.126);
				if(JESOption==6)tau21_SF=(1.031-0.126)*(1.031-0.126);
			}
			
			Int_t* FATjetHadronFlavor        = data.GetPtrInt("FATjetHadronFlavor");
			
			if(subjetSDCSV[0][0]<0.46 && subjetSDCSV[0][1]<0.46){
				if(FATjetHadronFlavor[0]==5 && FATsubjetSDHadronFlavor[0][0]==5 && FATsubjetSDHadronFlavor[0][1]==5)th5[260]->Fill(thisJet->Pt());
				else if(FATjetHadronFlavor[0]==5 && (FATsubjetSDHadronFlavor[0][0]==5 || FATsubjetSDHadronFlavor[0][1]==5))th5[261]->Fill(thisJet->Pt());
				else if(FATjetHadronFlavor[0]==4 && FATsubjetSDHadronFlavor[0][0]==4 && FATsubjetSDHadronFlavor[0][1]==4)th5[262]->Fill(thisJet->Pt());
				else th5[263]->Fill(thisJet->Pt());
			}
			else if (subjetSDCSV[0][0]>0.46 && subjetSDCSV[0][1]>0.46){
				if(FATjetHadronFlavor[0]==5 && FATsubjetSDHadronFlavor[0][0]==5 && FATsubjetSDHadronFlavor[0][1]==5)th5[268]->Fill(thisJet->Pt());
				else if(FATjetHadronFlavor[0]==5 && (FATsubjetSDHadronFlavor[0][0]==5 || FATsubjetSDHadronFlavor[0][1]==5))th5[269]->Fill(thisJet->Pt());
				else if (FATjetHadronFlavor[0]==4 && FATsubjetSDHadronFlavor[0][0]==4 && FATsubjetSDHadronFlavor[0][1]==4)th5[270]->Fill(thisJet->Pt());
				else th5[271]->Fill(thisJet->Pt());
			}
			else {
				if(FATjetHadronFlavor[0]==5 && FATsubjetSDHadronFlavor[0][0]==5 && FATsubjetSDHadronFlavor[0][1]==5)th5[264]->Fill(thisJet->Pt());
				else if(FATjetHadronFlavor[0]==5 && (FATsubjetSDHadronFlavor[0][0]==5 || FATsubjetSDHadronFlavor[0][1]==5))th5[265]->Fill(thisJet->Pt());
				else if (FATjetHadronFlavor[0]==4 && FATsubjetSDHadronFlavor[0][0]==4 && FATsubjetSDHadronFlavor[0][1]==4)th5[266]->Fill(thisJet->Pt());
				else th5[267]->Fill(thisJet->Pt());
			}
			
			
			
			
			
			
			//10.btag
			
			th4[8]->Fill(btaggingscaleFactor);
			if(nbtag2==2 && nbtag==2)th4[12]->Fill(btaggingscaleFactor);
			if(nbtag2==1 && nbtag==2)th4[13]->Fill(btaggingscaleFactor);
			
			//uncertainty -------------------------------------
			//double scaleFactor=btaggingscaleFactor*PU_weight[0]*tau21_SF;
			
			for(int i=0;i<3;i++){
				passPileup[i]+=btaggingscaleFactor*PU_weight[i]*tau21_SF;
				th7[i]->Fill(mjj,btaggingscaleFactor*PU_weight[i]*tau21_SF);
			}
			 for(int i=0;i<101;i++)passPDF[i]+=btaggingscaleFactor*PU_weight[0]*tau21_SF;
			 for(int i=5;i<14;i++)th7[i]->Fill(mjj,btaggingscaleFactor*PU_weight[0]*tau21_SF);
			 fixScaleNum[1]+=PU_weight[0];
			//--------------------------------------
			Float_t  FATjetPuppiSDmassThea[2] ={0};
			FATjetPuppiSDmassThea[0]=FATjetPuppiSDmass[0]*getPUPPIweight(thisJet->Pt(),thisJet->Eta(),tf1);
			FATjetPuppiSDmassThea[1]=FATjetPuppiSDmass[1]*getPUPPIweight(thatJet->Pt(),thatJet->Eta(),tf1);
			eta[0]=thisJet->Eta();
			eta[1]=thatJet->Eta();
			//Float_t*  FATjet_DoubleSV = data.GetPtrFloat("FATjet_DoubleSV");
			Float_t*  ADDjet_DoubleSV = data.GetPtrFloat("ADDjet_DoubleSV");
			Float_t  FATjet_DoubleSV[2]={0};
			
			Int_t ADDnJet        = data.GetInt("ADDnJet");
			bool matchThis=0,matchThat=0;
			TClonesArray* ADDjetP4 = (TClonesArray*) data.GetPtrTObject("ADDjetP4");
			for(int i=0;i<ADDnJet;i++){
				TLorentzVector* thisAddJet ;
				thisAddJet= (TLorentzVector*)ADDjetP4->At(i);
				if(!matchThis && thisAddJet->DeltaR(*thisJet)<0.8){
					matchThis=1;
					FATjet_DoubleSV[0]=ADDjet_DoubleSV[i];
					continue;
				}
				if(!matchThat && thisAddJet->DeltaR(*thatJet)<0.8){
					matchThat=1;
					FATjet_DoubleSV[1]=ADDjet_DoubleSV[i];
				}
				if(matchThis&& matchThat){
					//cout<<"match"<<FATjet_DoubleSV[0]<<","<<FATjet_DoubleSV[0]<<endl;
					break;
				}
				
			}
			
			if(FATjet_DoubleSV[0]>0.6){
				if(FATjetHadronFlavor[0]==5 && FATsubjetSDHadronFlavor[0][0]==5 && FATsubjetSDHadronFlavor[0][1]==5)th5[272]->Fill(thisJet->Pt());
				else if(FATjetHadronFlavor[0]==5 && (FATsubjetSDHadronFlavor[0][0]==5 || FATsubjetSDHadronFlavor[0][1]==5))th5[273]->Fill(thisJet->Pt());
				else if (FATjetHadronFlavor[0]==4 && FATsubjetSDHadronFlavor[0][0]==4 && FATsubjetSDHadronFlavor[0][1]==4)th5[274]->Fill(thisJet->Pt());
				else th5[275]->Fill(thisJet->Pt());
			}
			
			
			
			Float_t*  FATjetCISVV2 = data.GetPtrFloat("FATjetCISVV2");
			for(int i=0;i<2;i++){
				for(int j=0;j<2;j++){
					for(int k=0;k<5;k++){
						if(nbtag!=k)continue;
						th5[(i*2+j)*5+k]->Fill(subjetPt[i][j]);
						th5[(i*2+j)*5+k+20]->Fill(subjetEta[i][j]);
						th5[(i*2+j)*5+k+85]->Fill(subjetSDCSV[i][j]);
						th6[(i*2+j)*5+k]->Fill(subjetPt[i][j],scaleFactor);
						th6[(i*2+j)*5+k+20]->Fill(subjetEta[i][j],scaleFactor);
						th6[(i*2+j)*5+k+85]->Fill(subjetSDCSV[i][j],scaleFactor);
						if(subjetPt[i][j]>30 && fabs(subjetEta[i][j])<2.4)th5[(i*2+j)*5+k+194]->Fill(subjetSDCSV[i][j]);
						if(subjetPt[i][j]>30 && fabs(subjetEta[i][j])<2.4)th6[(i*2+j)*5+k+194]->Fill(subjetSDCSV[i][j],scaleFactor);
					}
				}
				for(int k=0;k<5;k++){
					if(nbtag!=k)continue;
					th5[i*5+k+40]->Fill(dr[i]);
					th5[i*5+k+50]->Fill(pt[i]);
					th5[i*5+k+60]->Fill(eta[i]);
					th5[i*5+k+70]->Fill(fatjetPRmassL2L3Corr[i]);
					th5[i*5+k+105]->Fill(tau21[i]);
					th6[i*5+k+40]->Fill(dr[i],scaleFactor);
					th6[i*5+k+50]->Fill(pt[i],scaleFactor);
					th6[i*5+k+60]->Fill(eta[i],scaleFactor);
					th6[i*5+k+70]->Fill(fatjetPRmassL2L3Corr[i],scaleFactor);
					th6[i*5+k+105]->Fill(tau21[i],scaleFactor);
					th5[i*5+k+120]->Fill(FATjetPuppiSDmassL2L3Corr[i]);
					th6[i*5+k+120]->Fill(FATjetPuppiSDmassL2L3Corr[i],scaleFactor);
					th5[i*5+k+130]->Fill(puppiTau21[i]);
					th6[i*5+k+130]->Fill(puppiTau21[i],scaleFactor);
					th5[i*5+k+140]->Fill(FATjetPRmass[i]);
					th6[i*5+k+140]->Fill(FATjetPRmass[i],scaleFactor);
					th5[i*5+k+150]->Fill(FATjetPuppiSDmass[i]);
					th6[i*5+k+150]->Fill(FATjetPuppiSDmass[i],scaleFactor);
					th5[i*5+k+170]->Fill(FATjet_DoubleSV[i]);
					th6[i*5+k+170]->Fill(FATjet_DoubleSV[i],PU_weight[0]);
					th5[i*5+k+184]->Fill(FATjetCISVV2[i]);
					th6[i*5+k+184]->Fill(FATjetCISVV2[i],scaleFactor);
					th5[i*5+k+250]->Fill(FATjetPuppiSDmassThea[i]);
					th6[i*5+k+250]->Fill(FATjetPuppiSDmassThea[i],scaleFactor);
				}
			}
			for(int k=0;k<5;k++){
				if(nbtag!=k)continue;
				th5[k+80]->Fill(mjj);
				th6[k+80]->Fill(mjj,scaleFactor);
				th5[k+115]->Fill(dEta);
				th6[k+115]->Fill(dEta,scaleFactor);
				th5[k+160]->Fill(log10(pt[0]));
				th6[k+160]->Fill(log10(pt[0]),scaleFactor);
				th5[k+160]->Fill(log10(pt[1]));
				th6[k+160]->Fill(log10(pt[1]),scaleFactor);
				th5[k+165]->Fill(mjjRed);
				th6[k+165]->Fill(mjjRed,scaleFactor);
				nPassB[k]+=scaleFactor;
				nPass[k+10]++;
				if(k<3)th4[9+k]->Fill(btaggingscaleFactor);
			}
			
			th5[182]->Fill(nVtx);
			th5[183]->Fill(ntrue);
			th6[182]->Fill(nVtx,scaleFactor);
			th6[183]->Fill(ntrue,scaleFactor);

		}//end event loop----------------------------------------------------------------------------------------
	}	//end ntuple loop----------------------------------------------------------------------------------------
	cout<<"entries="<<total<<endl;	
	for(int i=0;i<6;i++)cout<<"nPassB["<<i<<"]="<<nPassB[i]<<endl;
	for(int i=0;i<16;i++)cout<<"nPass["<<i<<"]="<<nPass[i]<<endl;
	
	for(int i=0;i<3;i++)th7[3]->SetBinContent(i+1,passPileup[i]/totalPileup[i]);
	for(int i=0;i<101;i++)th7[4]->SetBinContent(i+1,passPDF[i]/totalPDF[i]);
	
	TH1D * th2o=new TH1D("Nbtagjet","Nbtagjet",5,-0.5,4.5);
	for (int i=0;i<5;i++){
		if(nameRoot==2 && i>2)continue;
		th2o->SetBinContent(i+1,nPass[i+10]);
	}
	TH1D * th2ob=new TH1D("NbtagjetB","NbtagjetB",5,-0.5,4.5);
	for (int i=0;i<5;i++){
		if(nameRoot==2 && i>2)continue;
		th2ob->SetBinContent(i+1,nPassB[i]);
	}
	
	TH1D * fixScale=new TH1D("fixScale","fixScale",2,-0.5,1.5);
	fixScale->SetBinContent(1,fixScaleNum[0]);
	fixScale->SetBinContent(2,fixScaleNum[1]);
	TH1D * cutflow=new TH1D("cutflow","cutflow",16,0.5,16.5);
	cutflow->SetBinContent(1,total);
	if(nameRoot==2)for(int ii=1;ii<14;ii++)cutflow->SetBinContent(ii+1,nPass[ii-1]);
	else for(int ii=1;ii<16;ii++)cutflow->SetBinContent(ii+1,nPass[ii-1]);
	
	TFile* outFile ;
	if(JESOption==0)outFile= new TFile(Form("sf/%s.root",st2.data()),"recreate");
	else if(JESOption==1)outFile= new TFile(Form("sf/%s_JESUp.root",st2.data()),"recreate");
	else if(JESOption==2)outFile= new TFile(Form("sf/%s_JESDown.root",st2.data()),"recreate");
	else if(JESOption==3)outFile= new TFile(Form("sf/%s_BtagUp.root",st2.data()),"recreate");
	else if(JESOption==4)outFile= new TFile(Form("sf/%s_BtagDown.root",st2.data()),"recreate");
	else if(JESOption==5)outFile= new TFile(Form("sf/%s_tau21Up.root",st2.data()),"recreate");
	else if(JESOption==6)outFile= new TFile(Form("sf/%s_tau21Down.root",st2.data()),"recreate");
	th2o->Write();
	th2ob->Write();
	fixScale->Write();
	cutflow->Write();
	for(int i=0;i<6;i++)th3[i]->Write();
	for(int i=0;i<14;i++)th4[i]->Write();
	for(int i=0;i<276;i++){
		th5[i]->Write();
		th6[i]->Write();
	}
	for(int i=0;i<14;i++)th7[i]->Write();
	outFile->Close();
}