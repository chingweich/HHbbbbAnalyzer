
float getPUPPIweight(float puppipt, float puppieta ){

   TF1* puppisd_corrGEN      = new TF1("puppisd_corrGEN","[0]+[1]*pow(x*[2],-[3])");
  puppisd_corrGEN->SetParameters(
   				 1.00626,
   				 -1.06161,
   				 0.07999,
   				 1.20454
   				 );
  TF1* puppisd_corrRECO_cen = new TF1("puppisd_corrRECO_cen","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)");
  puppisd_corrRECO_cen->SetParameters(
   				      1.05807,
   				      -5.91971e-05,
   				      2.296e-07,
   				      -1.98795e-10,
   				      6.67382e-14,
   				      -7.80604e-18
   				      );

  TF1* puppisd_corrRECO_for = new TF1("puppisd_corrRECO_for","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)");
  puppisd_corrRECO_for->SetParameters(
   				      1.26638,
   				      -0.000658496,
   				      9.73779e-07,
   				      -5.93843e-10,
   				      1.61619e-13,
   				      -1.6272e-17);

  float genCorr  = 1.;
  float recoCorr = 1.;
  float totalWeight = 1.;

  genCorr =  puppisd_corrGEN->Eval( puppipt );
  if( fabs(puppieta)  <= 1.3 ) recoCorr = puppisd_corrRECO_cen->Eval( puppipt );
  else
    if( fabs(puppieta) > 1.3 ) recoCorr = puppisd_corrRECO_for->Eval( puppipt );

  totalWeight = genCorr * recoCorr;
  // file->Close();
  return totalWeight;
}


void HH4bBtagEffBase_80_v2(int wMs,int wM, string st,string st2,string option=""){	
cout<<Form("sf2/%s.root",st2.data())<<endl;
	//1=signal ,0=QCD ,2=data-----------------------------------------------------------
	int nameRoot=1;
	if((st2.find("QCD")!= std::string::npos)||
	(st2.find("bGen")!= std::string::npos)||
	(st2.find("bEnriched")!= std::string::npos))nameRoot=0;
	if(st2.find("data")!= std::string::npos)nameRoot=2;
	
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
	double nPass[30]={0},total=0;
	double fixScaleNum[2]={0};

	//saving variables----------------------------------------------------------------------------------------
	TH1D * th5[300],* th_flavor[4][300];
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
			th5[i*5+k+70]=new TH1D(Form("prMassL2L3_j%d_%db",i,k),Form("prMassL2L3_j%d_%db",i,k),40,40,200);
			th5[i*5+k+105]=new TH1D(Form("tau21_j%d_%db",i,k),Form("tau21_j%d_%db",i,k),25,0,1);
			th5[i*5+k+120]=new TH1D(Form("PuppiSDMassL2L3_j%d_%db",i,k),Form("PuppiSDMassL2L3_j%d_%db",i,k),40,40,200);
			th5[i*5+k+130]=new TH1D(Form("puppiTau21_j%d_%db",i,k),Form("puppiTau21_j%d_%db",i,k),25,0,1);
			th5[i*5+k+140]=new TH1D(Form("prMass_j%d_%db",i,k),Form("prMass_j%d_%db",i,k),15,90,150);
			th5[i*5+k+150]=new TH1D(Form("PuppiSDMass_j%d_%db",i,k),Form("PuppiSDMass_j%d_%db",i,k),40,40,200);
			th5[i*5+k+170]=new TH1D(Form("doubleSV_j%d_%db",i,k),Form("doubleSV_j%d_%db",i,k),40,-1,1);
			th5[i*5+k+184]=new TH1D(Form("AK8SV_j%d_%db",i,k),Form("AK8SV_j%d_%db",i,k),20,0,1);
			th5[i*5+k+250]=new TH1D(Form("PuppiSDMassThea_j%d_%db",i,k),Form("PuppiSDMassThea_j%d_%db",i,k),40,40,200);
		}
	}
	for(int k=0;k<5;k++){
		th5[k+80]=new TH1D(Form("totalMass_%db",k),Form("totalMass_%db",k),200,1000,5000);
		th5[k+115]=new TH1D(Form("deltaEta_%db",k),Form("deltaEta_%db",k),40,0,2);
		th5[k+160]=new TH1D(Form("logPt_%db",k),Form("logPt_%db",k),70,0,7);
		th5[k+165]=new TH1D(Form("totalMassRed_%db",k),Form("totalMassRed_%db",k),200,1000,5000);
		th5[k+260]=new TH1D(Form("HT_%db",k),Form("HT_%db",k),500,500,3000);
	}
	th5[180]= new TH1D("h_nvtx","h_nvtx",60,0,60);
	th5[181]= new TH1D("h_ntrue","h_ntrue",60,0,60);
	th5[182]= new TH1D("h_nvtx_cut","h_nvtx_cut",60,0,60);
	th5[183]= new TH1D("h_ntrue_cut","h_ntrue_cut",60,0,60);
	
	string prMass_no[5]={"prMass","prMassL2L3","PuppiSDMass","PuppiSDMassL2L3","PuppiSDMassThea"};
	string tau21_no[2]={"tau21","puppiTau21"};
	for (int i=0;i<4;i++){
		th5[214+i*2]=new TH1D(Form("%s_j0_noPr_noTau21",prMass_no[i].data()),Form("%s_j0_noPr_noTau21",prMass_no[i].data()),200,0,200);
		th5[215+i*2]=new TH1D(Form("%s_j1_noPr_noTau21",prMass_no[i].data()),Form("%s_j1_noPr_noTau21",prMass_no[i].data()),200,0,200);
		th5[222+i*2]=new TH1D(Form("%s_j0_noPr",prMass_no[i].data()),Form("%s_j0_noPr",prMass_no[i].data()),200,0,200);
		th5[223+i*2]=new TH1D(Form("%s_j1_noPr",prMass_no[i].data()),Form("%s_j1_noPr",prMass_no[i].data()),200,0,200);
		th5[230+i*2]=new TH1D(Form("%s_j0_noTau21",prMass_no[i].data()),Form("%s_j0_noTau21",prMass_no[i].data()),15,90,150);
		th5[231+i*2]=new TH1D(Form("%s_j1_noTau21",prMass_no[i].data()),Form("%s_j1_noTau21",prMass_no[i].data()),15,90,150);
	}
	th5[265]=new TH1D(Form("%s_j0_noPr",prMass_no[4].data()),Form("%s_j0_noPr",prMass_no[4].data()),200,0,200);
	th5[266]=new TH1D(Form("%s_j1_noPr",prMass_no[4].data()),Form("%s_j1_noPr",prMass_no[4].data()),200,0,200);
	for (int i=0;i<2;i++){
		th5[238+i*2]=new TH1D(Form("%s_j0_noPr_noTau21",tau21_no[i].data()),Form("%s_j0_noPr_noTau21",tau21_no[i].data()),25,0,1);
		th5[239+i*2]=new TH1D(Form("%s_j1_noPr_noTau21",tau21_no[i].data()),Form("%s_j1_noPr_noTau21",tau21_no[i].data()),25,0,1);
		th5[242+i*2]=new TH1D(Form("%s_j0_noPr",tau21_no[i].data()),Form("%s_j0_noPr",tau21_no[i].data()),25,0,1);
		th5[243+i*2]=new TH1D(Form("%s_j1_noPr",tau21_no[i].data()),Form("%s_j1_noPr",tau21_no[i].data()),25,0,1);
		th5[246+i*2]=new TH1D(Form("%s_j0_noTau21",tau21_no[i].data()),Form("%s_j0_noTau21",tau21_no[i].data()),25,0,1);
		th5[247+i*2]=new TH1D(Form("%s_j1_noTau21",tau21_no[i].data()),Form("%s_j1_noTau21",tau21_no[i].data()),25,0,1);
	}
	for(int i=0;i<267;i++){
		th5[i]->Sumw2();
		th_flavor[0][i]=(TH1D* )th5[i]->Clone(Form("%s_bb",th5[i]->GetTitle()));
		th_flavor[1][i]=(TH1D* )th5[i]->Clone(Form("%s_b",th5[i]->GetTitle()));
		th_flavor[2][i]=(TH1D* )th5[i]->Clone(Form("%s_cc",th5[i]->GetTitle()));
		th_flavor[3][i]=(TH1D* )th5[i]->Clone(Form("%s_udcsg",th5[i]->GetTitle()));
		//th_flavor[0][i]->Sumw2();
		//th_flavor[1][i]->Sumw2();
		//th_flavor[2][i]->Sumw2();
		//th_flavor[3][i]->Sumw2();
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
		cout<<w<<endl;
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
		
		
		
		for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){//event loop----------------------------------------------------------------------------------------
			data.GetEntry(jEntry);
			
			
			double PU_weight[3]={1,1,1};
			Float_t ntrue= data.GetFloat("pu_nTrueInt");
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
			fixScaleNum[0]+=PU_weight[0];
			
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
				if( ((thisTrig.find("HLT_PFHT800")!= std::string::npos
						) && results==1)){
					passTrigger=true;
					break;
				}
			}
			if(!passTrigger && nameRoot==2)continue;
			nPass[1]++;

			const int nAK8Jet=data.GetInt("AK8PuppinJet");
			//2.nJets
			if(nAK8Jet<2)continue;nPass[2]++;
			int* AK8PuppinSubSDJet=data.GetPtrInt("AK8PuppinSubSDJet");
			if(AK8PuppinSubSDJet[0]!=2||AK8PuppinSubSDJet[1]!=2)continue;
			TClonesArray* AK8PuppijetP4 = (TClonesArray*) data.GetPtrTObject("AK8PuppijetP4");
			float*  AK8PuppijetCorrUncUp = data.GetPtrFloat("AK8PuppijetCorrUncUp"); 
			float*  AK8PuppijetCorrUncDown = data.GetPtrFloat("AK8PuppijetCorrUncDown"); 
			TLorentzVector* thisJet ,* thatJet;
			
			if(JESOption==0){
				thisJet=(TLorentzVector*)AK8PuppijetP4->At(0);
				thatJet=(TLorentzVector*)AK8PuppijetP4->At(1);
			}
			else if (JESOption==1){
				TLorentzVector thisJetScaleUp= (*((TLorentzVector*)AK8PuppijetP4->At(0)))*(1+AK8PuppijetCorrUncUp[0]);
				TLorentzVector thatJetScaleUp= (*((TLorentzVector*)AK8PuppijetP4->At(1)))*(1+AK8PuppijetCorrUncUp[1]);
				thisJet= &thisJetScaleUp;
				thatJet= &thatJetScaleUp;
			}
			else if (JESOption==2){
				TLorentzVector thisJetScaleDown= (*((TLorentzVector*)AK8PuppijetP4->At(0)))*(1-AK8PuppijetCorrUncDown[0]);
				TLorentzVector thatJetScaleDown= (*((TLorentzVector*)AK8PuppijetP4->At(1)))*(1-AK8PuppijetCorrUncDown[1]);
				thisJet= &thisJetScaleDown;
				thatJet= &thatJetScaleDown;
			}
			
			//3. Pt 
			if(thisJet->Pt()>99998 ||thatJet->Pt()>99998 )continue;
			if(thisJet->Pt()<300)continue;
			if(thatJet->Pt()<300)continue;
			nPass[3]++;
			//4tightId-----------------------------------------
			vector<bool>    &AK8PuppijetPassIDTight = *((vector<bool>*) data.GetPtr("AK8PuppijetPassIDTight"));
			if(AK8PuppijetPassIDTight[0]==0)continue;
			if(AK8PuppijetPassIDTight[1]==0)continue;
			Float_t*  AK8PuppijetCEmEF = data.GetPtrFloat("AK8PuppijetCEmEF");
			Float_t*  AK8PuppijetMuoEF = data.GetPtrFloat("AK8PuppijetMuoEF");
			if(AK8PuppijetMuoEF[0]>0.8)continue;
			if(AK8PuppijetCEmEF[0]>0.9)continue;
			if(AK8PuppijetMuoEF[1]>0.8)continue;
			if(AK8PuppijetCEmEF[1]>0.9)continue;
			nPass[4]++;
			//5. Eta-----------------------------------------
			if(fabs(thisJet->Eta())>2.4)continue;
			if(fabs(thatJet->Eta())>2.4)continue;
			nPass[5]+=PU_weight[0];
			
			int event_flavor=-1;
			Int_t* AK8PuppijetHadronFlavor        = data.GetPtrInt("AK8PuppijetHadronFlavor");
			vector<Int_t>   *AK8PuppisubjetSDHadronFlavor =  data.GetPtrVectorInt("AK8PuppisubjetSDHadronFlavor");
			if(AK8PuppijetHadronFlavor[0]==5 && AK8PuppisubjetSDHadronFlavor[0][0]==5 && AK8PuppisubjetSDHadronFlavor[0][1]==5)event_flavor=0;
			else if(AK8PuppijetHadronFlavor[1]==5 && AK8PuppisubjetSDHadronFlavor[1][0]==5 && AK8PuppisubjetSDHadronFlavor[1][1]==5)event_flavor=0;
			else if(AK8PuppijetHadronFlavor[0]==5 && (AK8PuppisubjetSDHadronFlavor[0][0]==5 || AK8PuppisubjetSDHadronFlavor[0][1]==5))event_flavor=1;
			else if(AK8PuppijetHadronFlavor[1]==5 && (AK8PuppisubjetSDHadronFlavor[1][0]==5 || AK8PuppisubjetSDHadronFlavor[1][1]==5))event_flavor=1;
			else if(AK8PuppijetHadronFlavor[0]==4 && AK8PuppisubjetSDHadronFlavor[0][0]==4 && AK8PuppisubjetSDHadronFlavor[0][1]==4)event_flavor=2;
			else if(AK8PuppijetHadronFlavor[1]==4 && AK8PuppisubjetSDHadronFlavor[1][0]==4 && AK8PuppisubjetSDHadronFlavor[1][1]==4)event_flavor=2;
			else if(AK8PuppijetHadronFlavor[0]==4 && (AK8PuppisubjetSDHadronFlavor[0][0]==4 || AK8PuppisubjetSDHadronFlavor[0][1]==4))event_flavor=2;
			else if(AK8PuppijetHadronFlavor[1]==4 && (AK8PuppisubjetSDHadronFlavor[1][0]==4 || AK8PuppisubjetSDHadronFlavor[1][1]==4))event_flavor=2;
			else event_flavor=3;
			
			th5[180]->Fill(nVtx,PU_weight[0]);
			th5[181]->Fill(ntrue,PU_weight[0]);
			
			th_flavor[event_flavor][180]->Fill(nVtx,PU_weight[0]);
			th_flavor[event_flavor][181]->Fill(ntrue,PU_weight[0]);
			//6. DEta-----------------------------------------
			float dEta = fabs(thisJet->Eta()-thatJet->Eta());
			if(dEta>1.3)continue;
			nPass[6]+=PU_weight[0];
			//7. Mjj-----------------------------------------
			float mjj = (*thisJet+*thatJet).M();
			float mjjRed = (*thisJet+*thatJet).M()+250-thisJet->M()-thatJet->M();
			if(mjjRed<1000)continue;
			nPass[7]+=PU_weight[0];
			//8. fatjetPRmassL2L3Corr-----------------------------------------
			Float_t*  fatjetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
			Float_t*  FATjetPuppiSDmassL2L3Corr = data.GetPtrFloat("FATjetPuppiSDmassL2L3Corr");
			Float_t*  FATjetPRmass = data.GetPtrFloat("FATjetPRmass");
			Float_t*  AK8PuppijetSDmass = data.GetPtrFloat("AK8PuppijetSDmass");
			
			
			
			vector<float>   *subjetSDCSV =  data.GetPtrVectorFloat("AK8PuppisubjetSDCSV");
			
			vector<float>   *subjetSDPx  =  data.GetPtrVectorFloat("AK8PuppisubjetSDPx");
			vector<float>   *subjetSDPy  =  data.GetPtrVectorFloat("AK8PuppisubjetSDPy");
			vector<float>   *subjetSDPz  =  data.GetPtrVectorFloat("AK8PuppisubjetSDPz");
			vector<float>   *subjetSDE   =  data.GetPtrVectorFloat("AK8PuppisubjetSDE");
			
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
					subjetEta[i][j]=subjetP4[i][j]->Eta();
					if(subjetSDCSV[i][j]>0.46)nbtag++;
				}
				dr[i]=subjetP4[i][0]->DeltaR(*subjetP4[i][1]);
			}
			
			double scaleFactor=PU_weight[0];
			
			Float_t  FATjetPuppiSDmassThea[2] ={0};
			FATjetPuppiSDmassThea[0]=AK8PuppijetSDmass[0]*getPUPPIweight(thisJet->Pt(),thisJet->Eta());
			FATjetPuppiSDmassThea[1]=AK8PuppijetSDmass[1]*getPUPPIweight(thatJet->Pt(),thatJet->Eta());
			eta[0]=thisJet->Eta();
			eta[1]=thatJet->Eta();
			Float_t*  AK8Puppijet_DoubleSV = data.GetPtrFloat("AK8Puppijet_DoubleSV");
			
			
			Float_t*  fatjetTau1 = data.GetPtrFloat("FATjetTau1");
			Float_t*  fatjetTau2 = data.GetPtrFloat("FATjetTau2");
			Float_t*  AK8PuppijetTau1 = data.GetPtrFloat("AK8PuppijetTau1");
			Float_t*  AK8PuppijetTau2 = data.GetPtrFloat("AK8PuppijetTau2");
			double tau21[2];
			tau21[0]=(fatjetTau2[0]/fatjetTau1[0]),tau21[1]=(fatjetTau2[1]/fatjetTau1[1]);
			double puppiTau21[2];
			puppiTau21[0]=(AK8PuppijetTau2[0]/AK8PuppijetTau1[0]),puppiTau21[1]=(AK8PuppijetTau2[1]/AK8PuppijetTau1[1]);
			
			for (int i=0;i<2;i++){
				th5[214+i]->Fill(FATjetPRmass[i],scaleFactor);
				th5[216+i]->Fill(fatjetPRmassL2L3Corr[i],scaleFactor);
				th5[218+i]->Fill(AK8PuppijetSDmass[i],scaleFactor);
				th5[220+i]->Fill(FATjetPuppiSDmassL2L3Corr[i],scaleFactor);
				th5[238+i]->Fill(tau21[i],scaleFactor);
				th5[240+i]->Fill(puppiTau21[i],scaleFactor);
				
				th_flavor[event_flavor][214+i]->Fill(FATjetPRmass[i],scaleFactor);
				th_flavor[event_flavor][216+i]->Fill(fatjetPRmassL2L3Corr[i],scaleFactor);
				th_flavor[event_flavor][218+i]->Fill(AK8PuppijetSDmass[i],scaleFactor);
				th_flavor[event_flavor][220+i]->Fill(FATjetPuppiSDmassL2L3Corr[i],scaleFactor);
				th_flavor[event_flavor][238+i]->Fill(tau21[i],scaleFactor);
				th_flavor[event_flavor][240+i]->Fill(puppiTau21[i],scaleFactor);
			}
			/*
			if(!(puppiTau21[0]>0.6 || puppiTau21[1]>0.6)){
				for (int i=0;i<2;i++){
					
					th5[222+i]->Fill(FATjetPRmass[i],scaleFactor);
					th5[224+i]->Fill(fatjetPRmassL2L3Corr[i],scaleFactor);
					th5[226+i]->Fill(AK8PuppijetSDmass[i],scaleFactor);
					th5[228+i]->Fill(FATjetPuppiSDmassL2L3Corr[i],scaleFactor);
					th5[242+i]->Fill(tau21[i],scaleFactor);
					th5[244+i]->Fill(puppiTau21[i],scaleFactor);
					th5[265+i]->Fill(FATjetPuppiSDmassThea[i],scaleFactor);
					
					th_flavor[event_flavor][222+i]->Fill(FATjetPRmass[i],scaleFactor);
					th_flavor[event_flavor][224+i]->Fill(fatjetPRmassL2L3Corr[i],scaleFactor);
					th_flavor[event_flavor][226+i]->Fill(AK8PuppijetSDmass[i],scaleFactor);
					th_flavor[event_flavor][228+i]->Fill(FATjetPuppiSDmassL2L3Corr[i],scaleFactor);
					th_flavor[event_flavor][242+i]->Fill(tau21[i],scaleFactor);
					th_flavor[event_flavor][244+i]->Fill(puppiTau21[i],scaleFactor);
					th_flavor[event_flavor][265+i]->Fill(FATjetPuppiSDmassThea[i],scaleFactor);
				}
			}
			*/
		
		
			if(AK8PuppijetSDmass[0]<50)continue;
			if(AK8PuppijetSDmass[1]<50)continue;
			
			nPass[8]+=PU_weight[0];
			for (int i=0;i<2;i++){
				
				th5[230+i]->Fill(FATjetPRmass[i],scaleFactor);
				th5[232+i]->Fill(fatjetPRmassL2L3Corr[i],scaleFactor);
				th5[234+i]->Fill(AK8PuppijetSDmass[i],scaleFactor);
				th5[236+i]->Fill(FATjetPuppiSDmassL2L3Corr[i],scaleFactor);
				th5[246+i]->Fill(tau21[i],scaleFactor);
				th5[248+i]->Fill(puppiTau21[i],scaleFactor);
				
				th_flavor[event_flavor][230+i]->Fill(FATjetPRmass[i],scaleFactor);
				th_flavor[event_flavor][232+i]->Fill(fatjetPRmassL2L3Corr[i],scaleFactor);
				th_flavor[event_flavor][234+i]->Fill(AK8PuppijetSDmass[i],scaleFactor);
				th_flavor[event_flavor][236+i]->Fill(FATjetPuppiSDmassL2L3Corr[i],scaleFactor);
				th_flavor[event_flavor][246+i]->Fill(tau21[i],scaleFactor);
				th_flavor[event_flavor][248+i]->Fill(puppiTau21[i],scaleFactor);
			}
			
			//9.-----------------------------------------
			//if(puppiTau21[0]>0.6 || puppiTau21[1]>0.6) continue;
			nPass[9]+=PU_weight[0];
			
			
			
			
			
			//10.btag
			
			
			//uncertainty -------------------------------------
			//double scaleFactor=btaggingscaleFactor*PU_weight[0]*tau21_SF;
			
			for(int i=0;i<3;i++){
				passPileup[i]+=btaggingscaleFactor*PU_weight[i];
				th7[i]->Fill(mjj,btaggingscaleFactor*PU_weight[i]);
			}
			 for(int i=0;i<101;i++)passPDF[i]+=btaggingscaleFactor*PU_weight[0];
			 for(int i=5;i<14;i++)th7[i]->Fill(mjj,btaggingscaleFactor*PU_weight[0]);
			 fixScaleNum[1]+=PU_weight[0];
			//--------------------------------------
			
			
			
			
			
			Float_t*  AK8PuppijetCISVV2 = data.GetPtrFloat("AK8PuppijetCISVV2");
			for(int i=0;i<2;i++){
				for(int j=0;j<2;j++){
					
					th5[(i*2+j)*5+nbtag]->Fill(subjetPt[i][j],scaleFactor);
					th5[(i*2+j)*5+nbtag+20]->Fill(subjetEta[i][j],scaleFactor);
					th5[(i*2+j)*5+nbtag+85]->Fill(subjetSDCSV[i][j],scaleFactor);
					if(subjetPt[i][j]>30 && fabs(subjetEta[i][j])<2.4)th5[(i*2+j)*5+nbtag+194]->Fill(subjetSDCSV[i][j],scaleFactor);
					
					th_flavor[event_flavor][(i*2+j)*5+nbtag]->Fill(subjetPt[i][j],scaleFactor);
					th_flavor[event_flavor][(i*2+j)*5+nbtag+20]->Fill(subjetEta[i][j],scaleFactor);
					th_flavor[event_flavor][(i*2+j)*5+nbtag+85]->Fill(subjetSDCSV[i][j],scaleFactor);
					if(subjetPt[i][j]>30 && fabs(subjetEta[i][j])<2.4)th_flavor[event_flavor][(i*2+j)*5+nbtag+194]->Fill(subjetSDCSV[i][j],scaleFactor);
					
				}
				th5[i*5+nbtag+40]->Fill(dr[i],scaleFactor);
				th5[i*5+nbtag+50]->Fill(pt[i],scaleFactor);
				th5[i*5+nbtag+60]->Fill(eta[i],scaleFactor);
				th5[i*5+nbtag+70]->Fill(fatjetPRmassL2L3Corr[i],scaleFactor);
				th5[i*5+nbtag+105]->Fill(tau21[i],scaleFactor);
				th5[i*5+nbtag+120]->Fill(FATjetPuppiSDmassL2L3Corr[i],scaleFactor);
				th5[i*5+nbtag+130]->Fill(puppiTau21[i],scaleFactor);
				th5[i*5+nbtag+140]->Fill(FATjetPRmass[i],scaleFactor);
				th5[i*5+nbtag+150]->Fill(AK8PuppijetSDmass[i],scaleFactor);
				th5[i*5+nbtag+170]->Fill(AK8Puppijet_DoubleSV[i],PU_weight[0]);
				th5[i*5+nbtag+184]->Fill(AK8PuppijetCISVV2[i],scaleFactor);
				th5[i*5+nbtag+250]->Fill(FATjetPuppiSDmassThea[i],scaleFactor);
				
				th_flavor[event_flavor][i*5+nbtag+40]->Fill(dr[i],scaleFactor);
				th_flavor[event_flavor][i*5+nbtag+50]->Fill(pt[i],scaleFactor);
				th_flavor[event_flavor][i*5+nbtag+60]->Fill(eta[i],scaleFactor);
				th_flavor[event_flavor][i*5+nbtag+70]->Fill(fatjetPRmassL2L3Corr[i],scaleFactor);
				th_flavor[event_flavor][i*5+nbtag+105]->Fill(tau21[i],scaleFactor);
				th_flavor[event_flavor][i*5+nbtag+120]->Fill(FATjetPuppiSDmassL2L3Corr[i],scaleFactor);
				th_flavor[event_flavor][i*5+nbtag+130]->Fill(puppiTau21[i],scaleFactor);
				th_flavor[event_flavor][i*5+nbtag+140]->Fill(FATjetPRmass[i],scaleFactor);
				th_flavor[event_flavor][i*5+nbtag+150]->Fill(AK8PuppijetSDmass[i],scaleFactor);
				th_flavor[event_flavor][i*5+nbtag+170]->Fill(AK8Puppijet_DoubleSV[i],PU_weight[0]);
				th_flavor[event_flavor][i*5+nbtag+184]->Fill(AK8PuppijetCISVV2[i],scaleFactor);
				th_flavor[event_flavor][i*5+nbtag+250]->Fill(FATjetPuppiSDmassThea[i],scaleFactor);
				
			}
			th5[nbtag+80]->Fill(mjj,scaleFactor);
			th5[nbtag+115]->Fill(dEta,scaleFactor);
			th5[nbtag+165]->Fill(mjjRed,scaleFactor);
			
			Float_t HT =data.GetFloat("HT");
			
			th5[nbtag+260]->Fill(HT,scaleFactor);
			th_flavor[event_flavor][nbtag+260]->Fill(HT,scaleFactor);
			
			th_flavor[event_flavor][nbtag+80]->Fill(mjj,scaleFactor);
			th_flavor[event_flavor][nbtag+115]->Fill(dEta,scaleFactor);
			th_flavor[event_flavor][nbtag+165]->Fill(mjjRed,scaleFactor);
			
			nPass[nbtag+10]+=PU_weight[0]*btaggingscaleFactor;
				
			
			
			int nDoubleSV=0;
			if(AK8Puppijet_DoubleSV[0]>0.6)nDoubleSV++;
			if(AK8Puppijet_DoubleSV[1]>0.6)nDoubleSV++;
			
			if(nDoubleSV==0)nPass[15]+=PU_weight[0];
			else if(nDoubleSV==2)nPass[17]+=PU_weight[0];
			else {
				nPass[16]+=PU_weight[0];
				if((AK8Puppijet_DoubleSV[0]>0.6 && subjetSDCSV[1][0]<0.46 && subjetSDCSV[1][1]<0.46) ||
				     (AK8Puppijet_DoubleSV[1]>0.6 && subjetSDCSV[0][0]<0.46 && subjetSDCSV[0][1]<0.46))nPass[18]+=PU_weight[0];
				else if((AK8Puppijet_DoubleSV[0]>0.6 && subjetSDCSV[1][0]>0.46 && subjetSDCSV[1][1]>0.46) ||
				               (AK8Puppijet_DoubleSV[1]>0.6 && subjetSDCSV[0][0]>0.46 && subjetSDCSV[0][1]>0.46))nPass[20]+=PU_weight[0];
				else nPass[19]+=PU_weight[0];
			}
			
			
			th5[182]->Fill(nVtx,scaleFactor);
			th5[183]->Fill(ntrue,scaleFactor);
			th_flavor[event_flavor][182]->Fill(nVtx,scaleFactor);
			th_flavor[event_flavor][183]->Fill(ntrue,scaleFactor);

		}//end event loop----------------------------------------------------------------------------------------
	}	//end ntuple loop----------------------------------------------------------------------------------------
	cout<<"entries="<<total<<endl;	
	for(int i=0;i<22;i++)cout<<"nPass["<<i<<"]="<<nPass[i]<<endl;
	
	for(int i=0;i<3;i++)th7[3]->SetBinContent(i+1,passPileup[i]/totalPileup[i]);
	for(int i=0;i<101;i++)th7[4]->SetBinContent(i+1,passPDF[i]/totalPDF[i]);
	
	TH1D * th2o=new TH1D("Nbtagjet","Nbtagjet",5,-0.5,4.5);
	for (int i=0;i<5;i++){
		if(nameRoot==2 && i>2)continue;
		th2o->SetBinContent(i+1,nPass[i+10]);
	}
	
	TH1D * fixScale=new TH1D("fixScale","fixScale",2,-0.5,1.5);
	fixScale->SetBinContent(1,fixScaleNum[0]);
	fixScale->SetBinContent(2,fixScaleNum[1]);
	TH1D * cutflow=new TH1D("cutflow","cutflow",21,0.5,21.5);
	cutflow->SetBinContent(1,fixScaleNum[0]);
	if(nameRoot!=2)for(int ii=1;ii<22;ii++)cutflow->SetBinContent(ii+1,nPass[ii-1]);
	else for(int ii=1;ii<16;ii++)cutflow->SetBinContent(ii+1,nPass[ii-1]);
	
	TFile* outFile ;
	if(nameRoot==1)outFile= new TFile(Form("sf2/%s.root",st2.data()),"recreate");
	else if(JESOption==0)outFile= new TFile(Form("sf2/%s/%d.root",st2.data(),wMs),"recreate");
	else if(JESOption==1)outFile= new TFile(Form("sf2/%s_JESUp.root",st2.data()),"recreate");
	else if(JESOption==2)outFile= new TFile(Form("sf2/%s_JESDown.root",st2.data()),"recreate");
	else if(JESOption==3)outFile= new TFile(Form("sf2/%s_BtagUp.root",st2.data()),"recreate");
	else if(JESOption==4)outFile= new TFile(Form("sf2/%s_BtagDown.root",st2.data()),"recreate");
	else if(JESOption==5)outFile= new TFile(Form("sf2/%s_tau21Up.root",st2.data()),"recreate");
	else if(JESOption==6)outFile= new TFile(Form("sf2/%s_tau21Down.root",st2.data()),"recreate");
	th2o->Write();
	fixScale->Write();
	cutflow->Write();
	
	for(int i=0;i<267;i++){
		th5[i]->Write();
		th_flavor[0][i]->Write();
		th_flavor[1][i]->Write();
		th_flavor[2][i]->Write();
		th_flavor[3][i]->Write();
	}
	
	for(int i=0;i<14;i++)th7[i]->Write();
	outFile->Close();
}