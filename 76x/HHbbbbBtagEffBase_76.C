void HHbbbbBtagEffBase_76(int wMs,int wM, string st,string st2){	

	//0=signal ,1=QCD ,2=data
	int nameRoot=1;
	if(st2.find("QCD")!= std::string::npos)nameRoot=0;
	if(st2.find("data")!= std::string::npos)nameRoot=2;
	
	TFile *f;
	TTree *tree;
	int nPass[20]={0};int total=0;
	int dataPassingcsc=0;
	double nPassB[6]={0};
	
	TFile *f1;
	if(nameRoot==2)f1=TFile::Open("btagEffSource/data.root");
	else f1=TFile::Open(Form("btagEffSource/%s.root",st2.data()));
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
	
	TH1D * th5[100];
	
	for(int i=0;i<2;i++){
		for(int j=0;j<2;j++){
			for(int k=0;k<5;k++){
				th5[(i*2+j)*5+k]=new TH1D(Form("Pt_j%d_sj%d_%db",i,j,k),Form("Pt_j%d_sj%d_%db",i,j,k),200,0,2000);
				th5[(i*2+j)*5+k+20]=new TH1D(Form("Eta_j%d_sj%d_%db",i,j,k),Form("Eta_j%d_sj%d_%db",i,j,k),60,-3,3);
			}
		}
		for(int k=0;k<5;k++){
			th5[i*5+k+40]=new TH1D(Form("deltaR_j%d_%db",i,k),Form("deltaR_j%d_%db",i,k),20,0,1);
			th5[i*5+k+50]=new TH1D(Form("Pt_j%d_%db",i,k),Form("Pt_j%d_%db",i,k),200,0,2000);
			th5[i*5+k+60]=new TH1D(Form("Eta_j%d_%db",i,k),Form("Eta_j%d_%db",i,k),60,-3,3);
			th5[i*5+k+70]=new TH1D(Form("prMass_j%d_%db",i,k),Form("prMass_j%d_%db",i,k),15,90,150);
		}
	}
	for(int k=0;k<5;k++){
		th5[k+80]=new TH1D(Form("totalMass_%db",k),Form("totalMass_%db",k),200,1000,5000);
	}
	TH1D * th6[100];
	for(int i=0;i<85;i++){
		th6[i]=(TH1D* )th5[i]->Clone(Form("%ss",th5[i]->GetTitle()));
		th5[i]->Sumw2();
		th6[i]->Sumw2();
	}

	TH1D* th7[3];
	th7[0]=new TH1D("totalMass_pileup_up","totalMass_pileup_up",200,1000,5000);
	th7[1]=new TH1D("totalMass_pileup_down","totalMass_pileup_down",200,1000,5000);
	th7[2]=new TH1D("totalMass","totalMass",200,1000,5000);
	
	BTagCalibration calib("CSVv2L", "CSVv2_76.csv");
	BTagCalibrationReader LF(&calib,               // calibration instance
                             BTagEntry::OP_LOOSE,  // operating point
                             "incl",               // measurement type
                             "central");           // systematics type
	
	BTagCalibrationReader HF(&calib,               // calibration instance
                             BTagEntry::OP_LOOSE,  // operating point
                             "mujets",               // measurement type
                             "central");           // systematics type
	
	standalone_LumiReWeighting LumiWeights_central(0);
	standalone_LumiReWeighting LumiWeights_up(1);
	standalone_LumiReWeighting LumiWeights_down(-1);
	
	std::vector<TString> eventlist;                                                                                                                                            
	if(nameRoot==2){
		std::vector<TString> eventlist;                                                                                                                                            
		std::ifstream fin("csc2015.txt");
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
			
			Int_t nVtx        = data.GetInt("nVtx");
			//0. has a good vertex
			if(nVtx<1)continue;nPass[0]++;
			
			//1.trigger
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
			if(!passTrigger)continue;nPass[1]++;
		
			int nFATJet         = data.GetInt("FATnJet");
			const int nJets=nFATJet;
			TClonesArray* fatjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
			Float_t*  fatjetTau1 = data.GetPtrFloat("FATjetTau1");
			Float_t*  fatjetTau2 = data.GetPtrFloat("FATjetTau2");
			Float_t*  fatjetCISVV2 = data.GetPtrFloat("FATjetCISVV2");
			Float_t*  fatjetPRmass = data.GetPtrFloat("FATjetPRmass");
			Float_t*  fatjetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
			Float_t*  fatjetSDmass = data.GetPtrFloat("FATjetSDmass");
			Int_t*   nSubSoftDropJet = data.GetPtrInt("FATnSubSDJet");
			vector<float>   *subjetSDCSV =  data.GetPtrVectorFloat("FATsubjetSDCSV");
			vector<float>   *subjetSDPx  =  data.GetPtrVectorFloat("FATsubjetSDPx", nFATJet);
			vector<float>   *subjetSDPy  =  data.GetPtrVectorFloat("FATsubjetSDPy", nFATJet);
			vector<float>   *subjetSDPz  =  data.GetPtrVectorFloat("FATsubjetSDPz", nFATJet);
			vector<float>   *subjetSDE   =  data.GetPtrVectorFloat("FATsubjetSDE", nFATJet);
			vector<Int_t>   *FATsubjetSDHadronFlavor =  data.GetPtrVectorInt("FATsubjetSDHadronFlavor");
			vector<bool>    &FATjetPassIDTight = *((vector<bool>*) data.GetPtr("FATjetPassIDTight"));
			//2.nJets
			if(nJets<2)continue;nPass[2]++;
			TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(0);
			TLorentzVector* thatJet = (TLorentzVector*)fatjetP4->At(1);
			//3. Pt 
			if(thisJet->Pt()<300)continue;
			if(thatJet->Pt()<300)continue;
			nPass[3]++;
			//4tightId-----------------------------------------
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
			float mjj = (*thisJet+*thatJet).M();
			if(mjj<1000)continue;
			nPass[7]++;
			//8. fatjetPRmassL2L3Corr-----------------------------------------
			if(fatjetPRmassL2L3Corr[0]<105||fatjetPRmassL2L3Corr[0]>135)continue;
			if(fatjetPRmassL2L3Corr[1]<105||fatjetPRmassL2L3Corr[1]>135)continue;
			nPass[8]++;
			//9.-----------------------------------------
			double tau21_1=(fatjetTau2[0]/fatjetTau1[0]),tau21_2=(fatjetTau2[1]/fatjetTau1[1]);
			if(tau21_1>0.75||tau21_2>0.75)continue;
			if(tau21_1>0.6 &&tau21_2>0.6) continue;
			nPass[9]++;
			bool isHPHP=0;
			double tau21_SF=1.031*0.881;
			if(tau21_1<0.6 && tau21_2<0.6 ){
				isHPHP=1;
				tau21_SF=1.031*1.031;
			}
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

			TLorentzVector* thisSub1=new TLorentzVector(0,0,0,0);
			TLorentzVector* thisSub2=new TLorentzVector(0,0,0,0);
			TLorentzVector* thatSub1=new TLorentzVector(0,0,0,0);
			TLorentzVector* thatSub2=new TLorentzVector(0,0,0,0);
			thisSub1->SetPxPyPzE(subjetSDPx[0][0],subjetSDPy[0][0],subjetSDPz[0][0],subjetSDE[0][0]);
			thisSub2->SetPxPyPzE(subjetSDPx[0][1],subjetSDPy[0][1],subjetSDPz[0][1],subjetSDE[0][1]);
			thatSub1->SetPxPyPzE(subjetSDPx[1][0],subjetSDPy[1][0],subjetSDPz[1][0],subjetSDE[1][0]);
			thatSub2->SetPxPyPzE(subjetSDPx[1][1],subjetSDPy[1][1],subjetSDPz[1][1],subjetSDE[1][1]);
			
			float MaxBJetPt = 670., MaxLJetPt = 1000.;
			double sf[2][2];sf[0][0]=1;sf[0][1]=1;sf[1][0]=1;sf[1][1]=1;
			double subjetPt[2][2];subjetPt[0][0]=thisSub1->Pt();subjetPt[0][1]=thisSub2->Pt();subjetPt[1][0]=thatSub1->Pt();subjetPt[1][1]=thatSub2->Pt();
			double subjetEta[2][2];subjetEta[0][0]=thisSub1->Eta();subjetEta[0][1]=thisSub2->Eta();subjetEta[1][0]=thatSub1->Eta();subjetEta[1][1]=thatSub2->Eta();
			double eff[2][2];
			int nbtag=0;int nbtag2=0;
			double btaggingscaleFactor=1;
			
			for(int i=0;i<2;i++){
				for(int j=0;j<2;j++){
					//get btagging eff------------------------------------------------------------
					if(FATsubjetSDHadronFlavor[i][j]==5)eff[i][j]=th1[1]->GetBinContent(ceil(subjetPt[i][j]/10),ceil(subjetEta[i][j]/0.1)+30);
					else if(FATsubjetSDHadronFlavor[i][j]==4)eff[i][j]=th1[3]->GetBinContent(ceil(subjetPt[i][j]/10),ceil(subjetEta[i][j]/0.1)+30);
					else eff[i][j]=th1[5]->GetBinContent(ceil(subjetPt[i][j]/10),ceil(subjetEta[i][j]/0.1)+30);
					//check maxPt-------------------------------------------------------------
					if(FATsubjetSDHadronFlavor[i][j]!=0 && subjetPt[i][j]>MaxBJetPt )subjetPt[i][j]=MaxBJetPt-0.1;
					if(FATsubjetSDHadronFlavor[i][j]==0 && subjetPt[i][j]>MaxLJetPt )subjetPt[i][j]=MaxLJetPt-0.1;
					//Get SF from csv------------------------------------------------------------
					if(FATsubjetSDHadronFlavor[i][j]==5){
						sf[i][j]=HF.eval(BTagEntry::FLAV_B,subjetEta[i][j],subjetPt[i][j]); 
						th3[3]->Fill(subjetPt[i][j],sf[i][j]);
					}
					else if(FATsubjetSDHadronFlavor[i][j]==4){
						sf[i][j]=HF.eval(BTagEntry::FLAV_C,subjetEta[i][j],subjetPt[i][j]); 
						th3[4]->Fill(subjetPt[i][j],sf[i][j]);
					}
					else {
						sf[i][j]=LF.eval(BTagEntry::FLAV_UDSG,subjetEta[i][j],subjetPt[i][j]); 
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
					if(subjetSDCSV[i][j]>0.46 &&FATsubjetSDHadronFlavor[i][j]==0)nbtag2++;\
					if(subjetSDCSV[i][j]>0.46)th4[i*2+j]->Fill(sf[i][j]);
					else th4[i*2+j+4]->Fill(sf[i][j]);
					
				}
			}
			th4[8]->Fill(btaggingscaleFactor);
			if(nbtag2==2 && nbtag==2)th4[12]->Fill(btaggingscaleFactor);
			if(nbtag2==1 && nbtag==2)th4[13]->Fill(btaggingscaleFactor);
			
			double dr[2];dr[0]=thisSub1->DeltaR(*thisSub2),dr[1]=thatSub1->DeltaR(*thatSub2);
			double pt[2];pt[0]=thisJet->Pt(),pt[1]=thatJet->Pt();
			double eta[2];eta[0]=thisJet->Eta(),eta[1]=thatJet->Eta();
			
			Float_t ntrue= data.GetFloat("pu_nTrueInt");
			double PU_weight_central =1, PU_weight_up =1, PU_weight_down =1;
			if(nameRoot!=2){
				PU_weight_central = LumiWeights_central.weight(ntrue);
				PU_weight_up = LumiWeights_up.weight(ntrue);
				PU_weight_down = LumiWeights_down.weight(ntrue);
			}
			double scaleFactor=btaggingscaleFactor*PU_weight_central*tau21_SF;
			th7[0]->Fill(mjj,btaggingscaleFactor*PU_weight_up*tau21_SF);
			th7[1]->Fill(mjj,btaggingscaleFactor*PU_weight_down*tau21_SF);
			th7[2]->Fill(mjj,scaleFactor);
			
			subjetPt[0][0]=thisSub1->Pt();subjetPt[0][1]=thisSub2->Pt();subjetPt[1][0]=thatSub1->Pt();subjetPt[1][1]=thatSub2->Pt();
			subjetEta[0][0]=thisSub1->Eta();subjetEta[0][1]=thisSub2->Eta();subjetEta[1][0]=thatSub1->Eta();subjetEta[1][1]=thatSub2->Eta();
			
			for(int i=0;i<2;i++){
				for(int j=0;j<2;j++){
					for(int k=0;k<5;k++){
						if(nbtag!=k)continue;
						th5[(i*2+j)*5+k]->Fill(subjetPt[i][j]);
						th5[(i*2+j)*5+k+20]->Fill(subjetEta[i][j]);
						th6[(i*2+j)*5+k]->Fill(subjetPt[i][j],scaleFactor);
						th6[(i*2+j)*5+k+20]->Fill(subjetEta[i][j],scaleFactor);
					}
				}
				for(int k=0;k<5;k++){
					if(nbtag!=k)continue;
					th5[i*5+k+40]->Fill(dr[i]);
					th5[i*5+k+50]->Fill(pt[i]);
					th5[i*5+k+60]->Fill(eta[i]);
					th5[i*5+k+70]->Fill(fatjetPRmassL2L3Corr[i]);
					th6[i*5+k+40]->Fill(dr[i],scaleFactor);
					th6[i*5+k+50]->Fill(pt[i],scaleFactor);
					th6[i*5+k+60]->Fill(eta[i],scaleFactor);
					th6[i*5+k+70]->Fill(fatjetPRmassL2L3Corr[i],scaleFactor);
				}
			}
			for(int k=0;k<5;k++){
				if(nbtag!=k)continue;
				th5[k+80]->Fill(mjj);
				th6[k+80]->Fill(mjj,scaleFactor);
				nPassB[k]+=scaleFactor;
				nPass[k+10]++;
				if(k<3)th4[9+k]->Fill(btaggingscaleFactor);
			}
			if(nbtag==3){
				if(isHPHP)nPassB[5]+=scaleFactor;
				if(isHPHP)nPass[15]++;
			}

		}
	}	
	cout<<"entries="<<total<<endl;	
	
	for(int i=0;i<6;i++)cout<<"nPassB["<<i<<"]="<<nPassB[i]<<endl;
	for(int i=0;i<16;i++)cout<<"nPass["<<i<<"]="<<nPass[i]<<endl;
	
	TH1D * th2o=new TH1D("Nbtagjet","Nbtagjet",6,-0.5,5.5);
	th2o->SetBinContent(1,nPassB[0]);
	th2o->SetBinContent(2,nPassB[1]);
	th2o->SetBinContent(3,nPassB[2]);
	if(nameRoot!=2)th2o->SetBinContent(4,nPassB[3]);
	if(nameRoot!=2)th2o->SetBinContent(5,nPassB[4]);
	if(nameRoot!=2)th2o->SetBinContent(6,nPassB[5]);
	
	TH1D * cutflow=new TH1D("cutflow","cutflow",15,0.5,15.5);
	cutflow->SetBinContent(1,total);
	if(nameRoot==2)for(int ii=1;ii<14;ii++)cutflow->SetBinContent(ii+1,nPass[ii-1]);
	else for(int ii=1;ii<17;ii++)cutflow->SetBinContent(ii+1,nPass[ii-1]);
	
	TFile* outFile = new TFile(Form("sf/%s.root",st2.data()),"recreate");
	th2o->Write();
	cutflow->Write();
	for(int i=0;i<6;i++)th3[i]->Write();
	for(int i=0;i<14;i++)th4[i]->Write();
	for(int i=0;i<85;i++){
		th5[i]->Write();
		th6[i]->Write();
	}
	for(int i=0;i<3;i++)th7[i]->Write();
	outFile->Close();
}