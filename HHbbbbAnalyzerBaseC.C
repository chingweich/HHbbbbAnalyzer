void HHbbbbAnalyzerBaseC(int wMs,int wM, string st,string st2,double Xsec,int nameRoot=0){	
	TFile *f;
	TTree *tree;
	int nPass[20]={0};int total=0;
	int dataPassingcsc=0;
	
	TH1D * th1[50];
	
	th1[0]=new TH1D("Pt_j0_sj0_0b","Pt_j0_sj0_0b",200,0,2000);
	th1[1]=new TH1D("Pt_j0_sj1_0b","Pt_j0_sj1_0b",200,0,2000);
	th1[2]=new TH1D("Pt_j1_sj0_0b","Pt_j1_sj0_0b",200,0,2000);
	th1[3]=new TH1D("Pt_j1_sj1_0b","Pt_j1_sj1_0b",200,0,2000);
	
	th1[4]=new TH1D("Pt_j0_sj0_1b","Pt_j0_sj0_1b",200,0,2000);
	th1[5]=new TH1D("Pt_j0_sj1_1b","Pt_j0_sj1_1b",200,0,2000);
	th1[6]=new TH1D("Pt_j1_sj0_1b","Pt_j1_sj0_1b",200,0,2000);
	th1[7]=new TH1D("Pt_j1_sj1_1b","Pt_j1_sj1_1b",200,0,2000);
	
	th1[8]=new TH1D("Pt_j0_sj0_2b","Pt_j0_sj0_2b",200,0,2000);
	th1[9]=new TH1D("Pt_j0_sj1_2b","Pt_j0_sj1_2b",200,0,2000);
	th1[10]=new TH1D("Pt_j1_sj0_2b","Pt_j1_sj0_2b",200,0,2000);
	th1[11]=new TH1D("Pt_j1_sj1_2b","Pt_j1_sj1_2b",200,0,2000);
	
	th1[12]=new TH1D("deltaR0_0b","deltaR0_0b",20,0,1);
	th1[13]=new TH1D("deltaR1_0b","deltaR1_0b",20,0,1);
	th1[14]=new TH1D("deltaR0_1b","deltaR0_1b",20,0,1);
	th1[15]=new TH1D("deltaR1_1b","deltaR1_1b",20,0,1);
	th1[16]=new TH1D("deltaR0_2b","deltaR0_2b",20,0,1);
	th1[17]=new TH1D("deltaR1_2b","deltaR1_2b",20,0,1);
	
	th1[18]=new TH1D("pt0","pt0",480,200,1400);
	th1[19]=new TH1D("pt1","pt1",480,200,1400);
	
	th1[20]=new TH1D("DeltaEta","DeltaEta",16,-0.1,1.5);
	th1[21]=new TH1D("eta_j0","eta_j0",60,-3,3);
	th1[22]=new TH1D("eta_j1","eta_j1",60,-3,3);
	th1[23]=new TH1D("HT12","HT12",80,600,2600);
	th1[24]=new TH1D("prMass_j0","prMass_j0",12,90,150);
	th1[25]=new TH1D("prMassCut_j0","prMassCut_j0",12,90,150);
	th1[26]=new TH1D("prMass_j1","prMass_j1",12,90,150);
	th1[27]=new TH1D("prMassCut_j1","prMassCut_j1",12,90,150);
	th1[28]=new TH1D("tau21_j0","tau21_j0",10,0,1);
	th1[29]=new TH1D("tau21_j1","tau21_j1",10,0,1);
	
	th1[30]=new TH1D("totalMass_0b","totalMass_0b",75,1000,4000);
	th1[31]=new TH1D("totalMass_1b","totalMass_1b",75,1000,4000);
	th1[32]=new TH1D("totalMass_2b","totalMass_2b",75,1000,4000);
	th1[33]=new TH1D("totalMass_3b","totalMass_3b",75,1000,4000);
	th1[34]=new TH1D("totalMass_4b","totalMass_4b",75,1000,4000);
	th1[35]=new TH1D("totalMass_all","totalMass_all",75,1000,4000);
	
	th1[36]=new TH1D("CSV_j0_sj0","CSV_j0_sj0",25,0,1);
	th1[37]=new TH1D("CSV_j0_sj1","CSV_j0_sj1",25,0,1);
	th1[38]=new TH1D("CSV_j1_sj0","CSV_j1_sj0",25,0,1);
	th1[39]=new TH1D("CSV_j1_sj1","CSV_j1_sj1",25,0,1);
	
	
	for(int i=0;i<40;i++)th1[i]->Sumw2();
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
			th1[24]->Fill(fatjetPRmassL2L3Corr[0]);
			th1[26]->Fill(fatjetPRmassL2L3Corr[1]);
			if(fatjetPRmassL2L3Corr[0]<105||fatjetPRmassL2L3Corr[0]>135)continue;
			if(fatjetPRmassL2L3Corr[1]<105||fatjetPRmassL2L3Corr[1]>135)continue;
			nPass[6]++;
			//7.-----------------------------------------
			double tau21_1=(fatjetTau2[0]/fatjetTau1[0]),
		           tau21_2=(fatjetTau2[1]/fatjetTau1[1]);
			th1[28]->Fill(tau21_1);	   
			th1[29]->Fill(tau21_2);	   
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
			
			th1[18]->Fill(thisJet->Pt());
			th1[19]->Fill(thatJet->Pt());
			th1[23]->Fill(thisJet->Pt()+thatJet->Pt());
			th1[20]->Fill(dEta);
			th1[21]->Fill(thisJet->Eta());
			th1[22]->Fill(thatJet->Eta());
			th1[25]->Fill(fatjetPRmassL2L3Corr[0]);
			th1[27]->Fill(fatjetPRmassL2L3Corr[1]);
			
			//8.btag
			int nbtag=0;
			if(subjetSDCSV[0][0]>0.605)nbtag++;
			if(subjetSDCSV[0][1]>0.605)nbtag++;
			if(subjetSDCSV[1][0]>0.605)nbtag++;
			if(subjetSDCSV[1][1]>0.605)nbtag++;
			
			th1[36]->Fill(subjetSDCSV[0][0]);
			th1[37]->Fill(subjetSDCSV[0][1]);
			th1[38]->Fill(subjetSDCSV[1][0]);
			th1[39]->Fill(subjetSDCSV[1][1]);
			
	
			TLorentzVector* thisSub1=new TLorentzVector(0,0,0,0);
			TLorentzVector* thisSub2=new TLorentzVector(0,0,0,0);
			TLorentzVector* thatSub1=new TLorentzVector(0,0,0,0);
			TLorentzVector* thatSub2=new TLorentzVector(0,0,0,0);
			thisSub1->SetPxPyPzE(subjetSDPx[0][0],subjetSDPy[0][0],subjetSDPz[0][0],subjetSDE[0][0]);
			thisSub2->SetPxPyPzE(subjetSDPx[0][1],subjetSDPy[0][1],subjetSDPz[0][1],subjetSDE[0][1]);
			thatSub1->SetPxPyPzE(subjetSDPx[1][0],subjetSDPy[1][0],subjetSDPz[1][0],subjetSDE[1][0]);
			thatSub2->SetPxPyPzE(subjetSDPx[1][1],subjetSDPy[1][1],subjetSDPz[1][1],subjetSDE[1][1]);
			double dr1=thisSub1->DeltaR(*thisSub2),dr2=thatSub1->DeltaR(*thatSub2);
			
			th1[35]->Fill((*thisJet+*thatJet).M());
			if(nbtag==0){
				nPass[8]++;
				th1[0]->Fill(thisSub1->Pt());
				th1[1]->Fill(thisSub2->Pt());
				th1[2]->Fill(thatSub1->Pt());
				th1[3]->Fill(thatSub2->Pt());
				th1[12]->Fill(dr1);
				th1[13]->Fill(dr2);
				th1[30]->Fill((*thisJet+*thatJet).M());
				
			}
			if(nbtag==1){
				nPass[9]++;
				th1[4]->Fill(thisSub1->Pt());
				th1[5]->Fill(thisSub2->Pt());
				th1[6]->Fill(thatSub1->Pt());
				th1[7]->Fill(thatSub2->Pt());
				th1[14]->Fill(dr1);
				th1[15]->Fill(dr2);
				th1[31]->Fill((*thisJet+*thatJet).M());
				
			}
			if(nbtag==2){
				nPass[10]++;
				th1[8]->Fill(thisSub1->Pt());
				th1[9]->Fill(thisSub2->Pt());
				th1[10]->Fill(thatSub1->Pt());
				th1[11]->Fill(thatSub2->Pt());
				th1[16]->Fill(dr1);
				th1[17]->Fill(dr2);
				th1[32]->Fill((*thisJet+*thatJet).M());
				
			}
			if(nbtag==3){
				nPass[11]++;
				th1[33]->Fill((*thisJet+*thatJet).M());
				if(isHPHP)nPass[13]++;
			}
			if(nbtag==4){
				nPass[12]++;
				th1[34]->Fill((*thisJet+*thatJet).M());
			}
		}
	}	
	cout<<"entries="<<total<<endl;	
	if(nameRoot==2)cout<<"dataPassingcsc="<<dataPassingcsc<<endl;
	for(int i=0;i<14;i++)cout<<"nPass["<<i<<"]="<<nPass[i]<<endl;
	
	TH1D * th2=new TH1D("Nbtagjet","Nbtagjet",5,-0.5,4.5);
	th2->SetBinContent(1,nPass[8]);
	th2->SetBinContent(2,nPass[9]);
	th2->SetBinContent(3,nPass[10]);
	if(nameRoot!=2)th2->SetBinContent(4,nPass[11]);
	if(nameRoot!=2)th2->SetBinContent(5,nPass[12]);
	
	TH1D * th2s=new TH1D("NbtagjetS","NbtagjetS",5,-0.5,4.5);
	th2s->SetBinContent(1,nPass[8]);
	th2s->SetBinContent(2,nPass[9]);
	th2s->SetBinContent(3,nPass[10]);
	if(nameRoot!=2)th2s->SetBinContent(4,nPass[11]);
	if(nameRoot!=2)th2s->SetBinContent(5,nPass[12]);
	
	TH1D * th4=new TH1D("cutflow","cutflow",15,0.5,15.5);
	th4->SetBinContent(1,total);
	if(nameRoot==2)for(int ii=1;ii<12;ii++)th4->SetBinContent(ii+1,nPass[ii-1]);
	else for(int ii=1;ii<15;ii++)th4->SetBinContent(ii+1,nPass[ii-1]);
	
	TH1D * th4s=new TH1D("cutflowS","cutflowS",15,0.5,15.5);
	th4s->SetBinContent(1,total);
	if(nameRoot==2)for(int ii=1;ii<12;ii++)th4s->SetBinContent(ii+1,nPass[ii-1]);
	else for(int ii=1;ii<15;ii++)th4s->SetBinContent(ii+1,nPass[ii-1]);
	
	TFile* outFile = new TFile(Form("root_files_7476/%s.root",st2.data()),"recreate");
	if(nameRoot==0){
		for(int i=0;i<40;i++){
			//th1[i]->Sumw2();
			th1[i]->Scale(2245.87*Xsec/total);
			
		}
	}
	for(int i=0;i<40;i++)th1[i]->Write();
	th2->Write();
	th4->Write();
	th2s->Scale(2245.87*Xsec/total);
	th2s->Write();
	th4s->Scale(2245.87*Xsec/total);
	th4s->Write();
	outFile->Close();
	
}
