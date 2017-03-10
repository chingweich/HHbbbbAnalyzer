#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include "untuplizer.h"


void TMVARegressionApplication( int wMs,int wM, string st,string st2s ) 
{
		TFile* f;
		f = TFile::Open(Form("dataG/skimmed_NCUGlobalTuples_%d.root",wMs));
		if (!f || !f->IsOpen()) cout<<"missing"<<wMs <<endl;
		
		TTree* tree;
		
		tree=(TTree*)f->Get("treeMaker");
		TreeReader data(tree);
		cout<<data.GetEntriesFast()<<endl;
}

void checkSkimmedJobs(int a,int b){

	string st1[40]={
		"/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/B/skimmed/NCUGlobalTuples_",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/C/skimmed/NCUGlobalTuples_",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/D/skimmed/NCUGlobalTuples_",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/E/skimmed/NCUGlobalTuples_",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/F/skimmed/NCUGlobalTuples_",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/G/skimmed/NCUGlobalTuples_",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/H/skimmed/NCUGlobalTuples_",
	};
	string  fileName[40]={
	"dataB",
	"dataC",
	"dataD",
	"dataE",
	"dataF",
	"dataG",
	"dataH",
	};
	int aa[40]={1045,348,584,497,361,854,24};
	//cout<<st1[a]<<b<<".root"<<endl;
	TMVARegressionApplication(b,b+1,st1[a],fileName[a]);	
	
}
