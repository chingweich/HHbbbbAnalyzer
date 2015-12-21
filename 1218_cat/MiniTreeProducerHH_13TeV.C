{


  double mgg, mjj,evWeight, mtot, normWeight;
 Int_t categories;

 evWeight = 1.0;
 normWeight = 1;


 // string sInFile = "input/MassPlotFineBins.root";
 string sInFile = "input/MassPlotFineBins_2.root";        
cout << sInFile.c_str() << endl;
 TFile file0(sInFile.c_str(), "read");
 

 string sOutFile("MiniTrees/Background_HH_13TeV/dijetHH_miniTree.root");
 TFile f1(sOutFile.c_str(), "recreate");
 f1.cd();
 
 TTree *TCVARS = new TTree("TCVARS", "HH selection");
 TCVARS->Branch("mgg",&mgg,"mgg/D");
 
 TCVARS->Branch("evWeight",&evWeight,"evWeight/D");
 TCVARS->Branch("normWeight",&normWeight,"normWeight/D");
 
 TCVARS->Branch("categories",&categories,"categories/I");
 
 
 
 
 
 
 for (Int_t iCat = 0; iCat < 2; iCat++){
   TH1D* hMass = (TH1D*) file0.Get(Form("QCD_cat%d;1",iCat));
   

   TAxis* Axis =   hMass->GetXaxis();
   for (int i = 1 ; i < hMass->GetNbinsX()+1; i++){
     double N = fabs(hMass->GetBinContent(i));
      
     mgg = Axis->GetBinCenter(i);
     
     normWeight = N;
     categories = iCat;
 
     if (i%10 == 0) cout << "i = " << i << " N = " << N << " binCenter = " << hMass->GetBinCenter(i)<< " categories = " << categories << endl;

    if (N > 1e-10 && mgg > 1000 && mgg < 5000) TCVARS->Fill();
   }
 }

 

 TCVARS->Write();
 f1.Close();
 file0.Close();
 

}


