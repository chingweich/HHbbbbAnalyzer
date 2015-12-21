{


  double mgg, mjj,evWeight, mtot, normWeight;
 int categories;

 evWeight = 1.0;
 normWeight = 1;

   
 // string inFile("input/MassPlotFineBins.root");
 string inFile("input/MassPlotFineBins_2.root");        
 string outFile("dijetHH_RadionHH");
 
 int massrange=36;

 TFile file0(inFile.c_str(), "read");
 
 for (int iMass = 0; iMass<massrange; iMass++){
   
   // if ( iMass == 1 || iMass == 3  || iMass == 5  || iMass == 7 || iMass == 9 ) continue;
   if ( iMass == 1 || iMass == 3  || iMass == 5  || iMass == 7 ||iMass == 9 || iMass == 8)continue;
   if (iMass > 10 && iMass < 15) continue;
   if (iMass > 15 && iMass < 20) continue;
   if (iMass > 20 && iMass < 25) continue;

   if (iMass > 25 && iMass < 30) continue;
   if (iMass > 30 && iMass < 35) continue;
   
   int Mass = 1000+iMass*100;   
   
   string sOutFile = "MiniTrees/Signal_HH_13TeV/" + outFile + Form("OUT%d_miniTree.root", Mass);
   TFile* f1 = new TFile(sOutFile.c_str(), "recreate");
   
   TTree *TCVARS = new TTree("TCVARS", "hh selection");
   TCVARS->Branch("mgg",&mgg,"mgg/D");
   
   TCVARS->Branch("evWeight",&evWeight,"evWeight/D");
   TCVARS->Branch("normWeight",&normWeight,"normWeight/D");
   
   TCVARS->Branch("categories",&categories,"categories/I");
   
   

   
   
   for (int iCat = 0; iCat < 2; iCat++){
     TH1D* hMass = (TH1D*) file0.Get(Form("Radion%d_cat%d;1",Mass,iCat));
     
     
     TAxis* Axis =   hMass->GetXaxis();
     for (int i = 1 ; i < hMass->GetNbinsX()+1; i++){
       int N = abs(hMass->GetBinContent(i));
       
       if (i%10 == 0) cout << "i = " << i << "N = " << N << " binCenter = " << hMass->GetBinCenter(i) << endl;
       
       mgg = Axis->GetBinCenter(i);
       
       normWeight = N;
       categories = iCat;
       if (N > 1e-10 && mgg > 999 && mgg < 5000) TCVARS->Fill();
     }
   }
   
   TCVARS->Write();
   f1->Close();
   
     
 }

 file0.Close();

}


