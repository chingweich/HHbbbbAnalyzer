
void run(int w, std::string fileName)
{
  TStopwatch* myWatch = new TStopwatch();
  myWatch->Start();
  gROOT->ProcessLine(".L skimTree_C.so");
  //  gROOT->ProcessLine(".L skimTree.C++");
  gROOT->ProcessLine(Form("skimTree(%d,\"%s\")",w,fileName.data()));

  myWatch->Stop();
  cout << myWatch->RealTime() << " seconds has passed! " << endl; 

}
