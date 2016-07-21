wget https://raw.githubusercontent.com/chingweich/HHbbbbAnalyzer/master/80x_v2/makeAna.sh
wget https://raw.githubusercontent.com/chingweich/HHbbbbAnalyzer/master/80x_v2/makeAna.csh
wget https://raw.githubusercontent.com/chingweich/HHbbbbAnalyzer/master/80x_v2/makeAna2.sh
wget https://raw.githubusercontent.com/chingweich/HHbbbbAnalyzer/master/80x_v2/makeAna2.csh
wget https://raw.githubusercontent.com/chingweich/HHbbbbAnalyzer/master/80x_v2/dataMCplots.C
wget https://raw.githubusercontent.com/chingweich/HHbbbbAnalyzer/master/80x_v2/HH4bBtagEffBase_80.C
wget https://raw.githubusercontent.com/chingweich/HHbbbbAnalyzer/master/80x_v2/HH4bBtagMakeEffBase_80.C
wget https://raw.githubusercontent.com/chingweich/HHbbbbAnalyzer/master/80x_v2/HH4bBtagMakeEff_80.C
wget https://raw.githubusercontent.com/chingweich/HHbbbbAnalyzer/master/80x_v2/HHbbbbBtagEff.C
wget https://raw.githubusercontent.com/chingweich/HHbbbbAnalyzer/master/80x_v2/standalone_LumiReWeighting.cc
wget https://raw.githubusercontent.com/chingweich/HHbbbbAnalyzer/master/80x_v2/standalone_LumiReWeighting.h
wget https://raw.githubusercontent.com/chingweich/HHbbbbAnalyzer/master/80x_v2/untuplizer.h
wget https://raw.githubusercontent.com/chingweich/HHbbbbAnalyzer/master/80x_v2/setNCUStyle.C

wget https://twiki.cern.ch/twiki/pub/CMS/BtagRecommendation80X/CSVv2_ichep.csv
mv CSVv2_ichep.csv CSVv2_4invfb.csv
wget https://raw.githubusercontent.com/cms-sw/cmssw/CMSSW_8_0_X/CondTools/BTau/test/BTagCalibrationStandalone.cpp
wget https://raw.githubusercontent.com/cms-sw/cmssw/CMSSW_8_0_X/CondTools/BTau/test/BTagCalibrationStandalone.h

mkdir btagEffSource
mkdir sf
bash
source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.06.06/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh

chmod u+x *sh

./makeAna.sh
cd btagEffSource
mv data1.root data.root
cd ..
./makeAna2.sh

mkdir dataMC
cd dataMC
mkdir 0b
mkdir 1b
mkdir 2b
mkdir all

cd ..
root -l dataMCplots.C++