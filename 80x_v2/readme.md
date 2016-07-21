wget https://raw.githubusercontent.com/chingweich/HHbbbbAnalyzer/master/80x_v2/makeAna.sh<br>
wget https://raw.githubusercontent.com/chingweich/HHbbbbAnalyzer/master/80x_v2/makeAna.csh<br>
wget https://raw.githubusercontent.com/chingweich/HHbbbbAnalyzer/master/80x_v2/makeAna2.sh<br>
wget https://raw.githubusercontent.com/chingweich/HHbbbbAnalyzer/master/80x_v2/makeAna2.csh<br>
wget https://raw.githubusercontent.com/chingweich/HHbbbbAnalyzer/master/80x_v2/dataMCplots.C<br>
wget https://raw.githubusercontent.com/chingweich/HHbbbbAnalyzer/master/80x_v2/HH4bBtagEffBase_80.C<br>
wget https://raw.githubusercontent.com/chingweich/HHbbbbAnalyzer/master/80x_v2/HH4bBtagMakeEffBase_80.C<br>
wget https://raw.githubusercontent.com/chingweich/HHbbbbAnalyzer/master/80x_v2/HH4bBtagMakeEff_80.C<br>
wget https://raw.githubusercontent.com/chingweich/HHbbbbAnalyzer/master/80x_v2/HHbbbbBtagEff.C<br>
wget https://raw.githubusercontent.com/chingweich/HHbbbbAnalyzer/master/80x_v2/standalone_LumiReWeighting.cc<br>
wget https://raw.githubusercontent.com/chingweich/HHbbbbAnalyzer/master/80x_v2/standalone_LumiReWeighting.h<br>
wget https://raw.githubusercontent.com/chingweich/HHbbbbAnalyzer/master/80x_v2/untuplizer.h<br>
wget https://raw.githubusercontent.com/chingweich/HHbbbbAnalyzer/master/80x_v2/setNCUStyle.C<br>

wget https://twiki.cern.ch/twiki/pub/CMS/BtagRecommendation80X/CSVv2_ichep.csv<br>
mv CSVv2_ichep.csv CSVv2_4invfb.csv<br>
wget https://raw.githubusercontent.com/cms-sw/cmssw/CMSSW_8_0_X/CondTools/BTau/test/BTagCalibrationStandalone.cpp<br>
wget https://raw.githubusercontent.com/cms-sw/cmssw/CMSSW_8_0_X/CondTools/BTau/test/BTagCalibrationStandalone.h<br>

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