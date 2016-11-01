#!/bin/bash
#for ((k=300; k<=; k=k+1 ))
#for ((k=27; k<=32; k=k+1 ))  
#do
#echo $k
#./makeAna2.csh $k
#done
a=(0 300 400 500 600 700 800 900 1000 1250)
for ((k=0;k<=10;k=k+1))
do
mkdir ${a[$k]}
cd ${a[$k]}
#rm -rf *root 
cd /afs/cern.ch/work/c/chchen/HHbbbbAnalyzer/80x_v4/TheaHCorr/pt2  
done
pwd
./makeAna2.csh 0 5000
./makeAna2.csh 300 400
./makeAna2.csh 400 500
./makeAna2.csh 500 600
./makeAna2.csh 600 700
./makeAna2.csh 700 800
./makeAna2.csh 800 900
./makeAna2.csh 900 1000
./makeAna2.csh 1000 1250
./makeAna2.csh 1250 1500
for ((k=0;k<=10;k=k+1))
do
cd ${a[$k]}
rm -rf *merge.root
hadd Bmerge.root *root
cd /afs/cern.ch/work/c/chchen/HHbbbbAnalyzer/80x_v4/TheaHCorr/pt2
done