#!/bin/tcsh

#for ((i=9;i<18;i++))
#do
#echo $i
#done

#root -q -b -l HHbbbbAnalyzer.C++\($1\)
#root -q -b -l HHbbbbTrigger.C++\($1\)
#root -q -b -l HHbbbbBtagEff.C++\($1\) 
root -q -b -l HH4bBtagMakeEff_80.C++\($1\)
#root -q -b -l HHbbbbBtagEff_76_my.C++\($1\)   
