#!/bin/bash

rm -rf opt_*
for i in 8 9 10 11 12 13 14 15;
do
   
		
		mkdir 'opt_'$i
		root -l -b -q ../HHbbbbDrawerNum.C++\($i\)
		#	for p in 1000 1200 1400 1600 2000 2500 3000 3500 4000 4500;
		#	do
		SCRIPTPATH='/home/nstar/HHbbbbAnalyzer/1129DEta'
		
		cd /home/nstar/RKGlobalAnalyzer/Limits/LimitScript
		
		python MakeDataCards.py $SCRIPTPATH'/opt_'$i'/cout.txt' $SCRIPTPATH'/opt_'$i'/'$p'/'
		cd $SCRIPTPATH
		
		#	cp 'opt_'$i'/'$p'/DataCard_M'$p'GeV_MonoHbb_13TeV.txt' 'opt_'$i'/'
		#	cd 'opt_'$i'/'
		#	rm -rf 'cout_'$p'.txt'
		#	rm -rf $p'/'
		#	cd ..
		#	done

done


