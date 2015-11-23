#!/bin/bash

rm -rf opt_*
for i in 20 25 30 35 40 45 50 55 60;
do
    for j in 90 95 100 105 110;
    do
	k=$((i+j))
	if test $k -gt 150;then
	        continue;
	fi
	if test $k -le 125;then
            continue;
        fi
	#for l in 20 25 30 35 40 45 50 55 60;
	for l in 30 ;
	do
	    #    for m in 90 95 100 105 110 ;
	    for m in 105 ;
	    do
		n=$((i+j))
		o=$((l+m))
		if test $o -gt 150;then
		    continue;
		fi
		if test $o -le 125;then
                    continue;
		fi
		
		mkdir 'opt_'$j'_'$n'_'$m'_'$o
		root -l -b -q ../HHbbbbDrawerNum.C++\($j,$n,$m,$o\)
			for p in 1000 1200 1400 1600 2000 2500 3000 3500 4000 4500;
			do
			SCRIPTPATH='/home/nstar/HHbbbbAnalyzer/1116Op'
			
			cd /home/nstar/RKGlobalAnalyzer/Limits/LimitScript
			#echo $SCRIPTPATH
			python MakeDataCards.py $SCRIPTPATH'/opt_'$j'_'$n'_'$m'_'$o'/cout_'$p'.txt' $SCRIPTPATH'/opt_'$j'_'$n'_'$m'_'$o'/'$p'/'
			cd $SCRIPTPATH
			#pwd
			cp 'opt_'$j'_'$n'_'$m'_'$o'/'$p'/DataCard_M'$p'GeV_MonoHbb_13TeV.txt' 'opt_'$j'_'$n'_'$m'_'$o'/'
			cd 'opt_'$j'_'$n'_'$m'_'$o'/'
			rm -rf 'cout_'$p'.txt'
			rm -rf $p'/'
			cd ..
			done
	    done
	done
    done
done


