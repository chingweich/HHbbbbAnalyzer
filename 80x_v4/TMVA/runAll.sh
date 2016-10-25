#!/bin/bash
for ((k=5; k>=1; k=k-1 ))
do
    for ((j=1; j<=11; j=j+1 ))
    do
	width=$((k*5+15+90+j*5))
	if [ $width -gt 166 ]; then
	    continue;
	fi
	
	bmin=$((j*5+90))
	if [ $bmin -gt 116 ]; then
	    continue;
	fi
	echo MassPlotFineBins_subtr_Moriond_Silver"$bmin"to"$width".root
	#mv MassPlotFineBins_subtr_Moriond_Silver"$bmin"to"$width".root MassPlotFineBins_subtr_Moriond_Silver.root
	#mv MassPlotFineBins_subtr_Moriond_Silver.root MassPlotFineBins_subtr_Moriond_Silver"$bmin"to"$width".root
    done
done


