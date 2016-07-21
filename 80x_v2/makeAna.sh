#!/bin/bash
for ((k=22; k<=30; k=k+1 ))
#for ((k=27; k<=32; k=k+1 ))  
do
echo $k
./makeAna.csh $k
done
