#!/bin/bash
for ((k=22; k<=34; k=k+1 ))
#for ((k=27; k<=32; k=k+1 ))  
do
echo $k
./makeAna2.csh $k
done