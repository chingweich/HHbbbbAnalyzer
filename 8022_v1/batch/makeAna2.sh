#!/bin/bash
a=(2636 870 1459 1240 905 2135 24)

for ((k=5; k<=6; k=k+1 ))
#for ((k=27; k<=32; k=k+1 ))  
do
echo $k
for ((j=1; j<=${a[$k]} ; j=j+1 ))
do
./makeAna2.csh $k $j
done
#./makeAna2.csh $k
done

