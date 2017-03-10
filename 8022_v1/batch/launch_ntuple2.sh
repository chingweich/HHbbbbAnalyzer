#!/bin/bash

scriptname=`basename $0`
EXPECTED_ARGS=3
if [ $# -ne $EXPECTED_ARGS ]
then
    echo "Usage: $scriptname inputTextFile"
    echo "Example: ./$scriptname JetMET_Run2016B"
    exit 1
fi

export PRODHOME=`pwd`
inputfile=$1

#iteration=2634
iteration=$2
#iteration=209
lastfile=$3
#lastfile=`cat $inputfile | wc -l`
echo "There are "$lastfile" root files"

while [ $iteration -lt $lastfile ]; 
do
  iteration=$(( iteration + 1 ))
  input=(`head -n $iteration $inputfile  | tail -1`)
  echo $iteration $input
bsub -q1nd -R "rusage[mem=3000]" $PWD/runNtuple.sh $PWD $iteration $input
done




