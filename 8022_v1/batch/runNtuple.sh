#!/bin/sh

echo $1
cd $1
export SCRAM_ARCH=slc6_amd64_gcc530; eval `scramv1 runtime -sh`
export X509_USER_PROXY=$HOME/private/grid.proxy
root -q -b run.C\($2,\"$3\"\)
rm -rf core*