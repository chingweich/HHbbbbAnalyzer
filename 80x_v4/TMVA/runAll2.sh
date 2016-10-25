#!/bin/bash
file='mass.txt'
exec < $file 

while read line
do
    echo $line
    
done
