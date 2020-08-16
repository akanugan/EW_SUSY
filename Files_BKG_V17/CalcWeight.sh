#!/bin/bash

inputFile="file_list.txt"

while IFS= read -r line
do
    echo "$line" 
    
#    root -l -b -q 'weight.C("$line")'
    echo "************"
done <"$inputFile"
