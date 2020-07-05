#!/bin/bash

inputFile="file_list.txt"

while IFS= read -r line
do
    root -l -q 'weight.C($inputFile)'
    echo "$line"

done <"$inputFile"