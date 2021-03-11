#!/bin/bash

inputFile=$1

filePath="$(dirname $inputFile)"
fileBase="$(basename $inputFile)"

extension="${fileBase##*.}"
fileBase="${fileBase%.*}"

outputFile="$filePath/$fileBase.outgroupFree.$extension"
numTreads=$(nproc)

seqkit grep -j $numTreads -v -r -p "Outgroup_" $inputFile > $outputFile

# Return the outfile
echo $outputFile
