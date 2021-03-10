#!/bin/bash

inputFile=$1

filePath="$(dirname $inputFile)"
fileBase="$(basename $inputFile)"

extension="${fileBase##*.}"
fileBase="${fileBase%.*}"

outputFile="$filePath/$fileBase.nr90.$extension"
numTreads=$(nproc)

cd-hit -i $inputFile -o $outputFile -c 0.9 -M 0 -d 0 -T $numTreads

# Return the outfile
echo $outputFile
