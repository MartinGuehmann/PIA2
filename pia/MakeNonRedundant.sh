#!/bin/bash

inputFile=$1

filePath="$(dirname $inputFile)"
fileBase="$(basename $inputFile)"

extension="${fileBase##*.}"
fileBase="${fileBase%.*}"

outputFile="$filePath/$fileBase.nr90.$extension"
numTreads=$(nproc)

# This should work but it does not
#cd-hit -i <(seqkit sort -j $numTreads $inputFile) -o $outputFile -c 0.9 -M 0 -d 0 -T $numTreads >&2

tmpFile="$filePath/$fileBase.tmp.$extension"
seqkit sort -j $numTreads -i $inputFile > $tmpFile

cd-hit -i $tmpFile -o $outputFile -c 0.9 -M 0 -d 0 -T $numTreads >&2

rm -f $tmpFile

# Return the outfile
echo $outputFile
