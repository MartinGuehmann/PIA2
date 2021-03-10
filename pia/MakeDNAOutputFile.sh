#!/bin/bash

inputFile=$1
dataBaseFile=$2

filePath="$(dirname $inputFile)"
fileBase="$(basename $inputFile)"

extension="${fileBase##*.}"
fileBase="${fileBase%.*}"

outputFile="$filePath/$fileBase.nt.$extension"
numTreads=$(nproc)

seqIDs=$(seqkit seq -j $numTreads -n -i $inputFile | sed 's/|.*$//g' | sort | uniq)

genePrefix="${fileBase%%.*}"

rm -f $outputFile

for id in $seqIDs
do
	completePrefix=$(echo $id | grep -o "^.*${genePrefix}_")
	pattern=$(echo $id | sed "s/^.*${genePrefix}_//")
	seqkit grep -j $numTreads -p $pattern $dataBaseFile | sed "s/>/>$completePrefix/" >> $outputFile
done
