#!/bin/bash

#
# Converts an amino acid sequence file to a DNA sequence file
# Argument 1 (string): Amino acid sequence file to convert.
# Argument 2 (string): DNA sequence file from that the amino acid sequences where translated
# Reduces the number of sequences if the amino acid sequences are derived from the same DNA sequence
#

inputFile=$1
dataBaseFile=$2

filePath="$(dirname $inputFile)"
fileBase="$(basename $inputFile)"

extension="${fileBase##*.}"
fileBase="${fileBase%.*}"

outputFile="$filePath/$fileBase.nt.$extension"
numTreads=$(nproc)

seqIDs=$(seqkit seq -j $numTreads -n -i $inputFile | sed 's/|ORF.*$//g' | sort | uniq)

genePrefix="${fileBase%%.*}"

rm -f $outputFile

for id in $seqIDs
do
	completePrefix=$(echo $id | grep -o "^.*${genePrefix}_")
	pattern=$(echo $id | sed "s/^.*${genePrefix}_//")
	seqkit grep -j $numTreads -p $pattern $dataBaseFile | sed "s/>/>$completePrefix/" >> $outputFile
done

# Return the outfile
echo $outputFile
