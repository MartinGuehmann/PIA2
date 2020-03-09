#!/bin/bash

# Get the directory where this script is
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

#$1 infile
#$2 outfile
#$3 tree type (ie phylogram)
#$4 yes|no exclude tips
#$5 yes|no label taxa
#$6 name of Rfile "Rfile" by default in xml
#$7 yes|no to label OTUs with QUERY in title
#$8 yes|no to conduct midpoint rooting

#First call perl script which reads trees and writes 
"$DIR"/makeRtrees.pl $1 $2 $3 $4 $5 $7 $8 > $6 2>log.txt

R --vanilla < $6 2>log.txt

