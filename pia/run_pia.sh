#!/bin/bash

# Get the directory where this script is
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
CURDIR="$(pwd)"
SUBDIR="pia"

# Idiomatic parameter and option handling in sh
# Adapted from https://superuser.com/questions/186272/check-if-any-of-the-parameters-to-a-bash-script-match-a-string
# And advanced version is here https://stackoverflow.com/questions/7069682/how-to-get-arguments-with-flags-in-bash/7069755#7069755
while test $# -gt 0
do
    case "$1" in
        --rebuildDatabases)
            echo "Existing BLAST databases will be rebuilt"
            rebuildDatabases="true" # Could be any value, we actually check, whether it has been defined
            ;;
        --*)
            echo "bad option $1 is ignored"
            ;;
        *)
            echo "bad option $1 is ignored"
            ;;
    esac
    shift
done

aalength="30"			#minimum aminoacid sequence length
search_type="single"	#single or set
gene="r_opsin"			#gene(set) name
evalue="0.00000000000000000001"	#E-value threshold for BLAST search
blasthits="100"			#Number of BLAST hits to retain for the analysis
numThreads=$(nproc)     # Get the number of the currently available processing units to this process, which maybe less than the number of online processors

for file in *.fasta ; do

	if [ ! -d "$SUBDIR" ]
	then
		mkdir "$SUBDIR"
	fi

	cd "$SUBDIR"

	RESULTS_FILE="results_${file}"

	if [ ! -d "$RESULTS_FILE" ]
	then
		mkdir "$RESULTS_FILE"
	fi

	cd "$RESULTS_FILE"

	filebase=$(basename ${file} ".fasta")
	ORF_FILE_BASE="${filebase}_ORF_${aalength}aa"
	ORF_FILE="${ORF_FILE_BASE}.fasta"

 	if [ ! -f "$ORF_FILE" -o "$rebuildDatabases" ] # Build the dataBase if it does not exits or the rebuild argument is supplied
	then
		# The file get_orfs_or_cdss.py is from the pico_galaxy repository. Install that next to your PIA2/pia/ directory. So this is two levels up of this script.
		python "$DIR"/../../pico_galaxy/tools/get_orfs_or_cdss/get_orfs_or_cdss.py -i <(tr -d '\000' < $CURDIR/$file) -e open -m all --min_len $aalength --op "$ORF_FILE"
	fi

	perl "$DIR"/pia.pl "$ORF_FILE" $search_type $gene mafft $evalue $blasthits $numThreads

	perl "$DIR"/phylographics/makeRtrees.pl treeout.tab trees.pdf phylogram no None Rfile yes no >tree.R

	R --vanilla < tree.R 2>log.txt

	"$DIR"/post_pia.sh ${ORF_FILE_BASE} ${gene} ${numThreads}
	cd ../../

done



