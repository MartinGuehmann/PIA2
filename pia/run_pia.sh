#!/bin/bash

# Get the directory where this script is
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

aalength="30"			#minimum aminoacid sequence length
search_type="single"	#single or set
gene="r_opsin"			#gene(set) name
evalue="0.00000000000000000001"	#E-value threshold for BLAST search
blasthits="100"			#Number of BLAST hits to retain for the analysis

for file in *.fasta ; do

	mkdir pia
	cd pia
	mkdir results_${file}
	cd results_${file}

	# The file get_orfs_or_cdss.py is from the pico_galaxy repository. Install that next to your PIA2/pia/ directory. So this is two levels up of this script.
	python "$DIR"/../../pico_galaxy/tools/get_orfs_or_cdss/get_orfs_or_cdss.py -i $file -e open -m top --min_len $aalength --op ORF_prot.fasta

	perl "$DIR"/pia.pl ORF_prot.fasta $search_type $gene mafft $evalue $blasthits

	perl "$DIR"/phylographics/makeRtrees.pl treeout.tab trees.pdf phylogram no None Rfile yes no >tree.R

	R --vanilla < tree.R 2>log.txt

	"$DIR"/post_pia.sh
	cd ../../

done



