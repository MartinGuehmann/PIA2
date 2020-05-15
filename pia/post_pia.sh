#!/bin/bash

# Get the directory where this script is
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

### This collection of scripts is written in bash, perl, and python.
###
### Make sure dependencies are installed and configured before running.
### Dependencies include: BioPerl, NumPy, SciPy, statsmodels (for python).
###
### With the exception of cleanhits.py and the customized long_branch_finder2.py, 
### all other scripts can be found as part of the osiris_phylogenetics galaxy module
### available here: https://bitbucket.org/osiris_phylogenetics/osiris_phylogenetics
### 
### Special thanks to Todd Oakley's lab at UCSB for the original galaxy scripts on 
### which we have based this pipeline.
###
###  Jorge L. Perez-Moreno, Danielle DeLeo
###  CRUSTOMICS Lab at Florida International University 04/05/2017



###  Identify tips longer than ___ median absolute deviations of the tree's branch lengths. Modify number to change MAD multiplier.
###  long_branch_finder2.py can be reverted to calculate standard deviations instead of MADs. More info as comments in the script.

python "$DIR"/../osiris_phylogenetics/phylogenies/long_branch_finder.py "${1}.${2}.treeout.csv" 4 > "${1}.${2}.hits_to_prune.list"


### Clean the output of the long branch finder to avoid conflicts downstream.

python "$DIR"/cleanhits.py "${1}.${2}.hits_to_prune.list" > "${1}.${2}.hits_to_prune.clean.list"


### Then Remove the | from old assemblies.

sed -ie "s/|/_/g" "${1}.${2}.allhits.csv"


### Remove entries from PIA results phytab file that match to a list.

python "$DIR"/../osiris_phylogenetics/phyloconversion/prune_phytab_using_list.py "${1}.${2}.allhits.csv" "${1}.${2}.hits_to_prune.clean.list" discard > "${1}.${2}.allhits.pruned.csv"


### Convert back to FASTA


awk '{print ">"$2"_"$3"_"$1"\n"$4}' "${1}.${2}.allhits.pruned.csv" > "${1}.${2}.allhits.pruned.fasta"


### Remove duplicated sequences resulting from translation of similar isoforms.
if [ -s "${1}.${2}.allhits.pruned.fasta" ]
then
	usearch -cluster_fast "${1}.${2}.allhits.pruned.fasta" -sort length -id 1.00 -threads ${3} -centroids "${1}.${2}.PIA.results.fasta"
else
	cp "${1}.${2}.allhits.fasta" "${1}.${2}.PIA.results.fasta"
fi
