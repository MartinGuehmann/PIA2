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

rebuildDatabases="0"    # Whether the BLAST database should be rebuild
rebuildTrees="0"        # Whether to rebuild trees
aalength="30"           # Minimum aminoacid sequence length
search_type="single"    # Single or set
gene="r_opsin"          # Gene(set) name
evalue="0.00000000000000000001"	#E-value threshold for BLAST search
blasthits="100"         # Number of BLAST hits to retain for the analysis
numThreads=$(nproc)     # Get the number of the currently available processing units to this process, which maybe less than the number of online processors

# Idiomatic parameter and option handling in sh
# Adapted from https://superuser.com/questions/186272/check-if-any-of-the-parameters-to-a-bash-script-match-a-string
# And advanced version is here https://stackoverflow.com/questions/7069682/how-to-get-arguments-with-flags-in-bash/7069755#7069755
while test $# -gt 0
do
    case "$1" in
        --aalength)
            shift
            aalength="$1"
            ;;
        --search_type)
            shift
            search_type="$1"
            ;;
        --gene)
            shift
            gene="$1"
            ;;
        --evalue)
            shift
            evalue="$1"
            ;;
        --blasthits)
            shift
            blasthits="$1"
            ;;
        --getORFsagain)
            echo "Get open reading frames again, rebuild the BLAST database(s), and the RAxML tree(s)"
            getORFsAgain="true"
            rebuildDatabases="1"
            rebuildTrees="1"
            ;;
        --rebuildDatabases)
            echo "Existing BLAST database(s) and RAxML tree(s) will be rebuilt"
            rebuildDatabases="true" # Could be any value, we actually check, whether it has been defined
            rebuildTrees="1"
            ;;
        --rebuildTrees)
            echo "Rebuild RAxML tree(s)"
            rebuildTrees="1"
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

file="All.fasta"

filebase=$(basename ${file} ".fasta")
ORF_FILE_BASE="${filebase}_ORF_${aalength}aa"
ORF_FILE="${ORF_FILE_BASE}.fasta"

perl "$DIR"/pia.pl "$ORF_FILE" $search_type $gene mafft $evalue $blasthits $numThreads $rebuildDatabases $rebuildTrees

perl "$DIR"/phylographics/makeRtrees.pl "${ORF_FILE_BASE}.$gene.treeout.csv" "${ORF_FILE_BASE}.$gene.trees.pdf" phylogram no None Rfile yes no >"${ORF_FILE_BASE}.$gene.tree.R"

R --vanilla --slave < "${ORF_FILE_BASE}.$gene.tree.R" 2>R-stderr-Log.txt

"$DIR"/post_pia.sh ${ORF_FILE_BASE} ${gene} ${numThreads}
