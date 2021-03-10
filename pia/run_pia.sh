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
RESULTS_PREFIX="results_"
ALL_DEST="All"
RESULTS_FILE_ALL="${RESULTS_PREFIX}${ALL_DEST}"

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

if [ ! -d "$SUBDIR" ]
then
	mkdir "$SUBDIR"
fi

if [ ! -d "$SUBDIR/$RESULTS_FILE_ALL" ]
then
	mkdir "$SUBDIR/$RESULTS_FILE_ALL"
fi

rm -f "$CURDIR/$SUBDIR/$RESULTS_FILE_ALL/${ALL_DEST}_ORF_${aalength}aa.eValue=$evalue.$gene.allhits.fasta"

for file in *.fasta ; do

	cd "$SUBDIR"

	RESULTS_FILE="${RESULTS_PREFIX}${file}"

	if [ ! -d "$RESULTS_FILE" ]
	then
		mkdir "$RESULTS_FILE"
	fi

	cd "$RESULTS_FILE"

	filebase=$(basename ${file} ".fasta")
	ORF_FILE_BASE="${filebase}_ORF_${aalength}aa"
	ORF_FILE="${ORF_FILE_BASE}.fasta"
	FINAL_FILE_BASE="${ORF_FILE_BASE}.eValue=$evalue.$gene"

	if [ ! -f "$ORF_FILE" -o "$getORFsAgain" ] # Build the dataBase if it does not exits or the rebuild argument is supplied
	then
		echo "Extract open reading frames from $file"
		# The file get_orfs_or_cdss.py is from the pico_galaxy repository. Install that next to your PIA2/pia/ directory. So this is two levels up of this script.
		python3 "$DIR"/../../pico_galaxy/tools/get_orfs_or_cdss/get_orfs_or_cdss.py -i <(tr -d '\000' < $CURDIR/$file) -e open -m all --min_len $aalength --op "$ORF_FILE"
	fi

	perl "$DIR"/pia.pl "$ORF_FILE" $search_type $gene mafft $evalue $blasthits $numThreads $rebuildDatabases $rebuildTrees

	perl "$DIR"/phylographics/makeRtrees.pl "${FINAL_FILE_BASE}.treeout.csv" "${FINAL_FILE_BASE}.trees.pdf" phylogram no None Rfile yes no >"${FINAL_FILE_BASE}.tree.R"

	R --vanilla --slave < "${FINAL_FILE_BASE}.tree.R" 2>R-stderr-Log.txt

	# There is some problem with this, but actually we do not need this.
#	"$DIR"/post_pia.sh ${FINAL_FILE_BASE} ${numThreads}
	cd ../../

	cat "$CURDIR/$SUBDIR/$RESULTS_FILE/${FINAL_FILE_BASE}.allhits.fasta" >> "$CURDIR/$SUBDIR/$RESULTS_FILE_ALL/${ALL_DEST}_ORF_${aalength}aa.eValue=$evalue.$gene.allhits.fasta"

	"$DIR/MakeNonRedundant.sh" "$CURDIR/$SUBDIR/$RESULTS_FILE/${FINAL_FILE_BASE}.allhits.fasta"
done

"$DIR/MakeNonRedundant.sh" "$CURDIR/$SUBDIR/$RESULTS_FILE_ALL/${ALL_DEST}_ORF_${aalength}aa.eValue=$evalue.$gene.allhits.fasta"

cd "$SUBDIR/$RESULTS_FILE_ALL"

file="${ALL_DEST}.fasta"

filebase=$(basename ${file} ".fasta")
ORF_FILE_BASE="${filebase}_ORF_${aalength}aa"
ORF_FILE="${ORF_FILE_BASE}.fasta"
FINAL_FILE_BASE="${ORF_FILE_BASE}.eValue=$evalue.$gene"

perl "$DIR"/pia.pl "$ORF_FILE" $search_type $gene mafft $evalue $blasthits $numThreads $rebuildDatabases $rebuildTrees

perl "$DIR"/phylographics/makeRtrees.pl "${FINAL_FILE_BASE}.treeout.csv" "${FINAL_FILE_BASE}.trees.pdf" phylogram no None Rfile yes no >"${FINAL_FILE_BASE}.tree.R"

R --vanilla --slave < "${FINAL_FILE_BASE}.tree.R" 2>R-stderr-Log.txt

# *.allhits.cvs is missing, but doing this is not necessary for now
#"$DIR"/post_pia.sh ${FINAL_FILE_BASE} ${numThreads}
