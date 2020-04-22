#!/usr/bin/perl -w
use strict;
use Bio::Perl;
use Bio::DB::GenBank;
use Bio::SeqIO;
use Bio::Seq;
use Getopt::Long;
use IO::File;
use FindBin;
use lib $FindBin::Bin;
use FastaDb;
use FastqDb;
use Path::Tiny;
use Scalar::Util qw(openhandle);

# Get the directory where PIA is. This is the directory of this script:
use File::Basename ();
my $PIADIR       = File::Basename::dirname($0);
my $SUBDIR       = "LIT_1.1";
my $dataFileName = "LIT_1.1.csv";

# The directory where precalculated genetrees and LIT_1.1.csv file are located must match this:
my $DATADIR      = path($0)->parent->child("$SUBDIR");
my $dataFile     = path("$DATADIR")->child($dataFileName);

# First, read in the data table with gene name, path, etc
# Place these data into a hash
open(BLASTFILE, "<$dataFile") or die "File $dataFileName file must be available check full path in pia.pl";

my @genedata;
my %HoH;

while(<BLASTFILE>) {
	my $currentinput = "$_";
	if($currentinput =~ /\t/){
		my @genedata = split(/\t/);
		my $genename = $genedata[0];
		$HoH{$genename}{set}         = $genedata[1];
		$HoH{$genename}{bait}        = $genedata[2];
		$HoH{$genename}{fastatag}    = $genedata[3];
		$HoH{$genename}{reftreename} = $genedata[4];
		$HoH{$genename}{outgroup}    = $genedata[5];
	}else{
		die "File named genelist.txt must be availalbe in tab delimited format with gene data\n";
	}
}

close(BLASTFILE);

#To print hash for debug purposes
#for my $loop (keys %HoH) {
#	print "$loop: ";
 #   	for my $role ( keys %{ $HoH{$loop} } ) {
  #  	     print "$role=$HoH{$loop}{$role} ";
   # 	}
    #	print "\n";
#}
#exit;

#Get arguments from xml in galaxy
#my $command       = shift(@ARGV);      # 0 pia.pl command
my $data_file      = shift(@ARGV);      # 1 data file to search
my $scope          = shift(@ARGV);      # 2 scope [all|single|set]
my $genefamily     = shift(@ARGV);      # 3 [NULL|name of set|name of family|
my $alignment      = shift(@ARGV);      # 4 Alignment algorithm, e.g. mafft
my $evalue         = shift(@ARGV);      # 5 E-value
my $maxkeep        = shift(@ARGV);      # 6 Maximum blast hits to retain
my $numThreads     = shift(@ARGV);      # 7 The number of CPU cores
my $rebuildBLASTdb = shift(@ARGV);      # 8 Whether to rebuild the BLAST database
my $rebuilTrees    = shift(@ARGV);      # 9 Whether to rebuild the RAxML trees

my @genes2analyze;

#if $scope = functional then loop through all genes and analyze those that match to the set
if($scope eq "single"){
	push(@genes2analyze, $genefamily);
}else{ #scope = functional
	# Add all genes from set into an array
	print "Analyzing functional gene class\n";
	for my $hashgene (keys %HoH) {
		if($HoH{$hashgene}{set} eq $genefamily) {
			push(@genes2analyze, $hashgene);
			print "\t".$hashgene."\n";
		}
	}
}

my($filename, $dirs, $suffix) = File::Basename::fileparse($data_file, qr/\.[^.]*/);
my $BLASTDB = $filename . ".blastDB";

# Files to retain all hits found from all gene families in the run
my $AllHitsCSVFile = $filename . "." . $genefamily . ".allhits.csv";
my $AllHitsFasFile = $filename . "." . $genefamily . ".allhits.fasta";
my $FinalTreeFile  = $filename . "." . $genefamily . ".treeout.csv";

if(-f $AllHitsCSVFile)
{
	unlink $AllHitsCSVFile;
}

if(-f $AllHitsFasFile)
{
	unlink  $AllHitsFasFile;
}

if(-f $FinalTreeFile)
{
	unlink $FinalTreeFile;
}

if( ! -f "$BLASTDB.phr" or $rebuildBLASTdb)
{
	print ".........Creating Blast Database: $BLASTDB\n";
	system "makeblastdb -in $data_file -out $BLASTDB -dbtype prot";

	$rebuilTrees = 1;
}
else
{
	print ".........Use existing Blast database: $BLASTDB\n";
}

my $thisgene;
while($thisgene = shift(@genes2analyze) ) {
	print "\n\n\n**************Analyzing $thisgene ..............\n";
	# Use get_gb to obtain sequence from genbank using accession(s) that are looked up based on gene name
	# bait written to tempfile.tmp using get_gb subscript

	my $fh;
	my $query = "";
	open($fh, ">", \$query);
	my $accession = $HoH{$thisgene}{bait};

	# No return value, $fh is passed by reference
	get_gb('protein', 'fasta', $fh, $accession, '', '');
	close($fh);

	open($fh, "<", \$query);

	# Next BLASTP data
	system "blastp -version";
	print "...Searching data with BLASTP using an evalue of $evalue and retaining a maximum of $maxkeep (in case of tie, all genes are retained)\n\n";

	my $command = "blastp"
	            . " -query <(echo -e \"$query\") "
	            . " -db "              . $BLASTDB
	            . " -outfmt "          . 6
	            . " -num_threads "     . $numThreads
	            . " -evalue "          . $evalue
	            . " -max_target_seqs " . $maxkeep;

	my $finds   = qx(bash -c '$command');
	close($fh);

	# Add some text to the file when no blast hits have been found
	if($finds eq "")
	{
		$finds = "EMPTY\tEMPTY\tEMPTY\t**No Hits found\n";
	}

	print $finds;

	open($fh, "<", \$finds);
	my $decoratedHits = "";
	my $hits;
	open($hits, ">", \$decoratedHits);

	# Grab sequences from the input database. hits is temporary file containing hits that get_seqs will write to
	get_seqs($data_file, $fh, "2", "", "", $hits, "", "", "", "");
	close($fh);
	# close($hits); # Already closed in sub of get_seqs.

	print $decoratedHits;

	my $fastaToAdd = "";
	open($fh,   ">", \$fastaToAdd);
	open($hits, "<", \$decoratedHits);
	#Next add string to identify gene family to beginning of hits
	addstring2fashead($hits, $HoH{$thisgene}{fastatag}, $fh);
	close($hits);
	close($fh);

	# Place the reads with raxml

	# Add the full path
	my $path = path("$DATADIR")->child($HoH{$thisgene}{set})->child($HoH{$thisgene}{reftreename});
	chomp($path);

	genetree_read_placement($fastaToAdd, $alignment, $path, $thisgene);

	open(my $myHandle, ">>$AllHitsFasFile");
	$myHandle->print($fastaToAdd);
	close($myHandle);
}

sub addstring2fashead
{
	my $infileHandle  = $_[0];
	my $addstring     = $_[1];
	my $addFileHanler = $_[2];
	my $in_obj        = Bio::SeqIO->new(-fh => $infileHandle, '-format' => 'fasta');
	my $currentinput  = $addstring;

	# Grab sequence objects
	my $csvFile;
	open($csvFile, ">>$AllHitsCSVFile");

	while (my $seq = $in_obj->next_seq() ) {
		my $seq_obj = $in_obj;

		my $seq_id = $currentinput.$seq->id; # Combined ID geneName_origID
		print $addFileHanler ">".$seq_id;
		my $seq_seq = $seq->seq;
		print $addFileHanler "\n".$seq_seq."\n";
		# Extract just the name from the 'X_hit' string
		$seq_id =~ s/\|ORF\d+// ;
		my $left = index($seq_id, "_hit_")+5;
		my $fragment = substr $seq_id, $left;
		print $csvFile $fragment;                        # Print originqal ID
		$fragment = substr $seq_id, 0, $left-1;
		print $csvFile "\t".$fragment;                   # Print name of putative gene
		print $csvFile "\t".$seq_seq."\n";               # Print sequence
	}

	close($csvFile);

	return();
}

sub get_gb
{
	my $datatype   = $_[0];
	my $outtype    = $_[1];
	my $fileHandle = $_[2];
	my $manual     = $_[3];
	my $mannames   = $_[4];
	my $genenames  = $_[5];

	my $accessions;
	my @accnums;
	my @newnames;
	my $manbin    = 0;
	my @genenames;
	my $genebin   = 0;

	unless($mannames eq ''){
		@newnames = split(/ /,$mannames);
		$manbin=1;
	}

	unless($genenames eq ''){
		@genenames = split(/ /,$genenames);
		$genebin=1;
	}

	@accnums = split(/,/,$manual);
	my $countnames = 0;

	# A bit ugly, but when $fh goes out of scope it also closes $fileHandle,
	# which is supposed to be closed from the caller, since it opened it.
	my $fh;

	if($outtype ne "phytab")
	{
		$fh = Bio::SeqIO->new(-fh => $fileHandle, -format => $outtype);
	}

	foreach (@accnums){
		#Should check input for one word per line and throw error if not, which is not done

		$accessions = $_;
		chomp;
		if($accessions eq ""){
			die "Put spaces between accession numbers\n";
		}

		my $qry_string .= $accessions."[accession]"." ";
		my $GBseq;
		my $gb = new Bio::DB::GenBank;
		my $query = Bio::DB::Query::GenBank->new
		        (-query   =>$qry_string,
		         -db      =>$datatype);

		my $count;
		my $species;

		if($outtype eq "phytab"){ #print phytab format, do not use bioperl as below.
			if( defined (my $seqio = $gb->get_Stream_by_query($query)) ){

				while( defined ($GBseq = $seqio->next_seq )) {
					my $sequence = $GBseq;   # read a sequence object
					if($manbin ==1){ #replace GenBank Names with Custom Names
						$sequence->id($newnames[$countnames]);
						$sequence->desc('');
						$species = $sequence->id;
						$countnames++;
					}else{
						$species = $sequence->species->binomial;
						$species =~ s/ /_/g ;
					}
					if(@genenames > 0){
						if(@genenames == 1){
							$fileHandle->print($species."\t".$genenames[0]."\t".$sequence->accession."\t".$sequence->seq."\n");
						}else{
							$fileHandle->print($species."\t".$genenames[$countnames-1]."\t".$sequence->accession."\t".$sequence->seq."\n");
						}
					}else{
						$fileHandle->print($species."\tNone\t".$sequence->accession."\t".$sequence->seq."\n");
					}
				}

			}else{
				print "Did not find $accessions\n";
			}
		}else{

			if( defined (my $seqio = $gb->get_Stream_by_query($query)) ){

				while( defined ($GBseq = $seqio->next_seq )) {
					my $sequence = $GBseq;   # read a sequence object
					if($manbin ==1){ # Replace GenBank Names with Custom Names
						$sequence->id($newnames[$countnames]);
						$sequence->desc('');
						$countnames++;
					}
					$fh->write_seq($sequence); # write a sequence object
				}
			}else{
				print "Did not find $accessions\n";
			}
		}
	}
}

sub get_seqs
{
	#
	# This script creates a Fasta/Qual/Fastq file of selected sequences, with optional filters.
	#
	# 02/24/10 : created by Ed Kirton
	# 12/07/10 : fixed Fastq bug
	# 04/18/13 : THO makes this a subroutine of pia.pl 

my $usage = <<'ENDHERE';
NAME:
    get_seqs.pl
PURPOSE:
    To extract a subset of sequences by ID.
INPUT:
    --db <*.fasta|fastq> : file containing sequences in Fasta or Fastq format
    --table <*.tsv> : file containing sequence IDs (optional; default=stdin)
    --col <int> : column of table containing sequence IDs (optional; default=1=first column)
OUTPUT:
    --selected <*.fasta|fastq> : file containing named sequences
    --unselected <*.fasta|fastq> : file containing unselected sequences
OPTIONS:
    --cosorted : uses faster algorithm if IDs appear in both files in the same order
    --paired : filter complete read-pair when one read is selected (requires Illumina-style read IDs; i.e. */1, */2)
    --ignore case : ignore case differences between IDs
    --gzip : compress all outfiles
OPTIONAL FILTERS:
    Optional filters, each in the form of column:condition:value.
    Where column is the column in the table (containing IDs)
    Condition is one of the following:
        String operators:
            s_eq
            s_ne
            s_contains
            s_notcontains
            s_startswith
            s_notstartswith
            s_endswith
            s_notendswith
        Numerical operators:
            n_eq
            n_ne
            n_gt
            n_lt
    Where value is a string or number as appropriate.
AUTHOR/SUPPORT:
    Edward Kirton (ESKirton@LBL.gov)
ENDHERE

	#
	# VALIDATE INPUT
	#
	my ($help, $dbfile, $tablefile, $id_col, $ignorecase, $cosorted, $selected, $unselected, $gzip, $paired);

	$dbfile     = $_[0]; # BLAST database
	$tablefile  = $_[1]; # BLAST output file
	$id_col     = $_[2]; # Get data from that column
	$ignorecase = $_[3];
	$cosorted   = $_[4];
	$selected   = $_[5]; # To be saved
	$unselected = $_[6];
	$gzip       = $_[7];
	$paired     = $_[8];

	#GetOptions(
	#    'd|db=s'           => \$dbfile,
	#    't|table=s'        => \$tablefile,
	#    'c|col=i'          => \$id_col,		
	#    'ignorecase'     => \$ignorecase,
	#    'cosorted'       => \$cosorted,
	#    's|selected=s'   => \$selected,
	#    'u|unselected=s' => \$unselected,
	#    'g|gzip' => \$gzip,
	#    'p|paired' => \$paired,
	#    'h|help'         => \$help
	#);

	if ($help) { print $usage; exit; }

	die("DB required\n")                      unless    $dbfile;
	die("DB file not found: $dbfile\n")       unless -f $dbfile;
	die("Table required\n")                   unless    $tablefile;
	die("Table file not found: $tablefile\n") unless -f $tablefile or $tablefile == openhandle($tablefile);
	$selected   = '' if !defined($selected)   or $selected   eq 'None';
	$unselected = '' if !defined($unselected) or $unselected eq 'None';
	$id_col=1 unless $id_col;
	die("Invalid id column, $id_col\n")       unless $id_col > 0;

	my $filters = [];
	#THO not using filters in subroutine
	#while (my $filter = shift @ARGV) {
	#    next unless $filter;
	#    my @a_filter = split(/:/, $filter);
	#    die("Invalid number of filter options: @a_filter") unless @a_filter == 3;
	#    push @$filters, \@a_filter;
	#}

	#
	# MAIN
	#
	my ($n_selected,$n_unselected);
	if ($cosorted) {
		# SEARCH IS FAST AND EASY IF INPUTS SIMILARLY SORTED!
		($n_selected,$n_unselected) = search_cosorted($dbfile, $tablefile, $id_col, $ignorecase, $selected, $unselected, $paired, $gzip, $filters);
	} else {
		# INPUT NOT CO-SORTED SO KEEP ALL IDS IN RAM
		($n_selected,$n_unselected) = search($dbfile, $tablefile, $id_col, $ignorecase, $selected, $unselected, $paired, $gzip, $filters);
	}

	#print "Selected = $n_selected; Unselected = $n_unselected\n"; 
	return();
}

#
# RETURNS TRUE ONLY IF RECORD MATCHES (OPTIONAL) SEARCH CRITERIA
#
sub match
{
	my ($filters, $row) = @_;
	foreach my $filterA (@$filters) {
		my ($condition, $col, $value) = @$filterA;
		my $x = $row->[ $col - 1 ];
		if    ($condition eq 's_eq')            { return 0 unless $x eq $value }
		elsif ($condition eq 's_ne')            { return 0 unless $x ne $value }
		elsif ($condition eq 's_contains')      { return 0 unless $x =~ /$value/ }
		elsif ($condition eq 's_notcontains')   { return 0 unless $x !~ /$value/ }
		elsif ($condition eq 's_startswith')    { return 0 unless $x =~ /^$value/ }
		elsif ($condition eq 's_notstartswith') { return 0 unless $x !~ /^$value/ }
		elsif ($condition eq 's_endswith')      { return 0 unless $x =~ /$value$/ }
		elsif ($condition eq 's_notendswith')   { return 0 unless $x !~ /$value$/ }
		elsif ($condition eq 'n_eq')            { return 0 unless $x == $value }
		elsif ($condition eq 'n_ne')            { return 0 unless $x != $value }
		elsif ($condition eq 'n_gt')            { return 0 unless $x >  $value }
		elsif ($condition eq 'n_lt')            { return 0 unless $x <  $value }
	}
	return 1;
}

#
# SIMULTANEOUSLY PARSE TWO STREAMS
#
sub search_cosorted
{
	my ($dbfile, $tablefile, $id_col, $ignorecase, $selected, $unselected, $paired, $gzip, $filters) = @_;
	my $sfh          = new IO::File;
	my $ufh          = new IO::File;
	my $table        = new IO::File;
	my $n_selected   = 0;
	my $n_unselected = 0;

	# OPEN FILES
	if ($tablefile) {
		if($tablefile == openhandle($tablefile))
		{
			$table = $tablefile;
		}
		else
		{
			open($table, "<$tablefile") or die("Unable to open file, $tablefile: $!\n");
		}
	} else {
		$table=*STDIN;
	}

	if ($selected) {
		if ($gzip) {
#			open($sfh, '>:gzip', $selected) or die("Unable to open file, $selected: $!\n");
		} else {
			if($selected == openhandle($selected))
			{
				$sfh = $selected;
			}
			else
			{
				open($sfh, ">$selected") or die("Unable to open file, $selected: $!\n");
			}
		}
	} else {
		open($sfh, ">/dev/null");
	}
	if ($unselected) {
		if ($gzip) {
#			open($ufh, '>:gzip', $unselected) or die("Unable to open file, $unselected: $!\n");
		} else {
			if($unselected == openhandle($unselected))
			{
				$ufh = $unselected;
			}
			else
			{
				open($ufh, ">$unselected") or die("Unable to open file, $unselected: $!\n");
			}
		}
	} else {
		open($ufh, ">/dev/null");
	}

	# GET FIRST MATCHING TARGET ID
	my $prev_target_id = '';
	my $target_id = '';
	get_next_matching_target_id($table,$id_col,$ignorecase,$filters,\$target_id,\$prev_target_id,$paired);
	unless ($target_id) {
		# no records match search criteria
		close $table;
		close $sfh if $selected;
		if ($unselected) {
			open(DB, "<$dbfile") or die("Unable to open file, $dbfile: $!\n");
			while (<DB>) {
				print $ufh $_;
				++$n_unselected;
			}
			close DB;
		}
		close $ufh;
		return 0;
	}

	# DETERMINE FILETYPE
	open(DB, "<$dbfile") or die("Unable to open file, $dbfile: $!\n");
	my $format;
	while (<DB>) {
		chomp;
		if (/^#/ or ! $_) { next }
		elsif (/^>/) { $format='fasta' }
		elsif (/^@/) { $format='fastq' }
		else { die "Invalid DB file format" }
		last;
	}
	close DB;

	# PARSE
	my $db = $format eq 'fasta' ? FastaDb->new($dbfile) : FastqDb->new($dbfile);
	while (my $rec=$db->next_seq ) {
		unless ($target_id) {
			last unless $unselected; # done if no more seqs to get
			# otherwise dump rest of seqs in unselected file
			print $ufh $rec->output;
			++$n_unselected;
			while ($rec=$db->next_seq ) {
				print $ufh $rec->output;
				++$n_unselected;
			}
			last;
		}
		my $id=$ignorecase ? uc($rec->id):$rec->id;
		if ($id eq $prev_target_id or $id eq $target_id) {
			# selected seq
			print $sfh $rec->output;
			++$n_selected;
			get_next_matching_target_id($table,$id_col,$ignorecase,$filters,\$target_id,\$prev_target_id,$paired);
		} else {
			# unselected seq
			print $ufh $rec->output;
			++$n_unselected;
		}
	}

	close $table;
	close $sfh;
	close $ufh;

	# If some target seqs not found, it's likely the files were not cosorted, so try unsorted search function.
	if ($target_id) {
		print "Files don't appear to be cosorted, trying unsorted search\n";
		return search($dbfile, $tablefile, $id_col, $ignorecase, $selected, $unselected, $filters);
	}
	return ($n_selected,$n_unselected);

	# Subfunction of search_cosorted
	sub get_next_matching_target_id {
		my ($table,$id_col,$ignorecase,$filters,$target_idR,$prev_target_idR,$paired)=@_;
		$$prev_target_idR = $$target_idR;
		$$target_idR = '';
		while (<$table>) {
			chomp;
			my @row = split(/\t/);
			die("Bad input table") unless @row >= $id_col;
			next unless match($filters, \@row);
			my $new_target_id = $ignorecase ? uc($row[ $id_col - 1 ]) : $row[ $id_col - 1 ];
			$new_target_id=$1 if $new_target_id =~ /^(\S+)/; # use first word only
			$new_target_id=$1 if $paired and $new_target_id =~ /^(\S+)\/[12]$/;
			next if $new_target_id eq $$prev_target_idR;
			$$target_idR=$new_target_id;
			last;    # return to parsing db file
		}
	}
}

#
# LOAD IDS INTO RAM THEN PARSE DB.
#
sub search
{
	my ($dbfile, $tablefile, $id_col, $ignorecase, $selected, $unselected, $paired, $gzip, $filters) = @_;
	my $sfh          = new IO::File;    # selected seqs
	my $ufh          = new IO::File;    # unselected seqs
	my $table        = new IO::File;
	my $n_selected   = 0;
	my $n_unselected = 0;
	my %ids          = ();
	open(DB,    "<$dbfile")    or die("Unable to open file, $dbfile: $!\n");
	if ($tablefile) {
		if($tablefile == openhandle($tablefile))
		{
			$table = $tablefile;
		}
		else
		{
			open($table, "<$tablefile") or die("Unable to open file, $tablefile: $!\n");
		}
	} else {
		$table=*STDIN;
	}
	if ($selected) {
		if ($gzip) {
#			open($sfh, '>:gzip', $selected) or die("Unable to open file, $selected: $!\n");
		} else {
			if($selected == openhandle($selected))
			{
				$sfh = $selected;
			}
			else
			{
				open($sfh, ">$selected") or die("Unable to open file, $selected: $!\n");
			}
		}
	} else {
		open($sfh, ">/dev/null");
	}
	if ($unselected) {
		if ($gzip) {
#			open($ufh, '>:gzip', $unselected) or die("Unable to open file, $unselected: $!\n");
		} else {
			if($unselected == openhandle($unselected))
			{
				$ufh = $unselected;
			}
			else
			{
				open($ufh, ">$unselected") or die("Unable to open file, $unselected: $!\n");
			}
		}
	} else {
		open($ufh, ">/dev/null");
	}

	# LOAD IDS OF MATCHING ROWS
	my $num_targets=0;
	while (<$table>) {
		next if /^#/;
		chomp;
		my @row = split(/\t/);
		my $id = $ignorecase ? uc($row[ $id_col - 1 ]) : $row[ $id_col - 1 ];
		$id=$1 if $id =~ /^(\S+)/;
		$id=$1 if $paired and $id =~ /^(\S+)\/[12]$/;
		if (match($filters, \@row)) {
			# remember this ID
			$ids{$id} = 0; # number of reads with this ID found (counter for paired option)
			++$num_targets;
		}
	}
	unless ($num_targets) {
		# no records match search criteria
		close $table;
		close $sfh if $selected;
		if ($unselected) {
			open(DB, "<$dbfile") or die("Unable to open file, $dbfile: $!\n");
			while (<DB>) {
				print $ufh $_;
				++$n_unselected;
			}
			close DB;
		}
		close $ufh;
		return 0;
	}


	# DETERMINE FILETYPE
	open(DB, "<$dbfile") or die("Unable to open file, $dbfile: $!\n");
	my $format;
	while (<DB>) {
		chomp;
		if (/^#/ or /^$/) { next }
		elsif (/^>/) { $format='fasta' }
		elsif (/^@/) { $format='fastq' }
		else { die "Invalid DB file format" }
		last;
	}
	close DB;

	# GET SEQS
	my $db = $format eq 'fasta' ? FastaDb->new($dbfile) : FastqDb->new($dbfile);
	while (my $rec=$db->next_seq ) {
		my $id = $ignorecase ? uc($rec->id) : $rec->id;
		$id = $1 if $paired and $id =~ /^(\S+)\/[12]$/;
		if (exists($ids{$id})) {
			# selected
			print $sfh $rec->output;
			++$n_selected;
			if (!$paired) {
				delete $ids{$id};
			} else {
				$ids{$id} += 1;
				delete $ids{$id} if $ids{$id} == 2;
			}
		} else {
			# unselected
			print $ufh $rec->output;
			++$n_unselected;
		}
	}
	close $table;
	close $sfh;
	close $ufh;

	# MAKE SURE ALL TARGETS WERE FOUND
	foreach my $id (keys %ids) {
		if ($ids{$id}) {
			delete($ids{$id}); # SOMETIMES INFILES CONTAIN ONLY ONE READ OF PAIR
		}elsif($id eq 'EMPTY') {  # ADDEDBY THO TO ALLOW EMPTY blast results 
			#in workflow using checkempty.pl
			return();
		} else {
			warn("Seq not found: $id\n");
		}
	}
	return ($n_selected,$n_unselected);
}

sub genetree_read_placement
{
	my $newgenes = $_[0];
	my $align    = $_[1];
	my $path     = $_[2];
	my $name     = $_[3];

	#my $newgenes = shift(@ARGV);       #0 new genes to align
	#my $align    = shift(@ARGV);       #1 alignment program to use
	#my $path     = shift(@ARGV);       #2 path to tree and gene data
	#my $name     = shift(@ARGV);       #3 name of gene family

	#define outgroup using hash defined with all 
	my $outgroup = $HoH{$thisgene}{outgroup};

	if($newgenes eq ""){
		# If $newgenes has no hits, do not place the reads, just write a tree with no hits
		print "No hits found. Skipping read placement\n Tree copied to output.\n";
		system "cp $path.tre RAxML_labelledTree.".$thisgene;
		system "cp $path.tre RAxML_originalLabelledTree.".$thisgene;
	}else{

		my $ToALignFile    = $filename . "." . $thisgene . ".toalign.fasta";
		my $AlignedFile    = $filename . "." . $thisgene . ".aligned.fasta";
		my $AlignedPhy     = $filename . "." . $thisgene . ".aligned.phy";

		print "Aligning Hits to known sequences using $align....\n\n";

		#First concatenate fasta files and align
		if($align eq "muscle"){
			print "Align with MUSCLE\n";
			system "echo \"$newgenes\" | cat - $path.fas > $ToALignFile";
			system "muscle -version";
			system "muscle -in $ToALignFile -out $AlignedFile -quiet";
		}
		# MAFFT is annoying, it cuts the sequence ID to 50 characters,
		# replaces hyphens and dots with underscores.
		elsif($align eq "mafft") {
			print "Using MAFFT ";
			qx(mafft --version);
			system "echo \"$newgenes\" | cat - $path.fas > $ToALignFile";
			system "mafft --thread $numThreads --quiet --auto $ToALignFile > $AlignedFile ";
		}
		elsif($align eq "mafftprofile") {
			print "Using MAFFT-profile ";
			qx(mafft --version);
#			system "echo \"$newgenes\" | cat - $path.fas.aligned > $ToALignFile";
			system "mafft --thread $numThreads --add $newgenes --reorder $path.fas.aligned > $AlignedFile";
		}
		elsif($align eq "prank") {
			system "echo \"$newgenes\" | cat - $path.fas > $ToALignFile";
			system "prank -d=$ToALignFile -o=aligned -f=fasta -F";
			system "mv aligned.2.fas $AlignedFile";
		}

		# Convert to phylip format, uses seqConverter.pl
		system path($0)->parent->child("seqConverterG.pl") . " -d$AlignedFile -ope -O$AlignedPhy";

		if( -f "RAxML_classification."                  . $thisgene             and $rebuilTrees) { unlink "RAxML_classification."                  . $thisgene;             }
		if( -f "RAxML_classificationLikelihoodWeights." . $thisgene             and $rebuilTrees) { unlink "RAxML_classificationLikelihoodWeights." . $thisgene;             }
		if( -f "RAxML_entropy."                         . $thisgene             and $rebuilTrees) { unlink "RAxML_entropy."                         . $thisgene;             }
		if( -f "RAxML_info."                            . $thisgene             and $rebuilTrees) { unlink "RAxML_info."                            . $thisgene;             }
		if( -f "RAxML_labelledTree."                    . $thisgene             and $rebuilTrees) { unlink "RAxML_labelledTree."                    . $thisgene;             }
		if( -f "RAxML_originalLabelledTree."            . $thisgene             and $rebuilTrees) { unlink "RAxML_originalLabelledTree."            . $thisgene;             }
		if( -f "RAxML_portableTree."                    . $thisgene . ".jplace" and $rebuilTrees) { unlink "RAxML_portableTree."                    . $thisgene . ".jplace"; }

		print "Placing Hits on gene tree with Maximum Likelihood using Evolutionary Placement Algorithm (EPA) of RAxML...\n";
		system "raxmlHPC -version";
		system "raxmlHPC -f v -s $AlignedPhy -m PROTGAMMALG -t $path.tre -n $thisgene -T $numThreads";
	}

	my $RootedTree = $filename . "." . $thisgene . ".RootedTree";

	#RAxML does not use outgroup information for EPA. Use phyutility to reroot using $outgroup
	my @outgroups = split(',', $outgroup);
	print "Using phyutility to root with OUTGROUPS determined from midpoint rooting\n: @outgroups\n";
	system path($0)->parent->child("phyutility")->child("phyutility.sh") . " -rr -in RAxML_labelledTree.".$thisgene." -out $RootedTree -names @outgroups";

	#Now make tab delimited file to use in tab2trees
	#open treefile to read tree line
	#system "cat RAxML_info.$thisgene";
	open(TREE, "<",$RootedTree) or die "Can't open Rooted Tree File! Rooting requires phyutility";
	my $finaltree;
	while (<TREE>){
		if($_ =~ /\;/m){
			$finaltree = $_;
			chomp($finaltree);
		}
	}
	close TREE;

	$name =~ s/ /_/g;
	chomp($name);
	#remove clade labels
	$finaltree =~ s/\[I\d+\]//g;
	open(TAB, ">>$FinalTreeFile") or die "Can't open File!";
	print TAB $name."\t".$finaltree."\n";
	#print $name."\t".$finaltree."\n";
	close TAB;
}
