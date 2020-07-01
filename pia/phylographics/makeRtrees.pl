#!/usr/bin/perl -w

#This script generates an R script to print trees to a pdf file
#input is a table with treename<tab>newick tree
use strict;

my $filename = $ARGV[0];
my $outfile = $ARGV[1];
open FILE, $filename or die $!;
my $treetype = $ARGV[2];
my $extiplabels = $ARGV[3];
my $options;
my $labeltaxfile = $ARGV[4];
my $query = $ARGV[5];
my $midroot = $ARGV[6];
my %labelhash;
my $genecount=0;
my @genes;
my $markquery;
my $midpoint;

if($query eq 'yes'){
	$markquery = 1;
}else{
	$markquery = 0;
}
unless($labeltaxfile eq 'None'){
	open LABELFILE, $labeltaxfile or die $!;
	while (<LABELFILE>) {
	        chomp;
	        #get a line from the data file
	        my $currentinput = "$_";
		if($currentinput =~ /\t/){ 
			my @splitline = split(/\t/);
			my $speciesname= $splitline[0];
			$speciesname = "'".$speciesname."'";
			my $treename = $splitline[1];
			if(exists $labelhash{$treename}){
				push @{ $labelhash{$treename} }, $speciesname;
			}else{
				push @{ $labelhash{$treename} }, $speciesname;
				#$labelhash{$treename} = $speciesname;
				$genecount ++;
				push @genes, $treename;
			}	
		}
	}

}#end unless


if($extiplabels eq 'yes'){
	$options = ", show.tip.label=FALSE";
}else{
	$options = ", show.tip.label=TRUE";
}

print "require(ape);\n";
print "require(phytools);\n";
#print "require(phangorn);\n";
print "pdf(file='$outfile', paper = 'a4', width = 0, height = 0);\n";



while (<FILE>) {
        chomp;
        #get a line from the data file
        my $currentinput = "$_";
	my @splitline = split(/\t/);
	my $treename= $splitline[0];
	my $tree = $splitline[1];
	my $labelsvector;

	#print the R commands to make tree graphics
        print "raw_tree <- read.tree(text = '$tree');\n";
	print "raw_tree\$edge.length[ is.na(raw_tree\$edge.length) ] <- 0 \n";

	#Check if large tree, then make text size smaller
	print "thetips <- raw_tree\$tip.label \n";
	print "numtips <- length(thetips) \n";

#Make text smaller for trees with many tips
        print "if(numtips<45){numtips <- 45} \n";
	print "cexval <- 233*numtips^-1.5 \n";
	print "cexval <- round(cexval, 2) \n";
        print "if(cexval<0.05){cexval <- 0.05} \n";

	print "plot(raw_tree, cex=cexval, edge.width=0.1, type='$treetype' $options) \n";

        print "title('Tree File: $treename');\n";

#Add taxon labels, if optional file present and if labels exist for tree
	if(exists $labelhash{$treename}){
		$labelsvector = join ",", @{ $labelhash{$treename} };
		$labelsvector = "tolabel <- c(".$labelsvector.")";
		print "thetips <- raw_tree\$tip.label \n";
		print $labelsvector."\n";
		print "labels <- match(tolabel,thetips) \n";
		print "tiplabels(tip=labels, pch=21, cex=cexval) \n";
	}


#Add taxon labels if gene name contains QUERY - for readplacement
	if($markquery == 1){
		print "thetips <- raw_tree\$tip.label \n";
		print "qlabels <- grep(\'QUERY\',thetips) \n";
		print "tiplabels(tip=qlabels, pch=21, cex=cexval) \n";
		print "l1labels <- grep(\'LANDMARK1\',thetips) \n";
		print "tiplabels(tip=l1labels, pch=15, cex=cexval, col='red') \n";
	}
}

print "dev.off();\n";
close FILE;

sub printMidpointRoot(){

print <<END;
midpoint.root<-function(tree){
	D<-cophenetic(tree)
	dd<-max(D)
	ii<-which(D==dd)[1]
	ii<-c(ceiling(ii/nrow(D)),ii%%nrow(D))
	if(ii[2]==0) ii[2]<-nrow(D)
	spp<-rownames(D)[ii]
	nn<-which(tree\$tip.label==spp[2])
	tree<-reroot(tree,nn,tree\$edge.length[which(tree\$edge[,2]==nn)])
	aa<-getAncestors(tree,which(tree\$tip.label==spp[1]))
	D<-c(0,dist.nodes(tree)[which(tree\$tip.label==spp[1]),aa])
	names(D)[1]<-which(tree\$tip.label==spp[1])
	i<-0
	while(D[i+1]<(dd/2)) i<-i+1
	tree<-reroot(tree,as.numeric(names(D)[i]),D[i+1]-dd/2)
	tree
}

## function gets ancestor node numbers
## written by Liam J. Revell

getAncestors<-function(tree,node,type=c("all","parent")){
	type<-type[1]
	if(type=="all"){
		aa<-vector()
		rt<-length(tree\$tip.label)+1
		currnode<-node
		while(currnode!=rt){
			currnode<-getAncestors(tree,currnode,"parent")
			aa<-c(aa,currnode)
		}
		return(aa)
	} else if(type=="parent"){
		aa<-tree\$edge[which(tree\$edge[,2]==node),1]
		return(aa)
	} else stop("do not recognize type")
}
END
}

#Testing hash arrays
#my %nums;
#my $test='odd';
#for my $n (4,5,6,10) {
#    if ($n % 2) {
#        push @{ $nums{$test} }, $n;
#    } else {
#        push @{ $nums{even} }, $n;
#    }
#}
#
#print join ', ', @{ $nums{even} };
#print "\n\n";
