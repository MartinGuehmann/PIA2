# PIA2
## Modified version of the Phylogenetically Informed Annotation tool (Speiser et al., 2014)


### This collection of scripts is written in bash, perl, and python.
###
### Make sure dependencies are installed and configured before running.
### Dependencies include: BLAST, MAFFT, R (w/ ape and phytools packages), Perl, BioPerl, NumPy, SciPy, statsmodels (for python), Java, EPA-ng, Gappa, USEARCH, … and any dependencies these might have.
### 
### Special thanks to Dr. Todd Oakley's lab at UCSB for the original galaxy scripts on 
### which we have based this pipeline. Original versions can be found at:
### https://bitbucket.org/osiris_phylogenetics/osiris_phylogenetics
### and
### https://bitbucket.org/osiris_phylogenetics/pia
###
###  Jorge L. Perez-Moreno, Danielle DeLeo, Heather D. Bracken-Grissom
###  CRUSTOMICS Lab at Florida International University 04/05/2017
------------------------------------------------------------------------------------------

## Setting PIA up is pretty simple:

	# - Clone this repository to /some/directory/
	# - Clone https://github.com/peterjc/pico_galaxy.git to /some/directory/
		- pico_galaxy is a dependency
	# - Clone and build https://github.com/lczech/gappa.git, copy gappa to /usr/bin/
	# - Clone and build https://github.com/Pbdas/epa-ng.git, copy epa-ng to /usr/bin/
	# - You may edit numThreads in run_pia.pl if you want to use less CTP cores than you have available on your system
	# - You may edit the -m option in the system call raxmlHPC in pia.pl if you need to customize the model of evolution for EPA-ng. 

	# Further things to install (list probably incomplete)
	# - Install ncbi blast+ (sudo apt-get install ncbi-blast+)
	# - Install statsmodels (https://www.statsmodels.org/stable/install.html)
	# - Install R
		- Install packets ape
		- Install phytools requires (Ubuntu, Debian)
			- sudo apt install libmagick++-dev
			- sudo apt install libcurl4-openssl-dev
	# - Install usearch (https://www.drive5.com/usearch/download.html)
	# - Install CD-Hit (https://github.com/weizhongli/cdhit)

## The script run_pia.sh will run PIA and post-PIA on all fasta files in the working directory. 
## Both run_pia.sh and post_pia.sh should be edited before running to adjust parameters.



------------------------------------------------------------------------------------------

The original PIA (Phylogenetically Informed Annotation) is a set of tools for the Galaxy Bioinformatics Platform. In general, PIA uses BLAST, an alignment program, and standard RAxML's read evolutionary placement algorithm (RAxML-EPA) to put unknown sequences into pre-calculated phylogenetic trees.
We provide 102 genes called LIT (Light Interaction Toolkit) - vision genes like phototransduction genets - for use in PIA.

This version of PIA is designed to run from Bash and replaces RAxML-EPA with EPA-ng to run faster.

License:
All original source code for PIA is available under the MIT license (http://opensource.org/licenses/mit-license.html). See below:
The MIT License (MIT)
Copyright (c) 2014 Speiser et al

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

