#!/usr/bin/python -tt
            
##usage: ./long_branch_finder.py <tabular treelist file> <number of stdevs as cut-off>  > outfile
#import modules
import sys, os, numpy, re

def read(filename):
  f = open(filename)
  lines = f.readlines()
  for eachline in lines:
    line = eachline.split('\t')
    gene = line[0]
    d1 = {}
    d1[gene] = line[1]   #matches genename with its tree
    treetips = re.findall('[a-zA-Z0-9]+(?:_[a-zA-Z0-9]+)?\w*:\d+\.\d+', line[1]) #should be more flexible in recognizing speciesnames in trees
    # treetips = re.findall('[a-zA-Z0-9]+(?:_[a-zA-Z0-9]+)?:\d+\.\d+', line[1])  # makes a list of items like 'spname:bl'
    #treetips = re.findall('[A-Z][a-z]+_[a-z]+:\d+\.\d+', line[1])  # makes a list of items like 'spname:bl'
    d2 = {}
    for i in treetips:
      spbl = i.split(':')
      d2.update({spbl[1] : spbl[0]})  #creates link betwn taxon and its BL
    tipbl = re.findall('\d+\.\d+', str(treetips))
    std = numpy.std([float(i) for i in tipbl])
#    numstd = 3*std
    numstd = int(sys.argv[2])*std
    for i in tipbl:
      if float(i) > float(numstd):
#        print(d2[str(i)] + '\t' + gene  + '\t' + i)
        print(d2[str(i)] + '\t' + gene)
  f.close()

def main():
  read(sys.argv[1])

if __name__ == '__main__':
  main()
