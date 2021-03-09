#!/usr/bin/python -tt
            
##usage: ./pullgoodseqs.py <originalphytab> <bls2remove> <keep|discard> > outfile
#import modules
import sys, os, numpy, re

def read(filename):
  f = open(filename)
  bad = open(sys.argv[2])
  lines = f.readlines()
  badlines = bad.readlines()
  badstripped = [line[:-1] for line in badlines]
  str1 = '|'.join(badstripped)
  str2 = '('+str1[:-1]+')'  
  pattern = re.compile(str2)
  count=0
  for line in lines:
#    line.strip()
     match = pattern.findall(line)
     if match and sys.argv[3] == 'keep':
      print(line)
     if not match and sys.argv[3] == 'discard':
      print(line)
  bad.close()
  f.close()

def main():
  read(sys.argv[1])

if __name__ == '__main__':
  main()
