#!/usr/bin/python -tt
            
##usage: ./prune_phytab_using_list.py <originalphytab> <bls2remove> <keep|discard> > outfile

#import modules
import sys, os, numpy, re

def read(filename):
  f = open(filename)
  lines = f.readlines()
  # for case where list is an empty file (here, under 20 bytes)
  if os.lstat(sys.argv[2]).st_size < 20:
    for line in lines:
      print(line)
  else:
    bad = open(sys.argv[2]) 
    badlines = bad.readlines()
    badstripped = [line[:-1] for line in badlines]
    str1 = '|'.join(badstripped)
    str2 = '('+str1[:-1]+')'
    pattern = re.compile(str2)
    count=0
    for line in lines:
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


