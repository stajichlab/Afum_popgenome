#!/usr/bin/env python3

import csv,sys,re
if len(sys.argv) < 2:
    exit("Expected 1 argument: .csv file")
infile=sys.argv[1]
with open(infile,"r") as csvfile:
        reader = csv.reader(csvfile,delimiter=",")
        for row in reader:
            print("\t".join(row))
