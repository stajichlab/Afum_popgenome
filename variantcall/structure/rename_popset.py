#!/usr/bin/env python3

import csv

lookup = {}

with open("../tip_names.tab","r") as names:
    for line in names:
        line = line.strip()
        row = line.split("\t")
        lookup[row[0]] = row[1]

with open("AfumAf293.Run2.popset","r") as popset:
    for line in popset:
        line = line.strip()
        if line in lookup:
            print("%s\t%s"%(lookup[line],line))
        else:
            print("%s\t%s"%(line,line))
