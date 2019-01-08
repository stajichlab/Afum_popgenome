#!/usr/bin/bash

perl -p -e 's/\t/,/g' sra_samples.tsv > SRA_samples.csv
