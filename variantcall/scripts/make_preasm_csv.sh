#!/usr/bin/bash
ls -1 input/pre-assembled  | perl -p -e 's/((\S+)\_1\.fq\.gz)\n/$2,$1,/' > preasm_samples.csv