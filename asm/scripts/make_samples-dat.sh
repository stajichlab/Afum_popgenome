#!/bin/bash
pushd input
ls *_R1.fq.gz | perl -p -e 's/_R1\.fq\.gz/\tAscomycota/' > ../samples2.dat
popd
