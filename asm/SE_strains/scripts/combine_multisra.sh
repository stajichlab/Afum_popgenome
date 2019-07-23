#!/usr/bin/bash

if [ ! -f CEA10-MCH_R1.fq.gz ]; then
   zcat ERR232426_R1.fq.gz ERR232428_R1.fq.gz | gzip -c > CEA10-MCH_R1.fq.gz
fi
rm ERR232426_R1.fq.gz ERR232428_R1.fq.gz
if [ ! -f CEA10-MCH_R2.fq.gz ]; then
   zcat ERR232426_R2.fq.gz ERR232428_R2.fq.gz | gzip -c > CEA10-MCH_R2.fq.gz
fi
rm ERR232426_R2.fq.gz ERR232428_R2.fq.gz
if [ ! -f AF300_R1.fq.gz ]; then
   zcat ERR232419_R1.fq.gz ERR232420_R1.fq.gz | gzip -c > AF300_R1.fq.gz
fi
rm ERR232419_R1.fq.gz ERR232420_R1.fq.gz
if [ ! -f AF300_R2.fq.gz ]; then
   zcat ERR232419_R2.fq.gz ERR232420_R2.fq.gz | gzip -c > AF300_R2.fq.gz
fi
rm ERR232419_R2.fq.gz ERR232420_R2.fq.gz
if [ ! -f CEA10-2_R1.fq.gz ]; then
   zcat ERR232423_R1.fq.gz ERR232427_R1.fq.gz | gzip -c > CEA10-2_R1.fq.gz
fi
rm ERR232423_R1.fq.gz ERR232427_R1.fq.gz
if [ ! -f CEA10-2_R2.fq.gz ]; then
   zcat ERR232423_R2.fq.gz ERR232427_R2.fq.gz | gzip -c > CEA10-2_R2.fq.gz
fi
rm ERR232423_R2.fq.gz ERR232427_R2.fq.gz
