#!/bin/bash
ALLACC=sra_samples.tsv
echo -e "Run\tProject\tSample\tStrain\tLibrary" > $ALLACC

for f in runs/*.tab
do
 RunCol=$(head -n 1 $f | awk -F'\t' ' {
      for(i=1;i < NF;i++) {
         if($i ~ /Run/) { print i }
      }
}')
 ProjCol=$(head -n 1 $f | awk -F'\t' ' {
      for(i=1;i < NF;i++) {
         if($i ~ /BioProject/) { print i }
      }
}')
BioSam=$(head -n 1 $f | awk -F'\t' ' {
      for(i=1;i < NF;i++) {
	IGNORECASE = 1
         if($i ~ /BioSample$/) { print i }
      }
}')
Sample=$(head -n 1 $f | awk -F'\t' ' {
      for(i=1;i < NF;i++) {
         if($i ~ /Sample_Name/) { print i }
      }
}')
Lib=$(head -n 1 $f | awk -F'\t' ' {
      for(i=1;i < NF;i++) {
	IGNORECASE = 1
         if($i ~ /Library_Name/) { print i }
      }
}')
if [ ! $Lib ]; then 
 Lib=$RunCol
fi
 echo "samples are \"-f${RunCol},${ProjCol},${BioSam},${Sample},${Lib}\" for $f"
 tail -n +2 $f | cut -f${RunCol},${ProjCol},${BioSam},${Sample},${Lib} >> $ALLACC
done
