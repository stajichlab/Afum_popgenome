#!/usr/bin/env python3

import csv, os, re

path="runs"
files = [f for f in os.listdir(path)
             if f.endswith(".table.tab") and
             os.path.isfile(os.path.join(path, f))]

sratable = {}
for f in files:
    with open(os.path.join(path,f)) as tabfile:
        reader = csv.reader(tabfile,delimiter="\t")
        header = next(reader)
        headerdict = {}
        i = 0
        for h in header:
            headerdict[h] = i
            i += 1

        for row in reader:
            ID = row[headerdict["Run"]]

            Center = ""
            if "Center_Name" in headerdict:
                Center = row[ headerdict["Center_Name"]]
                
            Experiment = row[ headerdict["Experiment"]]
            Project    = row[ headerdict["BioProject"]]

            if Project == "PRJNA390160":
                Center = "Vanderbilt"
            elif Project == "PRJNA319359":
                Center = "JPL"

            Library = ""
            if "Library_Name" in headerdict:
                Library = row[ headerdict["Library_Name"]]

            Strain = ""
            if "strain" in headerdict:
                Strain =  row[ headerdict["strain"]]
            elif "Sample_Name" in headerdict:
                Strain = row[ headerdict["Sample_Name"]]
                Strain = re.sub(r'^\d+\s+','',Strain)
                Strain = re.sub(r'Afu\s+','Afu_',Strain)
                Strain = re.sub(r'(A\.|Aspergillus)\s+fumigatus\s+','',Strain)
                Strain = re.sub(r'\s+Sample$','',Strain)
                strainsplit = Strain.split()
                if len(strainsplit) == 2:
                    Strain = strainsplit[0]

            else:
                Strain = Library

            
                
            Biosample = ""
            if "BioSample" in headerdict:
                Biosample = row[ headerdict["BioSample"]]
            else:
                print("no biosample in file %f" %(f))
                
            if ID in sratable:
                print("already seen %s checking file %s" % (ID,f))
            sratable[ID] = [ID,Strain,Biosample,
                            Center,Experiment,Project]


with open('sra_samples.tsv', 'w') as tsvfile:
    writer = csv.writer(tsvfile,delimiter="\t",quoting=csv.QUOTE_MINIMAL,lineterminator='\n')
    writer.writerow(["RunAcc","Strain","BioSample",
                     "Center","Experiment","Project"])
    for srarun in sorted(sratable):
        writer.writerow(sratable[srarun])
#Cols
#Run -> Run
#BioSample -> BioSample
#Library_Name -> Library
#strain -> Strain
#sample_name -> Strain
#Experiment -> Experiment
#Center -> Center_Name
#BioProject -> Project
