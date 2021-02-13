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
            if "Run" in headerdict:
                ID = row[headerdict["Run"]]
            else:
                print("cannot find Run in file {}".format(tabfile))
                exit()

            Center = ""
            if "Center_Name" in headerdict:
                Center = row[ headerdict["Center_Name"]]
                
            Experiment = row[ headerdict["Experiment"]]
            Project    = row[ headerdict["BioProject"]]

            if Project == "PRJNA390160":
                Center = "Vanderbilt"
            elif Project == "PRJNA319359":
                Center = "JPL"
            elif Project == "PRJNA477519":
                Center = "ICIII"
            elif Project == "PRJEB1497":
                Center = "Manchester"
            elif ID == "SRR2954803":
                Center = "LA_MOLINA"
            elif Project == "PRJNA528395":
                Center = "Radboud"
            elif Project == "PRJNA638646":
                Center = "UMassAmherst"

            Library = ""
            if "Library_Name" in headerdict:
                Library = row[ headerdict["Library_Name"]]

            Strain = ""
            if "strain" in headerdict:
                Strain =  row[ headerdict["strain"]]
            elif "Strain" in headerdict:
                Strain = row[ headerdict["Strain"]]
            elif "Sample_Name" in headerdict:
                Strain = row[ headerdict["Sample_Name"]]
                Strain = re.sub(r'^\d+\s+','',Strain)
                Strain = re.sub(r'Afu\s+','Afu_',Strain)
                Strain = re.sub(r'(A\.|Aspergillus)\s+fumigatus\s+','',Strain)
                Strain = re.sub(r'\s+Sample$','',Strain)
                Strain = re.sub(r'AF293\s+','AF293-',Strain)
                Strain = re.sub(r'Af10\s+.+','Af10',Strain)
                strainsplit = Strain.split()
                if len(strainsplit) == 2:
                    Strain = strainsplit[0]

            else:
                Strain = Library

            Strain = re.sub(r'AF293\s+','AF293-',Strain)

            Biosample = ""
            if "BioSample" in headerdict:
                Biosample = row[ headerdict["BioSample"]]
            else:
                print("no biosample in file %f" %(f))

            # hardcode one sample to deal with same name 2x
            if Biosample == "SAMEA2052099":
                Strain = Strain + "-MCH"
            #black LIST these strains
            if Biosample == "SAMEA2052103" or Biosample == "SAMEA2051883":  # these are A.nidulans 
                continue
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
