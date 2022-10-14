#Get the sequon on each glycosylated site
#Operate after getSiteInt-4 and arrangeSequence

import glob
import os

sites = glob.glob("./siteOnly_int/*")
files = glob.glob("./glycosylated_site_sequence/*")

for file in files:
    fileName = os.path.basename(file)
    fileName = fileName.replace(".fasta",".txt")
    glycoSites = []
    for site in glob.glob("./siteOnly_int/"+fileName):
        siteFile = open(site, "r")
        glycoSites = siteFile.readlines()
        siteFile.close()
     
    for i in range(0, len(glycoSites)):
        glycoSites[i] = glycoSites[i].removesuffix("\n")
        glycoSites[i] = int(glycoSites[i])
         
    theFile = open(file, "r")
    lines = theFile.readlines()
    theFile.close()
    
    currFile = open(file, "w")
    for idx, line in enumerate(lines):
        if idx == 0:
            currFile.write(line)
        else:
            for site in glycoSites:
                segment = line[site-1-3 : site-1+4]#-1 to find actual index
                currFile.write(segment)
                currFile.write("\n")
    currFile.close()
