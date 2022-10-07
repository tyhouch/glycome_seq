#re-arrange the sequence to get the site easier
import glob

files = glob.glob("./glycosylated_site_sequence/*")

for file in files:
    theFile = open(file, "r")
    lines = theFile.readlines()
    theFile.close()
    
    currFile = open(file, "w")
    for idx, line in enumerate(lines):
        line = line.replace("\n", "")
        if idx != 0:
            currFile.write(line)
    currFile.close