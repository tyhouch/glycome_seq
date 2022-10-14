#re-arrange the sequence to get the site easier
import glob

files = glob.glob("./glycosylated_site_sequence/*")

for file in files:
    theFile = open(file, "r")
    lines = theFile.readlines()
    theFile.close()
    
    currFile = open(file, "w")
    for idx, line in enumerate(lines):
        if idx == 0:
            currFile.write(line)
        if idx != 0:
            line = line.replace("\n", "")
            currFile.write(line)
    currFile.close