#only leave the glycosylated site line, without other information
#operate after 2-glycoInfoOnly/1-extract

import glob


files = glob.glob("./glycoSiteOnly/*")

for file in files:
    theFile = open(file, "r")
    lines = theFile.readlines()
    theFile.close()
    
    currFile = open(file, "w")
    for idx, line in enumerate(lines):
        if line.find('FT   CARBOHYD') != -1:
            currFile.write(lines[idx])
    currFile.close()