#leave only the glycosylated site information
#operate after 1-extract
import glob


files = glob.glob("./glycoInformationOnly/*")

for file in files:
    theFile = open(file, "r")
    lines = theFile.readlines()
    theFile.close()
    
    currFile = open(file, "w")
    for idx, line in enumerate(lines):
        if line.find('FT   CARBOHYD') != -1:
            currFile.write(lines[idx])
            currFile.write(lines[idx+1])
            currFile.write(lines[idx+2])
    currFile.close()