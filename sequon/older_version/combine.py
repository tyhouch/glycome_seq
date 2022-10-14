import glob

files = glob.glob("./glycosylated_site_sequence/*")

for file in files:
    theFile = open(file, "r")
    lines = theFile.readlines()
    theFile.close()
    for line in lines:
        line = line.replace("\n", "")
        print(line)