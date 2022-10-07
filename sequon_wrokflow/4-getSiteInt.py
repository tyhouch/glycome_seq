#leave only the site number
#Operate after 3-siteOnly
import glob

files = glob.glob("./siteOnly_int/*")

for file in files:
    theFile = open(file, "r")
    lines = theFile.readlines()
    theFile.close()
    
    currFile = open(file, "w")
    for line in lines:
        site = ''.join(x for x in line if x.isdigit())
        currFile.write(site)
        currFile.write("\n")
    currFile.close()