#extract the glycosylated file from the lectin data base
import glob
import os
import shutil

glycosylated = []

origin = "./lectin_database/"
target = "./glycosylated_lectin/"

files = glob.glob("./lectin_database/*")
for file in files:
    lectin = open(file)
    lines = lectin.readlines()
    for l in lines:
        if l.find('FT   CARBOHYD') != -1:
            glycosylated.append(file)
            break

glycosylatedName = []

for file in glycosylated:
    lectinName = os.path.basename(file)
    glycosylatedName.append(lectinName)
    
for file_name in glycosylatedName:
   shutil.copy(origin+file_name, target+file_name)