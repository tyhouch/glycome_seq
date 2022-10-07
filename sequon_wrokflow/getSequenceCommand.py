#generate wget command to get the sequence fasta file
import glob
import os

files = glob.glob("./siteOnly/*")

glycosylatedName = []

for file in files:
    lectinName = os.path.basename(file)
    glycosylatedName.append(lectinName)


res=[]
for sub in glycosylatedName:
    sub = sub.replace(".txt", "")
    sub = sub.replace("\n", "")
    res.append(sub) 
    
  
website1="https://rest.uniprot.org/uniprotkb/"
website2=".fasta"
for x in res:
    stringOut = "wget "+website1+x+website2
    print(stringOut)
