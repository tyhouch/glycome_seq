#generate wget command to get the lectin information txt file
with open('lectins.txt', 'r') as f:
    lectinsName = f.readlines()


res=[]
for sub in lectinsName:
    res.append(sub.replace("\n", ""))

website1="https://rest.uniprot.org/uniprotkb/"
website2=".txt"
for x in res:
    stringOut = "wget "+website1+x+website2
    print(stringOut)
