import tkinter as tk
from tkinter import filedialog
import os
import glob
import shutil
from contextlib import redirect_stdout
import pandas as pd
import numpy as np
from collections import defaultdict


def popMessage(message):
    tk.messagebox.showinfo("message",message)
    
def createFolder(folder_parent, folder_name):
    target = folder_name
    target_parent = folder_parent
    if os.path.exists(target_parent+target):
        shutil.rmtree(target_parent+target)#remove path
    path = os.path.join(target_parent, target)
    os.mkdir(path)
    return target_parent+target

    
def readLines(file):
    theFile = open(file, "r")
    lines = theFile.readlines()
    theFile.close()
    return lines
    
   
def extractGlycoLectins() :
    popMessage("please select the folder with all lectin informations")
    
    source = filedialog.askdirectory()
    target = "glycosylated_lectin"
    target_parent = "./"
    if os.path.exists(target_parent+target):
        shutil.rmtree(target_parent+target)#remove path
    path = os.path.join(target_parent, target)
    os.mkdir(path)
    targetFolder = target_parent+target+"/"
    
    glycosylated = []
    files = glob.glob(source+"/*")
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
       shutil.copy(source+"/"+file_name, targetFolder+file_name)
    
def getGlycoInfo():
    popMessage("Please select the folder with the lectin informations")
    source = filedialog.askdirectory()
    target = "glycoInformationOnly"
    target_parent = "./"
    if os.path.exists(target_parent+target):
        shutil.rmtree(target_parent+target)#remove path
    path = os.path.join(target_parent, target)
    os.mkdir(path)
    targetFolder = target_parent+target+"/"
    
    file_names =[]
    sourceFiles = glob.glob(source+"/*")
    for file in sourceFiles:
        currFileName = os.path.basename(file)
        file_names.append(currFileName)
    for file_name in file_names:
        shutil.copy(source+"/"+file_name, targetFolder+file_name)
    
    targetFiles = glob.glob(targetFolder+"*")
    for file in targetFiles:
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
    
def getGlycoSitesNum():
    popMessage("Please select the folder with lectin information")
    source = filedialog.askdirectory()
    target = "glycoSites_int"
    target_parent = "./"
    if os.path.exists(target_parent+target):
        shutil.rmtree(target_parent+target)#remove path
    path = os.path.join(target_parent, target)
    os.mkdir(path)
    targetFolder = target_parent+target+"/"
    
    file_names =[]
    sourceFiles = glob.glob(source+"/*")
    for file in sourceFiles:
        currFileName = os.path.basename(file)
        file_names.append(currFileName)
    for file_name in file_names:
        shutil.copy(source+"/"+file_name, targetFolder+file_name)
    
    targetFiles = glob.glob(targetFolder+"*")
    for file in targetFiles:
        theFile = open(file, "r")
        lines = theFile.readlines()
        theFile.close()
    
        currFile = open(file, "w")
        for idx, line in enumerate(lines):
            if line.find('FT   CARBOHYD') != -1:
                site = ''.join(x for x in line if x.isdigit())
                currFile.write(site)
                currFile.write("\n")
        currFile.close()
        
def getGlcoSequence():
    popMessage("please first select the folder with all lectin sequences, then select the folder with their glycosylated site numbers")
    source = filedialog.askdirectory()
    siteNum = filedialog.askdirectory()
    target = "glycoSitesSequence"
    target_parent = "./"
    if os.path.exists(target_parent+target):
        shutil.rmtree(target_parent+target)#remove path
    path = os.path.join(target_parent, target)
    os.mkdir(path)
    targetFolder = target_parent+target+"/"
    
    file_names =[]
    sourceFiles = glob.glob(source+"/*")
    for file in sourceFiles:
        currFileName = os.path.basename(file)
        file_names.append(currFileName)
    for file_name in file_names:
        shutil.copy(source+"/"+file_name, targetFolder+file_name)
    
    targetFiles = glob.glob(targetFolder+"*")
    for file in targetFiles:
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
    
    for file in targetFiles:
        fileName = os.path.basename(file)
        fileName = fileName.replace(".fasta",".txt")
        glycoSites = []
        for site in glob.glob(siteNum+"/"+fileName):
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
        
def combineSequence():
    popMessage("Please select the folder with the files you want to combine")
    source = filedialog.askdirectory()
    sourceFiles = glob.glob(source+"/*")
    
    for file in sourceFiles:
        theFile = open(file, "r")
        lines = theFile.readlines()
        theFile.close()
        f = open(source + "/../Combined"+os.path.basename(source)+".txt", "a")
        for line in lines:
            f.write(line)
        f.close()
    

def makeLocalDatabase():
    popMessage("please select the combined sequence .fasta file")
    source = filedialog.askopenfilename()
    sourceName = os.path.basename(source)
    
    
    makeDatabase = os.system("makeblastdb -in " + sourceName + " -out Database -dbtype prot -parse_seqids" )


def localBlast():
    popMessage("""Please make sure you have the local database in the folder select the 
               glycosylated_lectin_sequence folder and make sure the glycosylated_lectin_sequence folder in the current directory""")
    
    source = filedialog.askdirectory()
    sourceFiles = glob.glob(source+"/*")
    
    target = "blast_output_outfmt6"
    target_parent = "./"
    if os.path.exists(target_parent+target):
        shutil.rmtree(target_parent+target)#remove path
    path = os.path.join(target_parent, target)
    os.mkdir(path)
    targetFolder = target_parent+target+"/"
    
    for file in sourceFiles:
        currFileName = os.path.basename(file)
        currFileNameInTxt = currFileName.replace(".fasta", ".txt")
        doLocalBlast = os.system("blastp -query ./glycosylated_lectin_sequence/"+ currFileName + " -db Database -out " 
                                 +targetFolder+currFileNameInTxt +" -outfmt 6")

def defaultLectinName():
    popMessage("""Please select the glycosylated_lectin_sequence folder or any other folder containing all the glycosylated sequences""")
    source = filedialog.askdirectory()
    allFileName = []
    sourceFiles = glob.glob(source+"/*")
    for file in sourceFiles:
        currFileName = os.path.basename(file).removesuffix(".fasta")
        allFileName.append(currFileName)
    
    target = "default_lectin_names"
    target_parent = "./"
    if os.path.exists(target_parent+target):
        shutil.rmtree(target_parent+target)#remove path
    path = os.path.join(target_parent, target)
    os.mkdir(path)
    targetFolder = target_parent+target+"/"
    
    for i in allFileName:
        currFile = open(targetFolder+i+".txt", "w")
        for j in range(len(allFileName)):
            currFile.write(allFileName[j]+"\n")


def extractEValue():
    popMessage("Plase select the output folder from the local blast")
    source = filedialog.askdirectory()
    sourceFiles = glob.glob(source+"/*")
    
    target = "e-value_only"
    target_parent = "./"
    if os.path.exists(target_parent+target):
        shutil.rmtree(target_parent+target)#remove path
    path = os.path.join(target_parent, target)
    os.mkdir(path)
    targetFolder = target_parent+target+"/"
    
    file_names =[]
    for file in sourceFiles:
        currFileName = os.path.basename(file)
        file_names.append(currFileName)
    for file_name in file_names:
        shutil.copy(source+"/"+file_name, targetFolder+file_name)
    
    targetFiles = glob.glob(targetFolder+"*")
    for file in targetFiles:
        theFile = open(file, "r")
        lines = theFile.readlines()
        theFile.close()
    
        currFile = open(file, "w")
        for line in lines:
            elements = line.split(",")
            currFile.write(elements[1])
            currFile.write(" ")
            currFile.write(elements[10])
            currFile.write("\n")
        currFile.close()
          
def lectinEValueMap():
    popMessage("Plase first select the default lectin names folder then select the e-values folder")
    lectinNames = filedialog.askdirectory()
    eValues = filedialog.askdirectory()
    lectinNamesFiles = glob.glob(lectinNames+"/*")
    
    
    target = "lectin_e-value_map"
    target_parent = "./"
    if os.path.exists(target_parent+target):
        shutil.rmtree(target_parent+target)#remove path
    path = os.path.join(target_parent, target)
    os.mkdir(path)
    targetFolder = target_parent+target+"/"
    
    file_names =[]
    for file in lectinNamesFiles:
        currFileName = os.path.basename(file)
        file_names.append(currFileName)
    for file_name in file_names:
        shutil.copy(lectinNames+"/"+file_name, targetFolder+file_name)
        
    targetFiles = glob.glob(targetFolder+"*")
    for file in targetFiles:
        theFile = open(file, "r")
        lines = theFile.readlines()
        theFile.close()
        
        
        eValueFile = open("./e-value_only/"+os.path.basename(file), "r")
        eValueLines = eValueFile.readlines()
        eValueFile.close()
        
        eValues = {}
        for line in eValueLines:
            line = line.removesuffix("\n")
            elements = line.split(" ")
            eValues[elements[0]] = elements[1]
        currFile = open(file, "w")
        for line in lines:
            lectinName = line.removesuffix("\n")
            if lectinName in eValues:
                currFile.write(lectinName+" "+
                               eValues.get(lectinName)+"\n")
            else:
                currFile.write(lectinName+" "+
                               "-1"+"\n")
        currFile.close()

def initializeList(size):
    toReturn = []
    for i in range(size):
        toReturn.append("\0")
    return toReturn
    
    
def splitFiles():
    source = filedialog.askdirectory()
    sourceFiles = glob.glob(source+"/*")
    
    fileNames = initializeList(1000)
    for idx, file in enumerate(sourceFiles):
        fileNames[idx] = os.path.basename(file)
    
    splits= defaultdict(list)   
    for i in range(20):
        for j in range(50):
            if fileNames[i*50+j] == "\0":
                break
            else:
                splits[i].append(fileNames[i*50+j])
                
    for i in range(20):
        target = "splits"+str(i)
        target_parent = "./splits/"
        if os.path.exists(target_parent+target):
            shutil.rmtree(target_parent+target)#remove path
        path = os.path.join(target_parent, target)
        os.mkdir(path)
        targetFolder = target_parent+target+"/"
        
        for fileName in splits[i]:
            shutil.copy(source+"/"+fileName, targetFolder+fileName)
    
def combineMatrix():
    source = filedialog.askdirectory()
    sourceFiles = glob.glob(source+"/*")
    
    for file in sourceFiles:
        theFile = open(file, "r")
        lines = theFile.readlines()
        theFile.close()
        f = open("./Combined"+os.path.basename(source)+"Matrix.txt", "a")
        for idx, line in enumerate(lines):
            if idx == len(lines)-1: 
                elements = []
                elements = line.removesuffix("\n").split(" ")
                f.write(elements[1])
            else:
                elements = []
                elements = line.removesuffix("\n").split(" ")
                f.write(elements[1]+" ")
        f.write("\n")
        f.close()
               
def MCLResultProcess():
    popMessage("Please select the MCL output file")
    source = filedialog.askopenfilename()
    
    file  = open(source, "r")
    lines = file.readlines()
    file.close()
    
    f = open("./Processed_"+os.path.basename(source)+".txt", "a")
    for line in lines:
        lectinNames = []
        line = line.removesuffix("\n")
        elements = line.split("\t")
        for element in elements:
            #step1: get rid of any prefix/surfix
            if "|" in element:
                subStrings = element.split("|")
                element = subStrings[1]
            #step2: remove dupilicates     
            if lectinNames.count(element) == 0:
                lectinNames.append(element)
        #rewite in a new file    
        for idx, line in enumerate(lectinNames):
            if idx == len(lectinNames)-1: 
                f.write(lectinNames[idx])
            else:
                f.write(lectinNames[idx]+" ")
        f.write("\n")
    
    f.close()
                
    
def getclusters():
    popMessage("please select the cluster file")
    source = filedialog.askopenfilename()#cluster file
    popMessage("please select the folder with all glycosylated lectin sequence")
    sequenceFolder = filedialog.askdirectory()#lectin sequence folder
    
    sequenceFiles = glob.glob(sequenceFolder+"/*")
    
    clusterFile = open(source, "r")
    clusterLines = clusterFile.readlines()
    clusterFile.close()
    
    clusters = defaultdict(list)
    
    for idx, line in enumerate(clusterLines):
        line = line.removesuffix("\n")
        elements = line.split(" ")
        for element in elements:
            clusters[idx].append(element)
    
    createFolder("./", "cluster")
    for x in range(len(clusters)):
        targetFolder = createFolder("./cluster/", "cluster"+str(x))
        for lectin in clusters[x]:
            shutil.copy(sequenceFolder + "/" + lectin + ".fasta", 
                        targetFolder+ "/" + lectin + ".fasta")
        
        
    
def test():# to fit CLANS
    source = filedialog.askopenfilename()
    
    file  = open(source, "r")
    lines = file.readlines()
    file.close()
    
    f = open("./Processed"+os.path.basename(source)+".txt", "a")
    for line in lines:
        lectinNames = []
        line = line.removesuffix("\n")
        elements = line.split("\t")
        for element in elements:
            #step1: get rid of any prefix/surfix
            if "|" in element:
                subStrings = element.split("|")
                element = subStrings[1]
            #step2: remove dupilicates     
            lectinNames.append(element)
        #rewite in a new file    
        for idx, line in enumerate(lectinNames):
            if idx == len(lectinNames)-1: 
                f.write(lectinNames[idx])
            else:
                f.write(lectinNames[idx]+" ")
        f.write("\n")
    
    f.close()    


def test2():# to fit CLANS
    source = filedialog.askopenfilename()
    
    file  = open(source, "r")
    lines = file.readlines()
    file.close()
    
    f = open("./Processed"+os.path.basename(source)+".txt", "a")
    num = 0
    for line in lines:
        if ">" in line:
            line = ">seq"+ str(num)+"\n"
            num= num +1
        f.write(line)
    
    f.close()    




"""
This function counts the how many times the specif keyword appears in a given file
and also tells the specific ID of that lectin and the cluster it belongs

input1: the file containing the name(or keywords) we want to search for
input2: the file of combined glycosylated sequence
input3: the processed MCL cluster result file with only the UniProt id for each lectins

output: a file in format: name counts ['id cluster', ...]
                        i.e.,Concanavalin 2 ['P02866 2', 'P02866 3', 'P14894 2', 'P14894 3']
"""
def countLectins():
    popMessage("Please select the file with names you want to search")
    nameFile = filedialog.askopenfilename()
    popMessage("Please select the combined glycosylated sequence file")
    sequenceFile = filedialog.askopenfilename()
    popMessage("Please select the processed MCL result file(with only UniProt id of the lectin)")
    clusterFile = filedialog.askopenfilename()
    
    nameOpenFile = open(nameFile, "r")
    nameLines = nameOpenFile.readlines()
    nameOpenFile.close()
    
    sequenceOpenFile = open(sequenceFile, "r")
    sequenceLines = sequenceOpenFile.readlines()
    sequenceOpenFile.close()
    
    clusterOpenFile = open(clusterFile, "r")
    clusterLines = clusterOpenFile.readlines()
    clusterOpenFile.close()
    
    clusters = defaultdict(list)
    
    for idx, line in enumerate(clusterLines):
        line = line.removesuffix("\n")
        elements = line.split(" ")
        for element in elements:
            clusters[idx].append(element)
    
    f = open("./counts_"+os.path.basename(nameFile), "a")
    for name in nameLines:
        count = 0
        name = name.removesuffix("\n")
        id = []
        for line in sequenceLines:
            if name in line:
                count = count + 1
                elements = line.split("|")
                for key in range(len(clusters)):
                    if elements[1] in clusters[key]:
                        id.append(elements[1]+ " " + str(key))
        f.write(name + " " + str(count) + " "+ str(id) +"\n")
    f.close()



"""
This function will extract the sequence and glyco sites information with a given
list of names
    
input1: the file containing the names(keywords) we want to search for
input2: the folder containg the glycosylated lectin sequences
input3: the folder containg the glyco sites(in int)
    
output: a folder under current directory; each name has a separate folder with their sequence & glyco sites
"""
def folderOfSequenceAndSites():
    popMessage("Please select the file with the lectin names you want to extract the sequence and glyco sites")
    nameFile = filedialog.askopenfilename()
    popMessage("Please select the folder with all sequence files")
    sequenceFolder= filedialog.askdirectory()
    popMessage("please select the folder with glyco sites in integer(the glycoSite_int folder)")
    glycoSiteFolder = filedialog.askdirectory()
    
    sequenceFiles = glob.glob(sequenceFolder+"/*")
    glycoSiteFiles = glob.glob(glycoSiteFolder+"/*")
    
    nameLines = readLines(nameFile)
    
    createFolder("./", "57Lectins")
    targetFolder = "./57Lectins/"
    
    for name in nameLines:
        name = name.removesuffix("\n")
        
        id = []
        for sequenceFile in sequenceFiles:
            sequenceLines = readLines(sequenceFile)
            for line in sequenceLines:
                if "|" in line:
                    if name in line:
                        elements = line.split("|")
                        id.append(elements[1])
        
        if len(id) == 0:
            continue
        
        lectinFolder = createFolder(targetFolder, name)
        seqTargetFolder = createFolder(lectinFolder+"/", "sequence")
        glycoTargetFolder = createFolder(lectinFolder+"/", "glyco_sites")            
        if len(id) > 0:
            for i in range(len(id)):
                shutil.copy(sequenceFolder + "/" + id[i] + ".fasta", 
                            seqTargetFolder+"/"+ id[i] + ".fasta")
                shutil.copy(glycoSiteFolder + "/" + id[i] + ".txt", 
                            glycoTargetFolder+"/"+ id[i] + ".txt")
                
        
"""
Combine fileswith file name
"""            
def combineWFileName():
    popMessage("Please select the folder with the files you want to combine")
    source = filedialog.askdirectory()
    sourceFiles = glob.glob(source+"/*")
    
    f = open(source + "/../Combined"+os.path.basename(source)+".txt", "a")
    for file in sourceFiles:
        theFile = open(file, "r")
        lines = theFile.readlines()
        theFile.close()
        #file name
        fileName = os.path.basename(file)
        if "." in fileName:
            elements = fileName.split(".")
            fileName = elements[0]
        f.write(fileName)
        for line in lines:
            f.write(line)
        f.close()
                
        
"""
GUI for above functions 
"""        
root = tk.Tk()
root.title("Sequencing_Glycome")
root.geometry('480x700')  

preparation = tk.Label(root, text = "Preparation")
preparation.place(x= 0, y= 0)

extractLectin = tk.Button(root, text = "Extract Glycosylated Lectins",
                          padx = 20, pady = 10, fg= "#000000", bg = "white", command=extractGlycoLectins)
extractLectin.place(x = 0, y = 20)

sequonProcess = tk.Label(root, text = "Sequon")
sequonProcess.place(x = 0, y = 70)

GetGlycoInfo= tk.Button(root, text = "Get Glycosylated Site Inofrmation",
                          padx = 20, pady = 10, fg= "#000000", bg = "white", command=getGlycoInfo)
GetGlycoInfo.place(x = 0, y = 90)

GetGlycoSitesNum = tk.Button(root, text = "Get Glycosylated Site Number",
                          padx = 20, pady = 10, fg= "#000000", bg = "white", command=getGlycoSitesNum)
GetGlycoSitesNum.place(x = 250, y = 90)

GetGlycoSitesSequence = tk.Button(root, text = "Get Glycosylated Site Sequence",
                          padx = 20, pady = 10, fg= "#000000", bg = "white", command=getGlcoSequence)
GetGlycoSitesSequence.place(x = 0, y = 140)

blastCluster = tk.Label(root, text = "Blast")
blastCluster.place(x = 0, y = 190)

MakeLocalDatabase = tk.Button(root, text = "Make Local Database",
                          padx = 20, pady = 10, fg= "#000000", bg = "white", command=makeLocalDatabase)
MakeLocalDatabase.place(x = 0, y = 210)

DoLocalBlast = tk.Button(root, text = "Local Blast",
                          padx = 20, pady = 10, fg= "#000000", bg = "white", command=localBlast)
DoLocalBlast.place(x = 250, y = 210)

DefaultLectinNames = tk.Button(root, text = "Generate Default Lectin Names",
                          padx = 20, pady = 10, fg= "#000000", bg = "white", command=defaultLectinName)
DefaultLectinNames.place(x = 0, y = 260)

ExtractEVlaue = tk.Button(root, text = "Extract E-Values",
                          padx = 20, pady = 10, fg= "#000000", bg = "white", command=extractEValue)
ExtractEVlaue.place(x = 250, y = 260)

LectinEValueMap = tk.Button(root, text = "Lectin E-value Map ",
                          padx = 20, pady = 10, fg= "#000000", bg = "white", command=lectinEValueMap)
LectinEValueMap.place(x = 0, y = 310)

CombineMatrix = tk.Button(root, text = "Combine Matrix",
                          padx = 20, pady = 10, fg= "#000000", bg = "white", command=combineMatrix)
CombineMatrix.place(x = 250, y = 310)

MCLResultAnalysis = tk.Label(root, text = "MCL Result Analysis")
MCLResultAnalysis.place(x = 0, y = 360)

MCLResult = tk.Button(root, text = "Process MCL Result",
                      padx = 20, pady = 10, fg= "#000000", bg = "white", command=MCLResultProcess)
MCLResult.place(x = 0, y = 380)

Getclusters = tk.Button(root, text = "Get Clusters",
                        padx = 20, pady = 10, fg= "#000000", bg = "white", command=getclusters)
Getclusters.place(x = 250, y = 380)

tools = tk.Label(root, text = "Tools")
tools.place(x = 0, y = 430)

CombineSequence = tk.Button(root, text = "Combine Sequence",
                          padx = 20, pady = 10, fg= "#000000", bg = "white", command=combineSequence)
CombineSequence.place(x = 0, y = 450)

SplitFiles = tk.Button(root, text = "Split Files ",
                          padx = 20, pady = 10, fg= "#000000", bg = "white", command=splitFiles)
SplitFiles.place(x = 250, y = 450)

Test = tk.Button(root, text = "test ",
                          padx = 20, pady = 10, fg= "#000000", bg = "white", command=test)
Test.place(x = 0, y = 500)

Test2 = tk.Button(root, text = "test2 ",
                          padx = 20, pady = 10, fg= "#000000", bg = "white", command=test2)
Test2.place(x = 250, y = 500)

CountLectins = tk.Button(root, text = "Count Lectins",
                          padx = 20, pady = 10, fg= "#000000", bg = "white", command=countLectins)
CountLectins.place(x = 0, y = 550)

FolderOfSequenceAndSites = tk.Button(root, text = "get sequence and glycosites with given name",
                          padx = 20, pady = 10, fg= "#000000", bg = "white", command=folderOfSequenceAndSites)
FolderOfSequenceAndSites.place(x = 0, y = 600)


CombineWFileName = tk.Button(root, text = "Combine with file name",
                          padx = 20, pady = 10, fg= "#000000", bg = "white", command=combineWFileName)
CombineWFileName.place(x = 0, y = 650)
root.mainloop()