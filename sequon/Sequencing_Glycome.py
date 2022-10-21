import tkinter as tk
from tkinter import filedialog
import os
import glob
import shutil
from contextlib import redirect_stdout
import pandas as pd
import numpy as np


def popMessage(message):
    tk.messagebox.showinfo("message",message)
    
    
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
        f = open("./CombinedSequence.fasta", "a")
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
    
    target = "blast_output_test"
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
                                 +targetFolder+currFileNameInTxt +" -outfmt 10")

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
    
root = tk.Tk()
root.title("Sequencing_Glycome")

extractLectin = tk.Button(root, text = "Extract Glycosylated Lectins",
                          padx = 20, pady = 10, fg= "#000000", bg = "white", command=extractGlycoLectins)
extractLectin.pack()

GetGlycoInfo= tk.Button(root, text = "Get Glycosylated Site Inofrmation",
                          padx = 20, pady = 10, fg= "#000000", bg = "white", command=getGlycoInfo)
GetGlycoInfo.pack()

GetGlycoSitesNum = tk.Button(root, text = "Get Glycosylated Site Number",
                          padx = 20, pady = 10, fg= "#000000", bg = "white", command=getGlycoSitesNum)
GetGlycoSitesNum.pack()

GetGlycoSitesSequence = tk.Button(root, text = "Get Glycosylated Site Sequence",
                          padx = 20, pady = 10, fg= "#000000", bg = "white", command=getGlcoSequence)
GetGlycoSitesSequence.pack()

CombineSequence = tk.Button(root, text = "Combine Sequence",
                          padx = 20, pady = 10, fg= "#000000", bg = "white", command=combineSequence)
CombineSequence.pack()

MakeLocalDatabase = tk.Button(root, text = "Make Local Database",
                          padx = 20, pady = 10, fg= "#000000", bg = "white", command=makeLocalDatabase)
MakeLocalDatabase.pack()

DoLocalBlast = tk.Button(root, text = "Local Blast",
                          padx = 20, pady = 10, fg= "#000000", bg = "white", command=localBlast)
DoLocalBlast.pack()

DefaultLectinNames = tk.Button(root, text = "Generate Default Lectin Names",
                          padx = 20, pady = 10, fg= "#000000", bg = "white", command=defaultLectinName)
DefaultLectinNames.pack()

ExtractEVlaue = tk.Button(root, text = "Extract E-Values",
                          padx = 20, pady = 10, fg= "#000000", bg = "white", command=extractEValue)
ExtractEVlaue.pack()

LectinEValueMap = tk.Button(root, text = "Lectin E-value Map ",
                          padx = 20, pady = 10, fg= "#000000", bg = "white", command=lectinEValueMap)
LectinEValueMap .pack()

root.mainloop()