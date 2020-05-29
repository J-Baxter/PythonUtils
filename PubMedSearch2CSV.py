#Python scipt for generating csv and abstract list from PubMed database search
#J Baxter
#28/05/20
#29/05/20 - updated to improve DOI calling and output format

from Bio import Entrez, Medline
import pandas as pd
import numpy as np
import os
from datetime import date
import argparse

#define functions
def removeNA(input):
    for i in range(len(input)):
        b = input[i]
        if b.__class__ == float:
            input[i] = 'Data not returned'
        else:
            b = b
    return input

def findDOI(input):
    output=[]
    for i in range(input.__len__()):
        if "doi" in input[i]:
            b = input[i]
            DOI = b.split(' [')
            output = str(DOI[0])
        else:
            pass

    return output

# create output directory
today = date.today()
d1 = today.strftime("%d%b%y")
propName = "PubMedResults_" + d1

if os.path.isdir(propName) == True:
    newName = "PubMedResults_" + d1 + "_{}"
    counter = 1
    while os.path.isdir(newName.format(counter)):
        counter += 1
    dirName = newName.format(counter)
else:
    dirName = propName

os.mkdir(dirName)
os.chdir(dirName)

# define arguments
parser = argparse.ArgumentParser()
parser.add_argument("query", help="search string for PubMed query")
args = parser.parse_args()

print("\n","The search query you have defined is:", args.query)

print("\n")
email = input("Please enter a valid email address to continue: ")

# Initial Entrez search
Entrez.email = email
handle = Entrez.egquery(term=args.query)
record_0 = Entrez.read(handle)
for row in record_0["eGQueryResult"]:
    if row["DbName"] == "pubmed":
        count = (row["Count"])

print("\n","Your search query returned", count, "results from PubMed.")
print("\n","Now retrieving records from MedLine...")

handle_1 = Entrez.esearch(db="pubmed", term=args.query, retmax=count)
record_1 = Entrez.read(handle_1)
handle_1.close()
idlist = record_1["IdList"]

# Retrieve records from Medline
handle_2 = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="text")
record_2 = Medline.parse(handle_2)

records = list(record_2)
print("Done")

print("\n", len(records), "records returned from Medline.")

# export into pandas dataframe, then export to csv
AID_Array = []
DP_Array = []
TI_Array = []
PMID_Array = []
COM_Array = []
AU_Array = []
AB_Array = []

for record in records:
    PMID_Array.append(record.get("PMID", np.nan))
    AID_Array.append(record.get("AID", np.nan))
    DP_Array.append(record.get("DP", np.nan))
    AU_Array.append(record.get("AU", np.nan))
    TI_Array.append(record.get("TI", np.nan))
    COM_Array.append(record.get("EIN|EFR|CRI|CFR|ECI|ECF|RPI|RPF|RIN|ROF" , 'none'))
    AB_Array.append(record.get("AB", np.nan))

PMID_Array = removeNA(PMID_Array)
AID_Array = removeNA(AID_Array)
DP_Array = removeNA(DP_Array)
AU_Array = removeNA(AU_Array)
TI_Array = removeNA(TI_Array)
COM_Array = removeNA(COM_Array)
AB_Array = removeNA(AB_Array)

YE_Array = []
FAU_Array = []
DOI_Array =[]

#Select YOP only
for i in range(len(DP_Array)):
    DP = DP_Array[i]
    year = int(str(DP)[0:4])
    YE_Array.append(int(year))

#select DOI only
for i in range(len(AID_Array)):
    link = findDOI(AID_Array[i])
    DOI_Array.append(link)

#select first author only
for i in range(len(AU_Array)):
    AU = AU_Array[i]
    FAU = AU[0]
    FAU_Array.append(str(FAU))

print("\n", "Saving results to file...")

df = pd.DataFrame(zip(PMID_Array, TI_Array, FAU_Array,  YE_Array, DOI_Array, COM_Array))

df.columns = ["PMID", "Title", "First Author" , "Publication Year", "DOI", "Observations"]

filename = "queryresults_" + d1 + ".csv"
df.to_csv(filename, index=False)

filename2 = "searchabstracts" + d1 + ".txt"

file = open(filename2 , "w+")
for i in range(DOI_Array.__len__()):
    file.write(TI_Array[i])
    file.write('\n')
    file.write('\n')
    file.write(FAU_Array[i])
    file.write(' ')
    file.write(YE_Array[i])
    file.write('\n')
    file.write('\n')
    file.write(AB_Array[i])
    file.write('\n')
    file.write('\n')
    file.write('PMID: ')
    file.write(PMID_Array[i])
    file.write('\n')
    file.write('DOI: ')
    file.write(DOI_Array[i])
    file.write('\n')
    file.write('\n')
    file.write('\n')
    file.write('-----------------------------')
    file.write('\n')
file.close

print("\n", "Done.")
