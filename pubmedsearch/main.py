##Generating csv and abstract list from PubMed database search through command line interface. For complex search terms,
# remember to use ''.

# import modules

import pandas as pd
import numpy as np
import os
from datetime import date
import argparse
from records import Record

# define functions
def remove_na(input):
    for i in range(len(input)):
        b = input[i]
        if b.__class__ == float:
            input[i] = 'Data not returned'
        else:
            b = b
    return input


def find_doi(input):
    output = []
    for i in range(input.__len__()):
        if "doi" in input[i]:
            b = input[i]
            DOI = b.split(' [')
            output = str(DOI[0])
        else:
            pass

    return output


# create output directory
def get_dirname():
    today = date.today()
    d1 = today.strftime("%d%b%y")
    prop_name = "PubMedResults_" + d1

    if os.path.isdir(prop_name):
        new_name = "PubMedResults_" + d1 + "_{}"
        counter = 1

        while os.path.isdir(new_name.format(counter)):
            counter += 1
        dir_name = new_name.format(counter)

    else:
        dir_name = prop_name

    return dir_name




#refactor records into classes with attributes


def reformat_records(records):
    aid_array = []
    dp_array = []
    ti_array = []
    pmid_array = []
    com_array = []
    au_array = []
    ab_array = []

    for record in records:
        pmid_array.append(record.get("PMID", np.nan))
        aid_array.append(record.get("AID", np.nan))
        dp_array.append(record.get("DP", np.nan))
        au_array.append(record.get("AU", np.nan))
        ti_array.append(record.get("TI", np.nan))
        com_array.append(record.get("EIN|EFR|CRI|CFR|ECI|ECF|RPI|RPF|RIN|ROF", 'none'))
        ab_array.append(record.get("AB", np.nan))

    pmid_array = removeNA(pmid_array)
    aid_array = removeNA(aid_array)
    dp_array = removeNA(dp_array)
    au_array = removeNA(au_array)
    ti_array = removeNA(ti_array)
    com_array = removeNA(com_array)
    ab_array = removeNA(ab_array)

    YE_Array = []
    FAU_Array = []
    DOI_Array = []

    # Select YOP only
    for i in range(len(dp_array)):
        DP = dp_array[i]
        year = int(str(DP)[0:4])
        YE_Array.append(int(year))

    # select DOI only
    for i in range(len(aid_array)):
        link = findDOI(aid_array[i])
        DOI_Array.append(str(link))

    # select first author only
    for i in range(len(au_array)):
        AU = au_array[i]
        FAU = AU[0]
        FAU_Array.append(str(FAU))

    zipped_arrays = zip(pmid_array, ti_array, FAU_Array, YE_Array, DOI_Array, com_array, ab_array)

    return zipped_arrays


def output(zipped_arrays):

    unzipped = [list(t) for t in zip(*zipped_arrays)]
    unzipped.
    print("\n", "Saving results to file...")

    # write csv containing key information (title author Id etc)
    df = pd.DataFrame(zip(pmid_array, ti_array, FAU_Array, YE_Array, DOI_Array, com_array))
    df.columns = ["PMID", "Title", "First Author", "Publication Year", "DOI", "Observations"]

    filename = "queryresults_" + d1 + ".csv"
    df.to_csv(filename, index=False)

    # write abstracts to formatted text file, annotated with title, first author year and ID tags
    filename2 = "searchabstracts" + d1 + ".txt"

    file = open(filename2, "w+")
    for i in range(DOI_Array.__len__()):
        file.write(TI_Array[i])
        file.write('\n')
        file.write('\n')
        file.write(FAU_Array[i])
        file.write(' ')
        file.write(str(YE_Array[i]))
        file.write('\n')
        file.write('\n')
        file.write(AB_Array[i])
        file.write('\n')
        file.write('\n')
        file.write('PMID: ')
        file.write(str(PMID_Array[i]))
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


def main():
    # define arguments for PubMed search from terminal
    parser = argparse.ArgumentParser()
    parser.add_argument("query", help="search string for PubMed query")
    args = parser.parse_args()

    print("\n", "The search query you have defined is:", args.query)
    print("\n")

    email = input("Please enter a valid email address to continue: ")

    # mkdir
    workdir = get_dirname()
    os.mkdir(workdir)
    os.chdir(workdir)


# Initial Entrez search


# Retrieve records from Medline


# export into pandas dataframe, then export to csv



