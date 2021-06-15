# Written by Luuk Nolden and Philippe Bors
''' 
This code extracts the PCAWG concensus calls mentioned in Cortes-Ciriano et al. 2020.
Per genome in this set (located in dataset/sv and dataset/cn) it calls the ShatterSeek
tool to detect chromothripsis in the given genome. Its summary is returned to this code
and analysed. The result is put into a dataframe and serves as an input to the
create_neuralnetwork.py file.
'''
import pandas as pd
import numpy as np
import os
import subprocess

# File postfix
cn_postfix = ".consensus.20170119.somatic.cna.txt"
sv_postfix = ".pcawg_consensus_1.6.161116.somatic.sv.bedpe"

# Boolean to clean CN data
clean_cn_data = False

def Retrieve_entry_list(path, extension):
    # Retrieves all file entries from the given path with the given extension
    new_list = []
    
    for file in os.listdir(path):
        if file.endswith(extension):
            new_list.append(file)

    return new_list

def Retrieve_entry_names(_list):
    # Retrieves all names from the given list
    new_list = []
    
    for entry in _list:
        new_list.append(entry.split('.')[0])

    return new_list

# Create a candidate list of genomes
cn_path = "./dataset/cn/"
sv_path = "./dataset/sv/"

cn_entries_full = Retrieve_entry_list(cn_path, '.txt')
sv_entries_full = Retrieve_entry_list(sv_path, '.bedpe')

cn_entries = Retrieve_entry_names(cn_entries_full)
sv_entries = Retrieve_entry_names(sv_entries_full)

genome_list = []

# Put all entries in a list if we have a CN and a SV copy of the data
for entry in cn_entries:
    if entry in sv_entries:
        genome_list.append(entry)

# This piece of codes cleans the cn_data from entries ShatterSeek does not accept
if clean_cn_data:
    for entry in genome_list:
        cn_file = cn_path + entry + cn_postfix

        with open(cn_file, "r+") as f:
            d = f.readlines()
            f.seek(0)
            for i in d:
                i = i.replace("NA", "-1") #TODO: Is -1 the correct approach?
                if i[0] != "Y": # Delete all Y chromosome data (ShatterSeek issue)
                    f.write(i)
            f.truncate()
            f.close()

# Now for every genome we have, call ShatterSeek to create a nice dataset

# Load the dataframe (create if non existant)
df_model = pd.read_csv('./csv/dataframe.csv', sep=',')
genome_done_list = df_model['genome'].tolist()

# Go through all the genomes and add them to the dataframe. Keeps track of genomes already added
for entry in genome_list:

    # If we already done this genome, just skip it
    if entry in genome_done_list:
        continue

    cn_param = cn_path + entry + cn_postfix
    sv_param = sv_path + entry + sv_postfix

    print(cn_param, sv_param)

    # Call the ShatterSeek code to detect chromothripsis
    poll = subprocess.Popen(['/usr/bin/Rscript', './luuk.R', cn_param, sv_param])
    
    # Check if ShatterSeek is finished
    while True:
        if poll.poll() is not None:
            print('done!')
            break
        else:
            pass
    
    # Make results ready for the dataframe
    df = pd.read_csv('./csv/summary.csv', sep=',')
    df_filtered = df[['chrom', 'start', 'end', 'inter_other_chroms']]
    df_filtered['bool'] = 0
    
    # Create the Boolean data for Yes chromothripsis or No chromothripsis 
    df_filtered.loc[df_filtered['start'].isna(), 'bool'] = 0
    df_filtered.loc[df_filtered['start'].notna(), 'bool'] = 1

    df_row = {'genome': entry, '1': df_filtered}
    
    temp_list = []

    # Add all chromosomes to the genome in a single row
    for i in range(len(df_filtered)):
        temp_list.append(df_filtered['bool'][i])

    temp_list.insert(0,entry) # Now we have the row we need to append

    series = pd.Series(temp_list, index = df_model.columns)
    df_model = df_model.append(series, ignore_index=True)

    print('Cycle complete')

    # Save the appended dataframe
    df_model.to_csv('./csv/dataframe.csv', index = False, header=True)

    print(df_model)

exit(0)



