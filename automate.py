"""
Automation of Predicted vs Experimental Comparison
|Requires a specific layout of the filesystem before using
|otherwise it won't work. Main reason for this is because it's easier to work
|with this setup
TODO: Generalise this. Don't need a specific filesystem layout, just the pdb
files and their paths
"""

import os
from math import sqrt
import re
import matplotlib.pyplot as plt

DEST = "/Users/AbdullahAhmad/Desktop/Aspartic_Proteases_Automation"
FOLDERS = os.walk(DEST).next()[1]
OUTPUT = "/Users/AbdullahAhmad/Desktop/Aspartic_Proteases_Automation_Output2"
DIST_CUTOFF = 2
DENS_CUTOFF = 5
PROTEIN_COUNT = len(FOLDERS)
REFERENCE_STRUCT = '4CMS'

#def plot_func(data_list, ref):
#    """
#    THIS IS SO MESSY
#    """
#    combined = []
#    reference = []
#    for main, comp, percentage in zip(*[iter(data_list)]*3):
#        if comp == "Combined":
#            combined.append(percentage)
#        else:
#            reference.append(percentage)
#    ind = np.arange(PROTEIN_COUNT)
#    width = 0.25
#    fig, ax = plt.subplots(figsize=(20,10))
#    rects1 = ax.bar(ind, combined, width, color='r')
#    rects2 = ax.bar(ind+width, reference, width, color='y')
#    ax.set_ylabel('Percentage Correct with ' + str(DENS_CUTOFF) + ' Cut-off')
#    ax.set_title(ref)
#    ax.set_xticks(ind + width / 2)
#    ax.set_xticklabels(FOLDERS)
#    ax.legend((rects1[0], rects2[0]), ('Combined', 'Reference: ' + ref))
#    plt.show()
#    fig.savefig("Score" + ref + ".png")
#    
#        
#def compare(pids, tdir, odir, cutoff, ref):
#    """
#    Main Comparison function
#    STARTING REWRITE MY FRIEND
#    """
#    tot_list = []
#    for pid in pids:
#        current_ref = (tdir + "/" + pid + "/" + pid + "_CW.pdb")
#        if not os.path.exists(odir + "/" + pid):
#            os.mkdir(odir + "/" + pid)
#        for auto in pids:
#            conserved = []
#            predicted = []
#            density = []
#            distance = []
#            generated_file = (tdir + "/" + auto + "/" + "6-Comparison/O.pdb")
#            with open(current_ref, "r") as exp_file:
#                for line in exp_file:
#                    if re.match('(.*)HOH(.*)', line):
#                        els = line.split()
#                        x_1 = float(els[6])
#                        y_1 = float(els[7])
#                        z_1 = float(els[8])
#                        expid = int(els[5])
#                        with open(generated_file, 'r') as pred_file:
#                            for line in pred_file:
#                                if re.match('(.*)HETATM(.*)', line):
#                                    subels = line.split()
#                                    x_2 = float(subels[5])
#                                    y_2 = float(subels[6])
#                                    z_2 = float(subels[7])
#                                    predid = int(subels[1])
#                                    dens = float(subels[8])
#                                    dist = float(sqrt((float(x_1 - x_2) ** 2) + \
#                                    (float(y_1 - y_2) ** 2) + \
#                                    (float(z_1 - z_2) ** 2)))
#                                    if dist <= cutoff:
#                                        conserved.append(expid)
#                                        predicted.append(predid)
#                                        density.append(dens)
#                                        distance.append(dist)
#                                        elist, plist, denlist, dlist = remove_duplicates(conserved, predicted, density,\
#                                                distance, odir, pid, auto)
#            a = (float(len(plist)) / float(max(plist))) * 100
#            if a != 0 and ("Combined" in pred_file.name or ref in pred_file.name): 
#                tot_list.append(pid), tot_list.append(auto), tot_list.append(a)
#    return tot_list

def init_data(input_path):
    """
    Initialise a dictionary of protein IDs and their directories for reading
    |Arguments
    |input_path - a directory containing the relevant proteins. Currently going 
    |            to use the numbered-folder setup.
    """
    #Scan folders in the input directory
    experiment_data = []
    protein_ids = os.walk(input_path).next()[1]
    for prot_id in protein_ids:
        experiment_data.append(prot_id)
        experiment_data.append(input_path + "/" + prot_id + "/" + prot_id + "_CW.pdb")
        experiment_data.append(input_path+ "/" + prot_id + "/" + "6-Comparison/O.pdb")
    return experiment_data
        

def vector_distance_total(input_list, cutoff):
    """
    Calculating distance, using a generator
    |Arguments
    |input_list - an list of floats where every 6 values correspond to 2 vectors
    """
    for x_1, y_1, z_1, x_2, y_2, z_2 in zip(*[iter(input_list)]*6):
        x_m_x, y_m_y, z_m_z = x_1 - x_2, y_1 - y_2, z_1 - z_2
        distance = sqrt(x_m_x ** 2 + y_m_y ** 2 + z_m_z ** 2) 
        if distance <= cutoff:
            yield distance
        
def sortlist(input_list, col, direc):
    """
    Function for sorting rows in a table
    |Arguments
    |input_list - list of lists, make sure the elements of the sub lists are numbers
    |             not strings
    |col - which column to sort in the table
    |direc - ascending, descending
    """
    return sorted(input_list, key=lambda row: (row[col]), reverse=direc)

        
def duplicate_filter(input_list):
    """
    Pure python method to remove duplicates from a table of data.
    |Arguments
    |input_list - a list of lists, with each list being a row in a table
    """
    outputlist = []
    seen = set()
    for row in sortlist(input_list, 3, False):
        if row[0]in seen:
            continue
        seen.add(row[0])
        outputlist.append(row)
    for i in sortlist(outputlist, 1, False):
        yield i

def split_list(l, n):
    """
    Takes a list and breaks it up in to smaller lists
    |Arguments
    |l - A list
    |The number of elements in each sub-list
    """
    for i in range(0, len(l), n):
        yield l[i:i+n]
        
        
def process_data(proteins, output, cutoff):
    """
    Read relevant data from each protein and process. This is the main function.
    |Arguments
    |proteins - a dictionary of the PDB ID and filesystem location
    |output - some path to output all the data to, either as .txt or maybe .csv
    """
    for prot, reference, dummy in zip(*[iter(proteins)]*3):

        if not os.path.exists(output+ "/" + prot):
            os.mkdir(output + "/" + prot)
        for gen_id, dummy, generated in zip(*[iter(proteins)]*3):
            process_list = []
            with open(reference, 'r') as ref_file:
                    for eline in ref_file:
                        if re.match('(.*)HOH(.*)', eline):
                            e_elements = eline.split()
                            x_1 = (float(e_elements[6]))
                            y_1 = (float(e_elements[7]))
                            z_1 = (float(e_elements[8]))
                            expid = int(e_elements[5])
                            with open(generated, 'r') as gen_file:
                                for pline in gen_file:
                                    if re.match('(.*)HETATM(.*)', pline):
                                        p_elements = pline.split()
                                        x_2 = (float(p_elements[5]))
                                        y_2 = (float(p_elements[6]))
                                        z_2 = (float(p_elements[7]))
                                        predid = int(p_elements[1])
                                        dens = float(p_elements[8])
                                        x_m_x, y_m_y, z_m_z = x_1 - x_2, y_1 - y_2, z_1 - z_2
                                        distance = sqrt(x_m_x ** 2 + y_m_y ** 2 + z_m_z ** 2) 
                                        if distance <= cutoff:
                                            process_list.append(expid)
                                            process_list.append(predid)
                                            process_list.append(dens)
                                            process_list.append(distance)
            #----
            for i in duplicate_filter(split_list(process_list, 4)):
                if i[2] > 5:
                    print prot, gen_id, i

                            

e_dat = init_data(DEST)
process_data(e_dat, OUTPUT, DIST_CUTOFF)