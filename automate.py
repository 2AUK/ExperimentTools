"""
Automation of Predicted vs Experimental Comparison
|Requires a specific layout of the filesystem before using
|otherwise it won't work. Main reason for this is because it's easier to work 
|with this setup
TODO: Generalise this. Don't need a specific filesystem layout, just the pdb
files and their paths
"""

import os
import math
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

DEST = "/Users/AbdullahAhmad/Desktop/Aspartic_Proteases_Automation"
FOLDERS = os.walk(DEST).next()[1]
OUTPUT = "/Users/AbdullahAhmad/Desktop/Aspartic_Proteases_Automation_Output"
DIST_CUTOFF = 2
DENS_CUTOFF = 5
PROTEIN_COUNT = len(FOLDERS)
REFERENCE_STRUCT = '4CMS'

def remove_duplicates(elist, plist, denlist, dlist, odir, ref_id, struct_id):
    """
    Remove any duplicates from the data
    """
    data = (zip(elist, plist, denlist, dlist))
    duplicate_df = pd.DataFrame(data, columns=["Experimental ID", "Predicted ID",\
    "Density", "Distance"])
    duplicate_df.sort_values("Distance", inplace=True, ascending=False)
    duplicate_df.drop_duplicates(subset=["Experimental ID"], keep="last", inplace=True)
    duplicate_df.sort_values("Predicted ID", inplace=True, ascending=True)
    duplicate_df.to_csv(odir + "/" + ref_id + "/" + "Ref_" + ref_id + "_CF_" +\
            struct_id + ".txt", index=False, sep='\t')
    duplicate_df.query("Density > " + str(DENS_CUTOFF)).to_csv(odir + "/" + ref_id + "/" + "Ref_" + ref_id + "_CF_" +\
            struct_id + "_C_5.txt", index=False, sep='\t')
    return duplicate_df.query("Density > " + str(DENS_CUTOFF))["Experimental ID"].tolist(), \
           duplicate_df.query("Density > " + str(DENS_CUTOFF))["Predicted ID"].tolist(), \
           duplicate_df.query("Density > " + str(DENS_CUTOFF))["Density"].tolist(), \
           duplicate_df.query("Density > " + str(DENS_CUTOFF))["Distance"].tolist()

def plot_func(data_list, ref):
    combined = []
    reference = []
    for main, comp, percentage in zip(*[iter(data_list)]*3):
        if comp == "Combined":
            combined.append(percentage)
        else:
            reference.append(percentage)
    ind = np.arange(PROTEIN_COUNT)
    width = 0.25
    fig, ax = plt.subplots(figsize=(20,10))
    rects1 = ax.bar(ind, combined, width, color='r')
    rects2 = ax.bar(ind+width, reference, width, color='y')
    ax.set_ylabel('Percentage Correct with ' + str(DENS_CUTOFF) + ' Cut-off')
    ax.set_title('Test_Metric')
    ax.set_xticks(ind + width / 2)
    ax.set_xticklabels(FOLDERS)
    ax.legend((rects1[0], rects2[0]), ('Combined', 'Reference: ' + ref))
    plt.show()
    fig.savefig("Score" + ref + ".png")
        
def compare(pids, tdir, odir, cutoff, ref):
    """
    Main Comparison function
    """
    tot_list = []
    for pid in pids:
        current_ref = (tdir + "/" + pid + "/" + pid + "_CW.pdb")
        if not os.path.exists(odir + "/" + pid):
            os.mkdir(odir + "/" + pid)
        for auto in pids:
            conserved = []
            predicted = []
            density = []
            distance = []
            generated_file = (tdir + "/" + auto + "/" + "6-Comparison/O.pdb")
            with open(current_ref, "r") as exp_file:
                for line in exp_file:
                    if re.match('(.*)HOH(.*)', line):
                        els = line.split()
                        x_1 = float(els[6])
                        y_1 = float(els[7])
                        z_1 = float(els[8])
                        expid = int(els[5])
                        with open(generated_file, 'r') as pred_file:
                            for line in pred_file:
                                if re.match('(.*)HETATM(.*)', line):
                                    subels = line.split()
                                    x_2 = float(subels[5])
                                    y_2 = float(subels[6])
                                    z_2 = float(subels[7])
                                    predid = int(subels[1])
                                    dens = float(subels[8])
                                    dist = float(math.sqrt((float(x_1 - x_2) ** 2) + \
                                    (float(y_1 - y_2) ** 2) + \
                                    (float(z_1 - z_2) ** 2)))
                                    if dist <= cutoff:
                                        conserved.append(expid)
                                        predicted.append(predid)
                                        density.append(dens)
                                        distance.append(dist)
                                        elist, plist, denlist, dlist = remove_duplicates(conserved, predicted, density,\
                                                distance, odir, pid, auto)
            a = (float(len(plist)) / float(max(plist))) * 100
            if a != 0 and ("Combined" in pred_file.name or ref in pred_file.name): 
                tot_list.append(pid), tot_list.append(auto), tot_list.append(a)
    return tot_list
                                             
for i in FOLDERS:
    output = compare(FOLDERS, DEST, OUTPUT, DIST_CUTOFF, i)
    plot_func(output, i)