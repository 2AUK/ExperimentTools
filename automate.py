"""
Automation of Predicted vs Experimental Comparison
"""

import os
import math
import re
import pandas as pd
import matplotlib as plt

DEST = "/Users/AbdullahAhmad/Desktop/Aspartic_Proteases_Automation"
FOLDERS = os.walk(DEST).next()[1]
OUTPUT = "/Users/AbdullahAhmad/Desktop/Aspartic_Proteases_Automation_Output"
CUTOFF = 2

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


def compare(pids, tdir, odir, cutoff):
    """
    Main Comparison function
    """
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
                                    x_2 = float(el[5])
                                    y_2 = float(el[6])
                                    z_2 = float(el[7])
                                    predid = int(el[1])
                                    dens = float(el[8])
                                    dist = float(math.sqrt((float(x_1 - x_2) ** 2) + \
                                    (float(y_1 - y_2) ** 2) + \
                                    (float(z_1 - z_2) ** 2)))
                                    if dist <= cutoff:
                                        conserved.append(expid)
                                        predicted.append(predid)
                                        density.append(dens)
                                        distance.append(dist)
                                        remove_duplicates(conserved, predicted, density,\
                                                distance, odir, pid, auto)

compare(FOLDERS, DEST, OUTPUT, CUTOFF)
