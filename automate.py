"""
Automation of Predicted vs Experimental Comparison
"""

import os
import math
import re

ROOT = "/Users/AbdullahAhmad/Desktop/Aspartic_Proteases_Automation"
FOLDERS = os.walk(ROOT).next()[1]
OUTPUT = "/Users/AbdullahAhmad/Desktop/Aspartic_Proteases_Automation_Output"
CUTOFF = 2

for pid in FOLDERS:
    current_ref = (ROOT + "/" + pid + "/" + pid + "_CW.pdb")
    if not os.path.exists(OUTPUT + "/" + pid):
        os.mkdir(OUTPUT + "/" + pid)
    for auto in FOLDERS:
        generated_file = (ROOT + "/" + auto + "/" + "6-Comparison/O.pdb")
        ofile = open(OUTPUT + "/" + pid + "/" +"Ref_" + pid + "_Gen_" + auto + ".txt", 'w')
        ofile.write("Crystal ID\tPredicted ID\tDensity\tDistance\n")
        with open(current_ref, "r") as exp_file:
            for line in exp_file:
                if re.match('(.*)HOH(.*)', line):
                    els = line.split()
                    x_1 = float(els[6])
                    y_1 = float(els[7])
                    z_1 = float(els[8])
                    Exp_ID = els[5]

                    with open(generated_file, 'r') as pred_file:
                        for line in pred_file:
                            if re.match('(.*)HETATM(.*)', line):
                                el = line.split()
                                x_2 = float(el[5])
                                y_2 = float(el[6])
                                z_2 = float(el[7])
                                Pred_ID = el[1]
                                Density = el[8]

                                distance = float(math.sqrt((float(x_1 - x_2) ** 2) + \
                                    (float(y_1 - y_2) ** 2) + \
                                    (float(z_1 - z_2) ** 2)))

                                if distance <= CUTOFF:
                                    ofile.write(Exp_ID + "\t\t")
                                    ofile.write(Pred_ID + "\t\t")
                                    ofile.write(Density + "\t\t")
                                    ofile.write(str(distance) + "\n")

ofile.close()
def remove_duplicates(elist, plist, denlist, dlist):
    """
    Remove any duplicates from the data
    """
    data = (zip(elist, plist, denlist, dlist))
    duplicate_df = pd.DataFrame(data, columns=["Experimental ID", "Predicted ID",\
    "Density", "Distance"])
    duplicate_df.sort_values("Distance", inplace=True, ascending=False)
    duplicate_df.drop_duplicates(subset=["Conserved ID"], keep="last", inplace=True)
    return duplicate_df["Conserved ID"].tolist(), duplicate_df["Reference ID"].tolist(), \
    duplicate_df["Density"].tolist(), duplicate_df["Distance"].tolist()

def write_text(filewrite, elist, plist, denlist, dlist):
    """
    Writing all the text files in a convenient and contained manner.
    """
    for i in zip(elist, plis, denlist, dlist):
        filewrite.write(i[0])
        filewrite.write(i[1])
        filewrite.write(i[2])
        filewrite.write(i[3])
