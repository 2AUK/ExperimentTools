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
        ofile.write("Crystal ID\tPredicted ID\tDistance\n")
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

                                distance = float(math.sqrt((float(x_1 - x_2) ** 2) + \
                                    (float(y_1 - y_2) ** 2) + \
                                    (float(z_1 - z_2) ** 2)))

                                if distance <= CUTOFF:
                                    ofile.write(Exp_ID + "\t\t")
                                    ofile.write(Pred_ID + "\t\t")
                                    ofile.write(str(distance) + "\n")
ofile.close()
