#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 14:05:50 2017

@author: AbdullahAhmad
"""


import os
import argparse
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
parser.add_argument("path")
args = parser.parse_args()
print args.path

def Total_Prediction_Gen(max_val):
    pred_list = []
    for i in np.arange(1, int(max_val) + 1):
        pred_list.append(float(i))
    return pred_list

def True_Prediction_Gen(prediction_list, true_prediction_list):
    true_pred_list = []
    i = 0.0
    for pred in prediction_list:
        if pred in true_prediction_list:
            i += 1.0
            true_pred_list.append(i)
        else:
            true_pred_list.append(i)
    return true_pred_list

def load_data(path):
    dataframes = []
    for f in path:
        element = []
        all_value_list = []
        pred_list = []
        with open(f, 'r') as input_file:
            input_file.readline()
            name = os.path.splitext(os.path.basename(f))[0]
            name = name.replace("_cvp", "")
            for i, line in enumerate(input_file):
                sub_element = []
                line = line.split()
                pred_list.append(float(line[1]))
                sub_element.extend((line[0], float(line[1]), float(line[2]), float(line[3])))
                all_value_list.append(sub_element)
            data = pd.DataFrame(all_value_list, columns=["Crystal Water ID", \
                                                         "Predicted Water ID", \
                                                         "Solvent Density", \
                                                         "Distance"])
            data.sort_values(["Predicted Water ID"], inplace=True)
            pred_list = sorted(pred_list)
            max_list = 50
            max_prediction = max_list
            tot_plist = Total_Prediction_Gen(max_prediction)
            true_plist = True_Prediction_Gen(tot_plist, pred_list)
            plist = zip(tot_plist, true_plist)
            df = pd.DataFrame(plist, columns=["N", "TP"])
            df = df.assign(TP_N=df.TP / df.N)
            df = df.assign(inverseTP_N=(1 - df.TP_N) / df.TP_N)
            df = df.assign(inverseN=(1 - df.N) / df.N)
            element.append(name)
            element.append(df)
            element.append(data)
            dataframes.append(element)
            return dataframes

def latexify(datalist):
    for data in datalist:
        name = data[0]
        df = data[2]
        print name
        print df.to_latex(index=False)

def visualise(datalist, path):
    for data in datalist:
        name = data[0]
#        df = data[1]
#        ax = df.plot.line(x="N", y="TP_N", legend=False)
#        actual_ax = df.plot(x="N", y="TP_N", title=name, kind='scatter', ax=ax)
#        df.plot(drawstyle="steps")
#        axinverse = df.plot.line(x="inverseN", y="inverseTP_N", legend=False)
#        df.plot(x="inverseN", y="inverseTP_N", title=name + "Test", kind='scatter', ax=axinverse)
#        fig = actual_ax.get_figure()
#        fig.savefig(path + name + ".eps")
        N = data[1]["N"].tolist()
        TP_N = data[1]["TP_N"].tolist()
        plt.figure()
        plt.xlabel("Number of Predictions")
        plt.ylabel("Normalised True Predictions")
        plt.title(name)
        plt.plot(N, TP_N, 'ro')
#        plt.figure()
        plt.step(N, TP_N)

path_list = ["/Users/AbdullahAhmad/Desktop/Parvalbumin_Data/Parvalbumins_Centered/Excel_Data/", \
             "/Users/AbdullahAhmad/Desktop/A_P_Data/", \
             "/Users/AbdullahAhmad/Desktop/FABP_Data/Fatty_Acid_Binding_Proteins_Centered/Excel_Data/", \
             "/Users/AbdullahAhmad/Desktop/Phospho_Datasets_Actual/Phospho_Extended_Centered/Excel_Data/"]
for path in path_list:
    graph_data = load_data(glob.glob(path + "*.txt"))
#    latexify(graph_data)
    visualise(graph_data, path)
