import os
import gc
import subprocess
import re

# file handling
import pandas as pd
import csv
# numerical
import numpy as np
import scipy
from scipy.signal import argrelextrema
from decimal import Decimal

# visualizing
import matplotlib.pyplot as plt


# this creates a representation of the mass spectrum as a numerical series of relevant (non neg. massdefect)
# intensities plotted on the y axis and m/z values plotted on the x axis
# this serves as input for the periodicityalgorithm:
# it also removes neg. massdefect and noise data

def FunctionMSMass(df):
    ''' this excludes noise i.e. if first decimal is between 0 and 100, then only take if first decimal is below 1
        remainder of integer(mass_i) will give the value to check for and the value should be relevant for the specific mass range'''
    relevantmass = []
    relevantintensities = []
    for i in range(len(df)):
        mass_i = df.iloc[i,1]
        if (mass_i >= 0 and mass_i < 100) and (mass_i%1 < 0.05):
            relevantmass.append(df.iloc[i,[1]])
            relevantintensities.append(df.iloc[i,[3]])
        elif (mass_i >= 100 and mass_i < 1000) and (mass_i%1 < (mass_i/1000)):
            relevantmass.append(df.iloc[i,[1]])
            relevantintensities.append(df.iloc[i,[3]])
        elif (mass_i >= 1000 and mass_i < 10000) and (mass_i%1 < (mass_i/10000)):
            relevantmass.append(df.iloc[i,[1]])
            relevantintensities.append(df.iloc[i,[3]])
    return relevantmass, relevantintensities

def FunctionMSFreq(df):
    ''' this excludes noise i.e. if first decimal is between 0 and 100, then only take if first decimal is below 1
        remainder of integer(mass_i) will give the value to check for and the value should be relevant for the specific mass range'''
    relevantfreq = []
    relevantintensities = []
    for i in range(len(df)):
        mass_i = df.iloc[i,1]
        if (mass_i >= 0 and mass_i < 100) and (mass_i%1 < 0.05):
            relevantfreq.append(df.iloc[i,[2]])
            relevantintensities.append(df.iloc[i,[3]])
        elif (mass_i >= 100 and mass_i < 1000) and (mass_i%1 < (mass_i/1000)):
            relevantfreq.append(df.iloc[i,[2]])
            relevantintensities.append(df.iloc[i,[3]])
        elif (mass_i >= 1000 and mass_i < 10000) and (mass_i%1 < (mass_i/10000)):
            relevantfreq.append(df.iloc[i,[2]])
            relevantintensities.append(df.iloc[i,[3]])
    return relevantfreq, relevantintensities


# # to file can be any kind of array dynamical naming
def resultstofile(filename, array):
    "testfunction prints to file so i can continue computing"
    keys = array[0][0].keys()
    with open(filename, 'wb') as output_file:
        dict_writer = csv.DictWriter(output_file, keys)
        dict_writer.writeheader()
        for i in array:
            dict_writer.writerows(i)
