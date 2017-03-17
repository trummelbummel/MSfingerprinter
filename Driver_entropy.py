import os
import pandas as pd
import numpy as np
# from MSfingerprinter import *
import MSfingerprinter as MSF
import MSfingerprinter.decoder as decoder
import MSfingerprinter.preprocessing as preprocessing
import MSfingerprinter.clustering as clustering
import MSfingerprinter.pysax as SAX
import MSfingerprinter.periodicityfinder as periodicityfinder
import MSfingerprinter.maxsubpatterntree as maxsubpatterntree
import MSfingerprinter.miningperiodicpatterns as miningperiodicpatterns
import MSfingerprinter.maxsubpatternhitset as maxsubpatternhitset
import MSfingerprinter.postprocessing as postprocessing
import warnings
import json
import subprocess
# import MySQLdb
import sys
import csv
import matplotlib.pyplot as plt

extensionscsv = ['.csv']
# extensionsasc = ['.asc']
cwd = os.getcwd()
samplingrate = 1000.0
# instantiate MSfingerprinter
fingerprinter = MSF.MSfingerprinter()


#################################################
#
# preprocessing
#
#################################################

# # # need to have rawdata in folder rawdata within MSfingerprinter module folder
rawdatadirectory = cwd + '/MSfingerprinter/rawdata/featureranking/'


##################################3
#
# Feature ranking
#
############################################

# gets the minimum number of datapoints in data to resample all other MS spectra to a vector of minimum datapoints
# this creates a common m/z or freq vector i.e. alignment
maximum, minimum = fingerprinter.getmaxminmass(rawdatadirectory, extensionscsv)
print('minimum mass in all dataframes')
print(minimum)
print('maximum mass in all dataframes')
print(maximum)
# resample to maximum number of datapoints in dataframe
sampling = maximum
#sampling = minimum
# regular mass and frequency space clustering
originaldataframecsvmass = fingerprinter.preprocess_directoryCSVmasssampling(rawdatadirectory, extensionscsv, sampling)
originaldataframecsvfreq = fingerprinter.preprocess_directoryCSVfreqsampling(rawdatadirectory, extensionscsv, sampling)

# Entropy based methods to investigate freq vs. mass space
bins = 10
rangenum = 5
# rank features in mass space
rankedfeaturesmass = fingerprinter.RankFeatures(originaldataframecsvmass, rangenum)
listsrankedfeatures = map(list, zip(*rankedfeaturesmass))
featureindizes = listsrankedfeatures[0]
featureentropies = listsrankedfeatures[1]
x = featureindizes
y = featureentropies
plt.title('Entropy of best features in mass space evaluated for each ' + str(rangenum) +' features' )
plt.scatter(x, y)
plt.xlabel('feature indices')
plt.ylabel('entropy')
plt.show()
sortedmassfeatures = sorted(rankedfeaturesmass, key=lambda x: x[1])
print('rankedfeaturesmass mean mass space')
print(np.mean(np.array(y)))
print('rankedfeaturesmass variance mass space')
print(np.var(np.array(y)))

# rank features for frequency space
rankedfeaturesfreq = fingerprinter.RankFeatures(originaldataframecsvfreq, rangenum)
listrankedfeaturesfreq = map(list, zip(*rankedfeaturesfreq))
sortedfreqfeatures = sorted(rankedfeaturesfreq, key=lambda x: x[1])
print('rankedfeaturesfreq mean frequency space')
print(np.mean(np.array(listrankedfeaturesfreq[1])))
print('rankedfeaturesfreq variance frequency space')
print(np.var(np.array(listrankedfeaturesfreq[1])))
# print('ten best in freq space')
# print(sortedfreqfeatures[0:10])

featureindizesfreq = listrankedfeaturesfreq[0]
featureentropiesfreq = listrankedfeaturesfreq[1]
x = featureindizesfreq
y = featureentropiesfreq
plt.title('Entropy of best features in frequency space evaluated for each ' + str(rangenum) +' features' )
plt.scatter(x, y)
plt.xlabel('feature indices')
plt.ylabel('entropy')
plt.show()
