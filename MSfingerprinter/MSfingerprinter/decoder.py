import os
import fnmatch
import numpy as np
import csv
import sys
import pandas as pd
import re
from sklearn import preprocessing
from scipy import signal
from scipy import stats

def readinputclustering(filename, preprocessingmode):
    df = pd.read_csv(filename, header=None)
    X = df.ix[:, 1::].astype(float)
    X.fillna(0, inplace=True)
    labels = df.ix[:, 0]
    if preprocessing == 'log':
        # log transform the dataframe cause values differ by orders of magnitude
        X = np.log(X)
        X[~np.isfinite(X)] = 0
        labels = df.ix[:, 0]
    else:
        min_max_scaler = preprocessing.MinMaxScaler()
        x_scaled = min_max_scaler.fit_transform(X)
        X = pd.DataFrame(x_scaled)
        labels = df.ix[:, 0]
    return X, labels

# reads input for SAX time series discretization
def readMStimedomaintransient(filename):
    MS = pd.read_csv(filename, sep =',')
    return MS

def crawl_folders(path, extensions):
    directories = []
    for dirpath, dirnames, files in os.walk(path):
        for directory in dirnames:
            if directory != 'rawdata' and directory != 'spectrograms' and directory != 'spectrogrampics' and directory != 'results':
                p = os.path.join(dirpath, directory)
                directories.append(p)
    return directories


# find files path, reads csv files only unless specified differently in extensions
def find_files(path, extensions):
    # Allow both with ".csv" and without "csv" to be used for extensions
    extensions = [e.replace(".", "") for e in extensions]
    for dirpath, dirnames, files in os.walk(path):
        for extension in extensions:
            for f in fnmatch.filter(files, "*.%s" % extension):
                p = os.path.join(dirpath, f)
                yield (p, extension)

# maybe specify a limit parameter such that optionally
# only part of the spectrogram is examined for now leave whole
# spectrogram
# to make comparisons between m/z series normalization within and between samples is necessary
def read(filename):
    spectrogram = pd.read_csv(filename, sep =',')
    min_max_scaler = preprocessing.MinMaxScaler()
    x_scaled = min_max_scaler.fit_transform(spectrogram)
    arr2D = x_scaled
    return arr2D

# read in original freq, mass, intensity data raw data from Basem
def readdataframe(filename):
    sampleID = os.path.basename(filename)
    originaldata = pd.read_table(filename, sep=',', header=0)
    colnames = pd.Series(['Index', 'mass', 'freq', 'intensity'])
    originaldata.columns = colnames
    return originaldata, sampleID


def readdataframecsvrelativefreq(filename, lowerbound, upperbound):
    sampleID = os.path.basename(filename).rstrip('.csv')
    originaldata = pd.read_table(filename, sep=',', header=0)
    # sample ID as placeholder variable for Intensity
    colnames = pd.Series(['Index', 'freq', 'mass', sampleID])
    originaldata.columns = colnames
    mask = originaldata['mass'].between(lowerbound, upperbound, inclusive=True)
    newdata = originaldata.loc[mask]
    del newdata['mass']
    del newdata['Index']
    return newdata, sampleID

def readdataframecsvrelativefreqminmax(filename, sampling):
    sampleID = os.path.basename(filename).rstrip('.csv')
    originaldata = pd.read_table(filename, sep=',', header=0)
    # sample ID as placeholder variable for Intensity
    colnames = pd.Series(['Index', 'freq', 'mass', sampleID])
    originaldata.columns = colnames

    return originaldata, sampleID

# only mass intensity space
def readdataframecsvmasssampling(filename, sampling):
    sampleID = os.path.basename(filename).rstrip('.csv')
    originaldata = pd.read_table(filename, sep=',', header=0)
    # sample ID as placeholder variable for Intensity
    colnames = pd.Series(['Index', 'mass', 'freq', sampleID])
    originaldata.columns = colnames
    del originaldata['freq']
    del originaldata['Index']
    del originaldata[sampleID]
    originaldata.columns = [sampleID]
    dataarray = np.array(originaldata.T)[0]
    originaldataresampled = pd.Series(signal.resample(dataarray, sampling))
    return originaldataresampled, sampleID


# only mass intensity space
def readdataframecsvmass(filename):
    sampleID = os.path.basename(filename).rstrip('.csv')
    originaldata = pd.read_table(filename, sep=',', header=0)
    # sample ID as placeholder variable for Intensity
    colnames = pd.Series(['Index', 'mass', 'freq', sampleID])
    originaldata.columns = colnames
    del originaldata['freq']
    del originaldata['Index']
    del originaldata[sampleID]
    originaldata.columns = [sampleID]
    dataarray = np.array(originaldata.T)[0]
    return originaldata, sampleID


# only mass intensity space
def readdataframecsvmassminmax(filename):
    sampleID = os.path.basename(filename).rstrip('.csv')
    originaldata = pd.read_table(filename, sep=',', header=0)
    # sample ID as placeholder variable for Intensity
    colnames = pd.Series(['Index', 'mass', 'freq', sampleID])
    originaldata.columns = colnames
    del originaldata['freq']
    del originaldata['Index']
    del originaldata[sampleID]
    originaldata.columns = [sampleID]
    return originaldata, sampleID


def readdataframecsvfreqsampling(filename, sampling):
    sampleID = os.path.basename(filename).rstrip('.csv')
    originaldata = pd.read_table(filename, sep=',', header=0)
    # sample ID as placeholder variable for Intensity
    colnames = pd.Series(['Index', 'mass', 'freq', sampleID])
    originaldata.columns = colnames
    del originaldata['mass']
    del originaldata['Index']
    del originaldata[sampleID]
    originaldata.columns = [sampleID]
    dataarray = np.array(originaldata.T)[0]
    originaldataresampled = pd.Series(signal.resample(dataarray, sampling))
    return originaldataresampled, sampleID


def readdataframecsvfreq(filename):
    sampleID = os.path.basename(filename).rstrip('.csv')
    originaldata = pd.read_table(filename, sep=',', header=0)
    # sample ID as placeholder variable for Intensity
    colnames = pd.Series(['Index', 'mass', 'freq', sampleID])
    originaldata.columns = colnames
    del originaldata['mass']
    del originaldata['Index']
    del originaldata[sampleID]
    originaldata.columns = [sampleID]
    dataarray = np.array(originaldata.T)[0]
    return originaldata, sampleID
