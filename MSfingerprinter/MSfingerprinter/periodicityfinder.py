from collections import Counter
from collections import OrderedDict
import os
import fnmatch
import numpy as np
import csv
import sys
import pandas as pd
import re
from sklearn import preprocessing
from scipy import signal
from pymatbridge import Matlab
import matplotlib.pyplot as plt
from detect_peaks import detect_peaks
import collections
from pymining import seqmining
from scipy.stats import itemfreq
from numpy.fft import fft, ifft

'''
minePeriodicities is finding period lengths for subsequent Hans maxsubpattern hitset algorithm applied for each candidate period
once periodic hints are found we have the
periods we want to use in order to look for periodicities in the time series
parameters a Timeseries, a specified period (low, high) i.e. candidate period length
and integer m = threshold indicating the ratio of the lengths of S and the patterns must be at least
m to be considered a significant periodicity
'''

# creates a binary vector of length symbolicTS for every letter in the alphabet
# if letter is present at position 1 is put, else 0
def createBinaryvectors(symbolicTS, alphabet):
    binaryvectorarray = []
    for i in alphabet:
        indexarray = []
        vector = np.zeros(len(symbolicTS))
        for j in range(len(symbolicTS)):
            if i == symbolicTS[j]:
                indexarray.append(j)
                np.put(vector, indexarray, 1)
        binaryvectorarray.append(vector)
    return binaryvectorarray


# finds circular normalized autocorrelation for periodicityhints
def circularAutocorrelation(binaryvectorarray, alphabet):
    circularAutocorrelationarray = []
    cwd = os.getcwd()
    for i in range(len(binaryvectorarray)):
        circularAutocorrelationarray.append(ifft(fft(binaryvectorarray[i]) * fft(binaryvectorarray[i]).conj()).real)
        # normalize CAF
    circularAutocorrelationarray = (circularAutocorrelationarray - np.mean(circularAutocorrelationarray)) / (np.std(circularAutocorrelationarray) * len(circularAutocorrelationarray))
    return circularAutocorrelationarray


# circular autocorrelation arrays
# first element in autocorrelation array is the dot product of the vector with itself since shifting lag is 0
# returns an array with possible period lengths of periods that are above minconfidence threshold
# the minconfidence threshold is a ratio of the the number of periods observed in the TS
# vs. the maximum possible periods given the period length

def findPeriodicHints(symbolicTS, circularautocorrelations, alphabet, minconfidence):
    periodicityhintsarray = []
    selectedpeaks = []
    indizespeaks = []
    possibleperiods = []
    possibleperiodpeaks = []
    candidates = []
    countedpeaksarray = []
    N = len(symbolicTS)
    # there should be as many circularautocorrelations as bins
    for i in range(len(circularautocorrelations)):
        # remove first element
        array = np.delete(np.array(circularautocorrelations[i]), 0)
        letter = alphabet[i]
        # select peaks that are greater than confidence threshold
        # which are half the size of max peaks, ignore first peak
        confidencethreshold = 0.5 * np.amax(array)
        indizes = detect_peaks(array, mph=confidencethreshold)
        selectedpeaks.append(np.take(array, indizes))
        indizespeaks.append(indizes)
    for j in range(len(selectedpeaks)):
        peaks = selectedpeaks[j]
        indizes = indizespeaks[j]
        periodpeaks, uniqueindizes, counts = np.unique(peaks, return_index=True, return_counts=True)
        countedpeaksarray.append(counts)
        possibleperiodpeaks.append(periodpeaks)
        possibleper = np.take(indizes, uniqueindizes)
        possibleperiods.append(possibleper)
    for k in range(len(possibleperiods)):
        pp = possibleperiods[k]
        cp = countedpeaksarray[k]
        for n in range(len(pp)):
            period = pp[n]
            maxpossibleperiods = float(N)/float(period)
            countedpeaks = cp[n]
            confidence = float(countedpeaks)/float(maxpossibleperiods)
            if minconfidence < confidence and pp[n] <= (len(symbolicTS)/2):
                # if minconfidence < confidence then the period is a candidate
                # print(pp[n]+1)
                candidates.append(pp[n])
    periodicityhintsarray.append(candidates)
    periodicityhintsarray = np.unique(periodicityhintsarray)
    return periodicityhintsarray
