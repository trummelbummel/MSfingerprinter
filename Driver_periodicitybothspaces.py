import os
import sys
import pandas as pd
import numpy as np
import MSfingerprinter as MSF
import MSfingerprinter.decoder as decoder
import MSfingerprinter.preprocessing as preprocessing
import MSfingerprinter.pysax as SAX
import MSfingerprinter.periodicityfinder as periodicityfinder
import MSfingerprinter.maxsubpatterntree as maxsubpatterntree
import MSfingerprinter.miningperiodicpatterns as miningperiodicpatterns
import MSfingerprinter.maxsubpatternhitset as maxsubpatternhitset
import MSfingerprinter.postprocessing as postprocessing
import json


# instantiate MSfingerprinter
fingerprinter = MSF.MSfingerprinter()

'''
uses datacube for more efficient implementation

This main program runs:
1. algorithm to filter periods that are potential periodic patterns
# steps:
 1. scan time series once and create binary vector of size N for every symbol of alphabet in time series
 2. for each symbol of the alphabet compute circular autocorrelation function
 3. scan only half the autocorrelation vector max period is N/2
 minimum confidence threshold (minconfidencecorrelation) for a periodicity to be considered a periodicityhint
Citation: On the Discovery of Weak Periodicities in Large Time Series (Berberidis et al. 2002)

# Hans algorithm for mining periodicities
# input for mining frequent partial periodic patterns in discretized mass spectrum
# 1. Single period apriori to mine frequent one cycles (minconfidencepattern)
Citation: Mining Segment-Wise Periodic Patterns in Time-Related Databases. (Han et al. 1998)
# 2. maximum subpattern hitset algorithm
Citation: Efficient mining of partial periodic patterns in time series database (Han et al. 1999)

IMPORTANT NOTE: This implementation only finds the maximum patterns (Cmax), for frequent subpatterns this implementation
needs to be modified
'''

############################################3
#
# mine periodicities
#
#######################################3

timeingautocorrelation = []
timeingperiodicitymining = []
periodhintsnum = []
datasetsize = []
numbinsarray = []

space_type = 'mass' # 'freq' or 'mass' as alternative spaces as argument
# from to defines the subsection of the spectrum that is used for periodicitysearch
counter = 0
for i in range(0,8):
    if counter == 0:
        from_i = 0
        to_i = 1000
        counter += 1
        datasize = to_i - from_i
        datasetsize.append(datasize)
    else:
        from_i = to_i
        to_i = from_i + 1000
        datasize = to_i - from_i
        datasetsize.append(datasize)

# threshold for autocorrelation to find periodicityhints

    # thresholds for frequent 1-cycles
    print('minconfidence for F1')
    minconfidencepattern = 0.8
    print(minconfidencepattern)
    # if set to 0 patterns will be retrieved even if they occur only once
    # setting threshold low makes sure all patterns based on periodically recurring one cycles are retrieved
    minconf = 0.0
    #  specify dimensions of the referencedatacube
    dimensions = 4
    # define how many bins for the y-axis (intensity) of the series there should be
    print('number of bins')
    numbins = 5
    print(numbins)
    numbinsarray.append(numbins)
    # min confidence value for selection of peaks in autocorrelationfunction
    print('minconfidencecorrelation')
    minconfidencecorrelation = 0.2
    print(minconfidencecorrelation)

    # input file for periodicity mining
    filename = os.getcwd() + '/MSfingerprinter/rawdata/6digitsprecision.csv'

    # logs both to terminal and to outputfile LOGGER
    filenameresults = 'resultsMSarray' + str(from_i) + str(to_i) + str(numbins) + 'bins' + str(minconfidencecorrelation) + 'percentcorrconf' + str(minconfidencepattern) +'percentconf' + space_type + 'bothspaces.txt'
    class Logger(object):
        def __init__(self, filename="Default.log"):
            self.terminal = sys.stdout
            self.log = open(filenameresults, "w")

        def write(self, message):
            self.terminal.write(message)
            self.log.write(message)

    sys.stdout = Logger(filenameresults)
    m = 0.8

    # now we combine freq and mass space patterns
    # FREQ Space, MSarray contains Intensity values
    MSarrayfreq, freqpoints = fingerprinter.preprocessMSSpectraFreq(filename)
    MSarrayfreq = MSarrayfreq[from_i:to_i]
    freqpoints = freqpoints[from_i:to_i]
    # m/z space, MSarray contains Intensity values
    MSarraymass, masspoints = fingerprinter.preprocessMSSpectraMass(filename)
    MSarraymass = MSarraymass[from_i:to_i]
    masspoints = masspoints[from_i:to_i]
    datasize = np.array(MSarraymass)
    masspoints = np.array(masspoints)
    spacepoints = masspoints
    normalizedTS, timepoints = fingerprinter.standardizeTimeSeries(MSarraymass, masspoints)
    cube, dictcube, dictcubemasspoints,dictcubefreqpoints = fingerprinter.constructDatacubebothspaces(normalizedTS, freqpoints, timepoints, dimensions, spacepoints)
    mode = 'auto'
    binarycube = fingerprinter.constructbinaryWorkingcube(normalizedTS, dictcube, cube, numbins)

    # possible reactions and masses NOM i.e. repeat units
    H = 1.007825
    C = 12.0
    O = 15.994915


    H2 = 2.0*H
    addH2 = H2
    subH2 = -(H2)

    # hydration
    H2O = 2.0*H + O
    addH2O = H2O
    subH2O = -(H2O)
    #
    CH2 = 2.0* H + C
    addCH2 = CH2
    subCH2 = -(CH2)

    CO2 = 2*O + C
    addCO2 = CO2
    subCO2 = -(CO2)

    addO = O
    subO = -(O)

    CO = C+O
    addCO = CO
    subCO = -(CO)

    # list of repeating units that are chemically relevant in CHO space
    nodelist = np.array([ H2O,   CH2,  H2,  CO,  CO2, O])

    # numbins + 1 because automatic discretization will add another value which is very sparse
    periodicityhintarray = fingerprinter.findPeriodicityHints(binarycube, normalizedTS, range(numbins + 1), minconfidencecorrelation)
    periodicityhintarray = list(periodicityhintarray)
    periodhintsnum.append(len(periodicityhintarray))
    print('\n')
    print('periodicityhints with periods of length: ' + str(periodicityhintarray) + ' found.')
    fingerprinter.doPeriodicitypatternsearch(minconfidencepattern, minconf, periodicityhintarray, normalizedTS, cube, dictcube, binarycube, mode, numbins, space_type, nodelist)
