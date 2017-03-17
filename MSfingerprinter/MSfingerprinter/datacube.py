import bitarray
import numpy as np
import pandas as pd
from sqlalchemy import create_engine
from cubes.tutorial.sql import create_table_from_csv
from cubes import Workspace
import subprocess
from collections import defaultdict, OrderedDict
import itertools


# datacube implemented as n dimensional array with periodindex, time index and time related m/z values
def createReferenceCube(TS, timepoints, dimensions, masspoints):
    TS = np.array(TS)
    timepoints = np.array(timepoints, dtype=float)
    masspoints = np.array(masspoints, dtype=float)

    dictionarycube = defaultdict(dict)
    # TODO: figure out a way to not already in the reference cube contain all dimensions
    # add another dimension to the datacube for a period index that has a period index for the specified period
    cube = np.ones((len(TS), dimensions))
    TS = TS.reshape((len(TS),))
    masspoints = masspoints.reshape((len(TS),))
    timepoints = timepoints.reshape((len(TS), ))
    cube[:,0] = TS
    cube[:,1] = timepoints
    cube[:,2] = masspoints

    dictionarycube['MSvalues'] = cube[:,0]
    dictionarycube['timepoints']= cube[:,1]
    dictionarycube['m/z'] = cube[:,2]
    return cube, dictionarycube

def createReferenceCubebothspaces(TSmass, freqpoints, timepoints, dimensions, masspoints):
    TSmass = np.array(TSmass)
    timepoints = np.array(timepoints, dtype=float)
    masspoints = np.array(masspoints, dtype=float)
    freqpoints = np.array(freqpoints, dtype=float)


    dictionarycube = defaultdict(dict)
    # TODO: figure out a way to not already in the reference cube contain all dimensions
    # add another dimension to the datacube for a period index that has a period index for the specified period
    cube = np.ones((len(TSmass), dimensions))
    TSmass = TSmass.reshape((len(TSmass),))

    masspoints = masspoints.reshape((len(TSmass),))
    timepoints = timepoints.reshape((len(TSmass), ))
    freqpoints = freqpoints.reshape((len(TSmass),))
    cube[:,0] = TSmass
    cube[:,1] = timepoints
    cube[:,2] = masspoints
    cube[:,3] = freqpoints

    dictionarycube['MSvalues'] = cube[:,0]
    dictionarycube['timepoints']= cube[:,1]
    dictionarycube['m/z'] = cube[:,2]
    dictionarycube['freq'] = cube[:,3]
    return cube, dictionarycube


# working cube discretizes the time dependent attribute
# and adds periodindizes for the periods found by the periodicityhintalgorithm
# i.e. now the Working cube is indexed with period indizes as well
# user has to specify the period index for which he wants to generate the working cube
# zeros in the working cube indicate an incomplete period i.e. period indizes start with 1
# furthermore depending on the chosen period_j the period index will be established

#TODO: make this work again!!
def createWorkingCubeperiodindex(TS, cube, periods, period_j):
    perioddict = OrderedDict()
    counter = 0
    periodindexarray = np.zeros((len(TS)))
    maxnumperiods = len(TS)/periods[period_j]
    periodarray = cube[:,2]
    valuesrange = range(maxnumperiods)
    for j in valuesrange:
        start = 0 + (j * periods[period_j])
        end = start + periods[period_j]
        subpart = np.array(periodarray[start:end])
        subpart[:] = j
        perioddict[str(j)] = subpart
        periodindexarray = np.insert(periodindexarray, start, subpart)
        counter += len(subpart)
    lastj = j + 1
    if counter < len(TS):
        diff = len(TS) - counter
        padding = np.full((diff), lastj)
        perioddict[str(lastj)] = padding
        counter += diff

    return perioddict

def createWorkingCubeperiodBinaryarrays(symbolicTS, workingcube, periods, period_j):
    binaryperiodsonevalue = OrderedDict()
    binaryperiods = defaultdict(dict)
    counter = 0
    maxnumperiods = (len(symbolicTS))/(periods[period_j])
    values = workingcube['binaryvector']
    valuesrange = range(maxnumperiods)
    for k,v in values.items():

        for j in valuesrange:
            start = 0 + (j * periods[period_j])
            end = start + periods[period_j]
            subpart = np.array(v[start:end])
            binaryperiodsonevalue[str(j)] = subpart
        # periodindexarray = np.insert(periodindexarray, start, subpart)
            counter += len(subpart)
        lastj = j + 1
        if counter < len(symbolicTS):
            diff = len(symbolicTS) - counter
            padding =  np.array(v[start:end])
            binaryperiodsonevalue[str(lastj)] = padding
            counter += diff
        binaryperiods[k] = binaryperiodsonevalue
    return binaryperiods


def createWorkingCubeperiodMSvalues(symbolicTS, cube, periods, period_j):
    symbolicdict = OrderedDict()
    counter = 0
    maxnumperiods = (len(symbolicTS))/(periods[period_j])
    # change here to 3 if not both spaces!
    values = cube[:,4]
    valuesrange = range(maxnumperiods)
    for j in valuesrange:
        start = 0 + (j * periods[period_j])
        end = start + periods[period_j]
        subpart = np.array(values[start:end])
        symbolicdict[str(j)] = subpart
        # periodindexarray = np.insert(periodindexarray, start, subpart)
        counter += len(subpart)
    lastj = j + 1
    if counter < len(symbolicTS):
        diff = len(symbolicTS) - counter
        padding =  np.array(values[start:end])
        symbolicdict[str(lastj)] = padding
        counter += diff
    # if counter < len(symbolicTS):

    return symbolicdict

# creates periodarrays for masses or freqs respectively depending on the space_type
def createWorkingCubeperiodMasses(masses, periods, period_j):
    symbolicdict = OrderedDict()
    counter = 0
    maxnumperiods = (len(masses))/(periods[period_j])
    values = masses
    valuesrange = range(maxnumperiods)
    for j in valuesrange:
        start = 0 + (j * periods[period_j])
        end = start + periods[period_j]
        subpart = np.array(values[start:end], dtype=float)
        symbolicdict[str(j)] = subpart
        # periodindexarray = np.insert(periodindexarray, start, subpart)
        counter += len(subpart)
    lastj = j + 1
    if counter < len(masses):
        diff = len(masses) - counter
        padding =  np.array(values[start:end], dtype=float)
        symbolicdict[str(lastj)] = padding
        counter += diff
    return symbolicdict

def createWorkingCubeperiodFreqs(freqs, periods, period_j):
    symbolicdict = OrderedDict()
    counter = 0
    maxnumperiods = (len(freqs))/(periods[period_j])
    values = freqs
    valuesrange = range(maxnumperiods)
    for j in valuesrange:
        start = 0 + (j * periods[period_j])
        end = start + periods[period_j]
        subpart = np.array(values[start:end], dtype=float)
        symbolicdict[str(j)] = subpart
        # periodindexarray = np.insert(periodindexarray, start, subpart)
        counter += len(subpart)
    lastj = j + 1
    if counter < len(freqs):
        diff = len(freqs) - counter
        padding =  np.array(values[start:end], dtype=float)
        symbolicdict[str(lastj)] = padding
        counter += diff
    return symbolicdict


# creates periodarrays for masses or freqs respectively depending on the space_type
def createWorkingCubeperiodTimepoints(timepoints, periods, period_j):
    symbolicdict = OrderedDict()
    counter = 0
    maxnumperiods = (len(timepoints))/(periods[period_j])
    values = timepoints
    valuesrange = range(maxnumperiods)
    for j in valuesrange:
        start = 0 + (j * periods[period_j])
        end = start + periods[period_j]
        subpart = np.array(values[start:end], dtype=float)
        symbolicdict[str(j)] = subpart
        # periodindexarray = np.insert(periodindexarray, start, subpart)
        counter += len(subpart)
    lastj = j + 1
    if counter < len(timepoints):
        diff = len(timepoints) - counter
        padding =  np.array(values[start:end], dtype=float)
        symbolicdict[str(lastj)] = padding
        counter += diff
    return symbolicdict

# creates a binary vector of lengthTS for every bin
# if value is present in specific bin put 1, else 0
def createBinaryvectors(TS, digitized, bins):
    counter = 0
    binaryvectordict = defaultdict(dict)
    for i in range(1,len(bins)):
        indexarray = []
        vector = np.zeros(len(TS))
        for j in range(len(TS)):
            if i == digitized[j]:
                indexarray.append(j)
                np.put(vector, indexarray, 1)
        if counter == 0:
            rangeforbins = str(i) +' ' + str(min(TS)) + ' - ' + str(bins[i])
            binaryvectordict[rangeforbins]= vector
            counter += 1
        if i == (len(bins)-1):
            rangeforbins = str(i) +' '+ str(bins[i-1]) + ' - ' + str(max(TS))
            binaryvectordict[rangeforbins]= vector
            pass
        else:
            rangeforbins = str(i)+' ' + str(bins[i-1]) + ' - ' + str(bins[i])
            binaryvectordict[rangeforbins]=vector

    return binaryvectordict


def createWorkingCubediscretizationbinning(TS, workingcubeperiodcube, numbins):
    workingcube = defaultdict(dict)
    # digitize returns indices of bins to which each value belongs
    datarange = (min(TS), max(TS))
    discretized = np.histogram(TS, bins=numbins, range=datarange)
    bins = discretized[1]
    digitized = np.digitize(TS, bins)
    digitized = np.reshape(digitized, (len(TS), 1))
    # adds another column to working cube
    workingcube = np.append(workingcubeperiodcube, digitized, axis=1)
    binarydict = createBinaryvectors(TS, digitized, bins)

    return workingcube, binarydict

# creates a symbolic representation based on binning.
def createWorkingCubediscreteseriesPattern(TS, workingcubeperiodcube, numbins):
    symbolicTS = defaultdict(dict)
    # digitize returns indices of bins to which each value belongs
    datarange = (min(TS), max(TS))
    discretized = np.histogram(TS, bins=numbins, range=datarange)
    bins = discretized[1]
    digitized = np.digitize(TS, bins)
    digitized = np.reshape(digitized, (len(TS), 1))
    # adds another column to working cube
    workingcube = np.append(workingcubeperiodcube, digitized, axis=1)
    symbolicTS = digitized
    return workingcube, symbolicTS


# digitized will give me the indizes i.e. into which bin a time point belongs
# from this we can create binary arrays with 1 if present
# i.e. eventually each row in the datacube represents one bin 1 if present, 0 if absent in bin
# the other rows, period and timepoints remain the same
def createWorkingCubediscretizationhistogram(TS, workingcubeperiodcube, mode):
    workingcube = defaultdict(list)
    length = len(TS) + 1
    datarange = (min(TS), max(TS))
    discretized = np.histogram(TS, bins=mode, range=datarange)
    # second return value of np.histogram is binedges
    bins = discretized[1]
    # digitize returns indices of bins to which each value belongs
    digitized = np.digitize(TS, bins)
    discretized = np.reshape(digitized, (len(TS), 1))
    workingcubeequalbin = np.append(workingcubeperiodcube, discretized, axis=1)
    return workingcubeequalbin

def createnestedDictCubeperiod(dictcube, perioddict):
    cubenesteddict = defaultdict(dict)
    cubenesteddict['period_indizes'] = perioddict
    return cubenesteddict

def createnestedDictCubesymbolic(dictcube, binaryvectordict, symbolicTS):
    cubenesteddict = defaultdict(dict)
    cubenesteddict['binaryvector'] = binaryvectordict
    cubenesteddict['symbolicTS'] = symbolicTS
    return cubenesteddict

def createnestedDictCubebinary(dictcube, binaryvectordict):
    dictcube['binaryvector'] = binaryvectordict
    return dictcube

def createnestedDictCubebinaryperiods(dictcube, binaryperiodsdict):
    dictcube['binaryperiods'] = binaryperiodsdict
    return dictcube


def makeTslice(workingcube, periods, period_j):
    # define periodindex all slice
    Tslice = workingcube
    masses = Tslice['m/z']
    freqs = Tslice['freq']
    symbolicTS = Tslice['symbolicTS']
    timepoints = Tslice['timepoints']
    Tslice['m/zPeriods'] = createWorkingCubeperiodMasses(masses, periods, period_j)
    Tslice['freqPeriods'] = createWorkingCubeperiodFreqs(freqs, periods, period_j)
    Tslice['timepointsPeriods'] = createWorkingCubeperiodTimepoints(timepoints, periods, period_j)
    Tslice['periodindexall'] = period_indexall(workingcube)
    # Tslice['periodicDiscreteMS'] = symbolicTS
    return Tslice

def getbinaryvectors(cubenesteddict):
    binaryvectorarray = []
    binarydict = cubenesteddict['binaryvector']
    for key, value in binarydict.iteritems():
        temp = value
        binaryvectorarray.append(temp)
    return binaryvectorarray

# use np.sum to fill periodindex all over all discrete bins
def period_indexall(workingcube):
    periodindexall = defaultdict(dict)
    firstoccurence = []
    values = workingcube['binaryperiods']
    for k, v in workingcube['binaryperiods'].iteritems():

        stackedbinaryperiodarray = np.array(list(v.values()))
        sumstackedarr = np.sum(stackedbinaryperiodarray, axis=0)
        periodindexall[k] = sumstackedarr
    return periodindexall
