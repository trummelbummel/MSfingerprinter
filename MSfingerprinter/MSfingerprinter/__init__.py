import MSfingerprinter.decoder as decoder
import MSfingerprinter.preprocessing as preprocessing
import MSfingerprinter.pysax as SAX
import MSfingerprinter.periodicityfinder as periodicityfinder
import MSfingerprinter.maxsubpatterntree as maxsubpatterntree
import MSfingerprinter.datacube as datacube
import MSfingerprinter.miningperiodicpatterns as miningperiodicpatterns
import MSfingerprinter.maxsubpatternhitset as maxsubpatternhitset
import MSfingerprinter.reactiontreecompletespace as reactiontreecompletespace
import MSfingerprinter.postprocessing as postprocessing
import MSfingerprinter.comparetrees as comparetrees
import matplotlib.pyplot as plt
import multiprocessing
import os
import pyms
import math
import traceback
import sys
import ntpath
import csv
import itertools
import re
import sys
import pandas as pd
import treelib
import numpy as np
from collections import defaultdict
import gc
import timeit

# np.set_printoptions(threshold=np.inf)

class MSfingerprinter(object):
# MSfingerprinter object instance, constructor method
        def __init__(self):
            super(MSfingerprinter, self).__init__()




############################################################
#
# Postprocessing
#
###############################################################




        def doPostProcessingSingleAlldiffs(self, resultfile1, repeatunitarray):
            cwd = os.getcwd()
            directory = cwd + '/MSfingerprinter/resultsPeriodicity/'
            # outputfiles postprocessing
            outfilemasspatterns = os.path.join(directory, 'postprocessedpatternsMASS' + resultfile1.rstrip('.txt') + 'MASSPATTERNS.csv')
            # outfile massdifferences between and within patterns Mass and Freq
            outfilepatternswithinmass = os.path.join('postprocessed' + resultfile1.rstrip('.txt') + 'DIFFWITHINPATTERNSMASS.txt')
            # save similarity/dissimilarity of freq vs. mass space to file
            periodicmasses, periodicfreqs = postprocessing.openresultssingle(resultfile1)
            # gets from results periodicmasses and periodicfreqs
            postperiodicmasses = postprocessing.retrieveperiodicmasses(periodicmasses)
            return postperiodicmasses



        def searchinTreesgroundtruthstart(self, meaningfuldiffpatterns, path, extensions):
            # for later analysis saves patterns and rootnodes of trees found in initiatortrees
            patternsplusrootnodes = []
            nodesdetected = 0
            treefilearray = list(comparetrees.find_files(path, extensions))
            # for each tree
            for tree in treefilearray:
                nametree = ntpath.basename(tree).rstrip('.json')
                data = comparetrees.getdata(tree)
                nodesfound = []
                counter = 0
                parent = None
                child = None
                reactiontreeinstance = treelib.Tree()
                reactiontreeinstance = comparetrees.retrievenodes(data, counter, parent, child, reactiontreeinstance)
                for i in range(len(meaningfuldiffpatterns)):

                    # get for each meaningfulpattern the root and target
                    rootsubtree = round(meaningfuldiffpatterns[i][3], 6)
                    targetsubtree = round(meaningfuldiffpatterns[i][1],6)
                    pattern = meaningfuldiffpatterns[i][4]
                    try:
                        rootnode, stoichiometformula = comparetrees.searchMasspatterngroundtruth(reactiontreeinstance, rootsubtree, targetsubtree, pattern, nametree)
                        if rootnode != None and stoichiometformula != None:
                            print('original root value without rounding')
                            print(meaningfuldiffpatterns[i][3])
                            print('original target value without rounding')
                            print(meaningfuldiffpatterns[i][1])
                            patternsplusrootnodes.append([rootnode, rootsubtree, targetsubtree, stoichiometformula])
                            nodesfound.append(pattern)
                    except:
                        continue
            print('all trees in directory : ' + path.rstrip('completeInitiatortrees/') + 'resultsubtrees/')
            return patternsplusrootnodes


        def searchinTreesgroundtruth(self, meaningfuldiffpatterns, path, extensions):
            # for later analysis saves patterns and rootnodes of trees found in initiatortrees
            patternsplusrootnodes = []
            nodesdetected = 0
            treefilearray = list(comparetrees.find_files(path, extensions))
            # for each tree
            for tree in treefilearray:
                nametree = ntpath.basename(tree).rstrip('.json')
                data = comparetrees.getdata(tree)
                nodesfound = []
                counter = 0
                parent = None
                child = None
                reactiontreeinstance = treelib.Tree()
                reactiontreeinstance = comparetrees.retrievenodes(data, counter, parent, child, reactiontreeinstance)
                for i in range(len(meaningfuldiffpatterns)):
                    # get for each meaningfulpattern the root and target
                    rootsubtree = round(meaningfuldiffpatterns[i][1], 6)
                    targetsubtree = round(meaningfuldiffpatterns[i][3],6)
                    pattern = meaningfuldiffpatterns[i][4]
                    try:
                        rootnode, stoichiometformula = comparetrees.searchMasspatterngroundtruth(reactiontreeinstance, rootsubtree, targetsubtree, pattern, nametree)
                        if rootnode != None and stoichiometformula != None:
                            print('original root value without rounding')
                            print(meaningfuldiffpatterns[i][1])
                            print('original target value without rounding')
                            print(meaningfuldiffpatterns[i][3])
                            patternsplusrootnodes.append([rootnode, rootsubtree, targetsubtree, stoichiometformula])
                            nodesfound.append(pattern)
                    except:
                        continue
            print('all trees in directory : ' + path.rstrip('completeInitiatortrees/') + 'resultsubtrees/')

            return patternsplusrootnodes




        def constructInitiatorTreescompletespace(self, Initiator, nodelist, massrange):
            counter = 0
            previousname = None
            currenttrees = []
            boolean = True
            initiatorname = Initiator[0]
            maxlevel, maxleveltwo = reactiontreecompletespace.getmaxtreelevel(nodelist, massrange)
            treelen = len(Initiator)
            for i in Initiator:
                trees = reactiontreecompletespace.createRootsInitiators(Initiator)
                counter += 1
                while boolean == True:
                    for i in trees:
                        if len(currenttrees) == treelen:
                            trees = currenttrees
                            currenttrees = []
                            continue
                        else:
                            reactiontreeinstance, counter, previousname = reactiontreecompletespace.createnextlevelinitiator(nodelist, i, counter, maxlevel, maxleveltwo, previousname, initiatorname)
                            if reactiontreeinstance != None and treelen > 0:
                                currenttrees.append(reactiontreeinstance)

                            else:
                                treelen = treelen - 1
                                return

            return



###############################################################################
#
# Preprocessing
#
#################################################################################3

    # preprocessing data for Periodicityfindingalgorithms, returns function (i.e. intensity values)
    # masspoints corresponding m/z values, Mass space

        def preprocessMSSpectraMass(self, filename):
            print('preprocessing raw data of ...'+ filename)
            originaldata, sampleID = decoder.readdataframe(filename)
            masspoints, function = preprocessing.FunctionMSMass(originaldata)
            return function, masspoints

        # preprocessing data for Periodicityfindingalgorithms, returns function (i.e. intensity values)
        # Frequencyspace
        def preprocessMSSpectraFreq(self, filename):
            print('preprocessing raw data of ...'+ filename)
            originaldata, sampleID = decoder.readdataframe(filename)
            freqpoints, function = preprocessing.FunctionMSFreq(originaldata)
            return function, freqpoints

        def preprocess_directoryCSVfreqsampling(self, path, extensions, sampling, nprocesses=None):
                filenames_to_preprocess = []
                colnames = []
                counter = 0
                #  Try to use the maximum amount of processes if not given.
                try:
                    nprocesses = nprocesses or multiprocessing.cpu_count()
                except NotImplementedError:
                    nprocesses = 1
                else:
                    nprocesses = 1 if nprocesses <= 0 else nprocesses
                for filename, _ in decoder.find_files(path, extensions):

                    filenames_to_preprocess.append(filename)
                    try:
                        filename = filename
                    except ValueError:
                        pass
                    # print('preprocessing raw data of ...'+ filename)
                    originaldata, sampleID = decoder.readdataframecsvfreqsampling(filename, sampling)
                    colnames.append(sampleID)
                    # three digits to make join possible
                    # originaldata['freq'] = originaldata['freq'].round(0)
                    if counter == 0:
                        originaldataframe = originaldata
                        counter+= 1
                    else:
                        originaldata = originaldata
                        merged = pd.concat([originaldataframe, originaldata], axis=1)
                        # merged = pd.merge(originaldataframe, originaldata, on='freq', how='outer')
                        originaldataframe = merged
                originaldataframe.fillna(value=np.nan, inplace=True)
                originaldataframe.columns = colnames
                # originaldataframe.sort(columns='freq', inplace=True)
                originaldataframe = originaldataframe.apply(np.log)

                return originaldataframe

        def preprocess_directoryCSVmasssampling(self, path, extensions, sampling, nprocesses=None):
                filenames_to_preprocess = []
                colnames = []

                counter = 0
                #  Try to use the maximum amount of processes if not given.
                try:
                    nprocesses = nprocesses or multiprocessing.cpu_count()
                except NotImplementedError:
                    nprocesses = 1
                else:
                    nprocesses = 1 if nprocesses <= 0 else nprocesses
                for filename, _ in decoder.find_files(path, extensions):

                    filenames_to_preprocess.append(filename)
                    try:
                        filename = filename
                    except ValueError:
                        pass
                    # print('preprocessing raw data of ...'+ filename)
                    originaldata, sampleID = decoder.readdataframecsvmasssampling(filename, sampling)
                    colnames.append(sampleID)
                    # three digits to make join possible
                    # originaldata['mass'] = originaldata['mass'].round(3)
                    if counter == 0:
                        originaldataframe = originaldata
                        counter+= 1
                    else:
                        originaldata = originaldata
                        merged = pd.concat([originaldataframe, originaldata], axis=1)
                        # merged = pd.merge(originaldataframe, originaldata, on='freq', how='outer')
                        originaldataframe = merged
                originaldataframe.fillna(value=999, inplace=True)
                originaldataframe.columns = colnames
                # originaldataframe.sort(columns='freq', inplace=True)
                originaldataframe = originaldataframe.apply(np.log)

                return originaldataframe


        def getmaxminmass(self, path, extensions, nprocesses=None):
                filenames_to_preprocess = []
                minimum = 10000000000000000000000
                maximum = 0
                counter = 0
                #  Try to use the maximum amount of processes if not given.
                try:
                    nprocesses = nprocesses or multiprocessing.cpu_count()
                except NotImplementedError:
                    nprocesses = 1
                else:
                    nprocesses = 1 if nprocesses <= 0 else nprocesses
                for filename, _ in decoder.find_files(path, extensions):
                    filenames_to_preprocess.append(filename)
                    try:
                        filename = filename
                    except ValueError:
                        pass
                    print('preprocessing raw data of ...'+ filename)
                    originaldata, sampleID = decoder.readdataframe(filename)
                    # three digits to make join possible
                    originaldata['mass'] = originaldata['mass'].round(3)
                    maxmass = max(originaldata['mass'])
                    minmass = min(originaldata['mass'])
                    if maxmass > maximum:
                        maximum = maxmass
                    if minmass < minimum:
                        minimum = minmass
                return maximum, minimum


#######################################################3
#
# Entropy based feature ranking
#
###############################################

            # each feature is removed and total entropy is calculated
        def RankFeatures(self, clusteringinputdata, rangenum):
            featurenames = np.arange(len(clusteringinputdata))
            bestfeatures = []

            # produces a matrix with each column being one feature
            X = np.array(clusteringinputdata)
            totallengthfeatures = range(len(featurenames)/rangenum)
            rangestart = 0
            for i in totallengthfeatures:
                rangeend = rangestart+rangenum
                Xslice = X[rangestart:rangeend,0:10]
                featurenamesslice = featurenames[rangestart:rangeend]
                rankedfeatures = clustering.doRankfeatures(Xslice, featurenamesslice, rangestart, rangeend)
                rangestart = rangeend
                bestfeatureofsubset = rankedfeatures[-1]
                bestfeatures.append(bestfeatureofsubset)
            return bestfeatures

        def doStats(self, entropyvectormass, entropyvectorfreq):
            freqstats = clustering.statsentropy(entropyvectorfreq, 'freq')
            massstats = clustering.statsentropy(entropyvectormass, 'mass')
            correlatefreqmassentropy(rankedentropyvectormass, rankedentropyvectorfreq)
            return freqstats, massstats



######################################################
#
# Time series standardization with SAX
#
########################################################



        def standardizeTimeSeries(self, MSarray, timepoints):
            # this uses the Symbolic aggregate approximation for time series discretization
            # https://github.com/dolaameng/pysax/blob/master/Tutorial-SAX%20(Symbolic%20Aggregate%20Approximation).ipynb
            # MS = decoder.readMStimedomaintransient(filename)
            # MSarray = np.asarray(MS.columns.tolist(), dtype=float)

            print('length original data')
            print(len(MSarray))
            # this does symbolization for each window, no overlap of windows
            sax = SAX.SAXModel()
            # standardizes the time series whiten accross windows
            normalizedMSarray = sax.whiten(MSarray)
            print('mean and standarddeviation of standardized MS: ')
            print(normalizedMSarray.mean(), normalizedMSarray.std())
            # saves MSarray (TS) to csv for later cube generation
            dforiginalTS = pd.DataFrame(normalizedMSarray)
            # time points are m/z values of the function
            timepoints = range(len(normalizedMSarray))
            dftime = pd.DataFrame(timepoints)
            frames = [dforiginalTS, dftime]
            concatenated = pd.concat(frames, axis=1)
            concatenated.columns = ['MSvalues', 'timepoints']
            cwd = os.getcwd()
            outpath = cwd + '/MSfingerprinter/MStimeseriesnormalized.csv'
            concatenated.to_csv(outpath, index=False)
            return normalizedMSarray, timepoints



############################################################
#
# Periodicity Hints finder (calls FFT based circularAutocorrelation)
#
###########################################################

        def findPeriodicityHints(self, workingcube, TS, bins, minconfidence):
            periodicityhintarray = []
            binaryvectorarray = datacube.getbinaryvectors(workingcube)
            circularautocorrelations = periodicityfinder.circularAutocorrelation(binaryvectorarray, bins)
            periodicityhintarray = periodicityfinder.findPeriodicHints(TS, circularautocorrelations, bins, minconfidence)
            return periodicityhintarray


##################################################
#
# Datacube Functions
#
################################################3

# construction of a reference cube
        def constructDatacube(self, TS, timepoints, dimensions, masspoints):
            cube, dictcube = datacube.createReferenceCube(TS, timepoints, dimensions, masspoints)

            return cube, dictcube, dictcube['m/z']

        def constructDatacubebothspaces(self, TSmass, TSfreq, timepoints, dimensions, masspoints):
            cube, dictcube = datacube.createReferenceCubebothspaces(TSmass, TSfreq, timepoints, dimensions, masspoints)
            return cube, dictcube, dictcube['m/z'], dictcube['freq']

        def constructbinaryWorkingcube(self, TS, dictcube, cube, numbins):
            # creates binaryTS arrays
            workingcubediscretizedbins, binarydict = datacube.createWorkingCubediscretizationbinning(TS, cube, numbins)
            binarydictcube = datacube.createnestedDictCubebinary(dictcube, binarydict)
            return binarydictcube

        def constructbinaryWorkingcubewithPeriods(self, TS, binarydictcube, periods, period_j):
            # creates periodarrays from binaryTS arrays
            binaryperioddict = datacube.createWorkingCubeperiodBinaryarrays(TS, binarydictcube, periods, period_j)
            binarydictcube = datacube.createnestedDictCubebinaryperiods(binarydictcube, binaryperioddict)
            return binarydictcube

# this constructs the discretized time dependent attribute and will save them as binary arrays
        def constructWorkingcubediscretized(self, TS, cube, dictcube, periods, period_j, mode, numbins):
            mode = mode
            # discretizes working cube
            workingcubediscretizedbins, binarydict = datacube.createWorkingCubediscretizationbinning(TS, cube, numbins)
            # creates a discretized version of the time series mapping discrete bins to symbols
            newcube, symbolicTS = datacube.createWorkingCubediscreteseriesPattern(TS, cube, numbins)
            symbolicperioddict = datacube.createWorkingCubeperiodMSvalues(symbolicTS, newcube, periods, period_j)
            # workingcubediscretizedhistogramequalwidth = datacube.createWorkingCubediscretizationhistogram(TS, workingcubeperiod, mode)
            workingdictcube = datacube.createnestedDictCubesymbolic(dictcube, binarydict, symbolicperioddict)
            return workingdictcube

            # this will add one periodindex to dictcube
        def constructWorkingcubewithperiods(self, TS, cube, dictcube, periods, period_j):
            periodindexdict = defaultdict(dict)
            # constructs period index for given workingcube
            periodindexdict = datacube.createWorkingCubeperiodindex(TS, cube, periods, period_j)
            workingdictcube = datacube.createnestedDictCubeperiod(dictcube, periodindexdict)
            return workingdictcube

        def constructWorkingcubeComplete(self, TS, cube, dictcube, binarycube, periods, period_j, mode, numbins):
            discretizedCube = self.constructWorkingcubediscretized(TS, cube, dictcube,  periods, period_j, mode, numbins)
            periodindexedCube = self.constructWorkingcubewithperiods(TS, cube, dictcube, periods, period_j)
            binaryworkingCube = self.constructbinaryWorkingcubewithPeriods(TS,  dictcube, periods, period_j)
            workingcube = dict(periodindexedCube.items() + discretizedCube.items() + binaryworkingCube.items())

            return workingcube


# T-slice does not contain the original MSvalues anymore, and not the original Timepoints
        def constructTslice(self, workingcube, periods, period_j):
            Tslice = datacube.makeTslice(workingcube, periods, period_j)
            return workingcube

######################################################
#
# Mining periodic patterns from Mining Segment wise periodic patterns in time related databases
#
#####################################################

        def doPeriodicitypatternsearch(self, minconfidencepattern, minconf, periodicityhintarray, normalizedTS, cube, dictcube, binarycube, mode, numbins, space_type, nodelist):
            if len(periodicityhintarray) > 0:
                # regex approach only supports 100 groups maximum
                periodicityhintarray = [i for i in periodicityhintarray if i < 100]
                for i in range(len(periodicityhintarray)):
                    period_j = i
                    print('periodicitysearch for period with length : ' + str(periodicityhintarray[period_j]))
                # workingcube = fingerprinter.constructWorkingcubediscretized(normalizedTS,cube, dictcube,  periodicityhintarray, period_j, mode, numbins)
                    workingcube = self.constructWorkingcubeComplete(normalizedTS, cube, dictcube, binarycube, periodicityhintarray, period_j, mode, numbins)
                    Tslice = self.constructTslice(workingcube, periodicityhintarray, period_j)
                # # apriori mining of F1 for maxsubpatternstree construction
                    F1, F1counts = self.mineF1(Tslice, minconfidencepattern, periodicityhintarray, period_j)

                    if len(F1) != 0:
                        maxsubpatterntree = self.mineMaxSubpatternTree(Tslice, F1, F1counts, periodicityhintarray, period_j, numbins, nodelist)
                        if space_type == 'mass':
                            print('doing massrun')
                            frequentperiods = self.mineFrequentPeriodicitiesMass(Tslice, maxsubpatterntree, minconf, periodicityhintarray, period_j, normalizedTS)

                        elif space_type == 'freq':
                            print('doing FREQ run')
                            frequentperiods = self.mineFrequentPeriodicitiesFreqs(Tslice, maxsubpatterntree, minconf, periodicityhintarray, period_j, normalizedTS)

                    else:
                        print('no 1-cyclepatterns found with these parameters for periodlength ' + str(periodicityhintarray[period_j]) + '. Try to reduct number of bins or reduce minconfidencepattern or minconf thresholds.')
                        continue


############################################################
#
# Mine frequent 1-cyclepatterns from Time series
#
###########################################################

        def mineF1(self, Tslice, minconfidence, periodicityhintarray, period_j):
            # finds onecycles that are corresponding to minconfidence value
            F1, F1counts = miningperiodicpatterns.findOneCyclePatterns(Tslice, minconfidence, periodicityhintarray, period_j) # F1 in paper
            return F1, F1counts


############################################################
#
# Retrieve patterns from pattern Tree
#
###########################################################

        def mineFrequentPeriodicitiesMass(self, Tslice, maxsubpatterntreeinstance, minconfidence, periods, period_j, normalizedTS):
            numperiods = float(len(normalizedTS))/float(periods[period_j])
            frequencies = maxsubpatterntree.getFrequencies(maxsubpatterntreeinstance, minconfidence, numperiods)
            masses = Tslice['m/z']
            if len(frequencies) == 0:
                print('no frequent patterns found for periodlength ' + str(periods[period_j]) + ' try different parametersettings.' )
            else:
                for i in frequencies:
                    # exclude maxpattern Cmax incidences from results
                    if type(i[2]) != dict:
                        print('confidence of pattern is: ' + str(i[0]))
                        print('frequency of pattern is: ' + str(i[1]))
                        print('first period occurence starts at index : ' + str(i[3]))
                        print('occurences at indizes :' + str(i[5]))
                        # retrieve masses where the pattern occurs i.e. where there are non* letters
                        print('masses :' + str(i[4]))
                        print('periodic intenities : ' + str(i[2]))
                        print('periodlength is : ' + str(periods[period_j]))



############################################################
#
# Construct max subpattern tree - only first level because we want to find the maximum subpattern
# can be extended in order to find periodic subpatterns
#
###########################################################


        def mineMaxSubpatternTree(self, Tslice, F1, F1counts, periods, period_j, numbins, nodelist):
            print('constructing Tree...')
            counter = 1
            previoushitset = None
            previousnodes = []
            counterindexpassed = 0
            boolean = True
            dictcounter = 0
            while boolean == True:
                if counter == 1:
                    print('constructing root and first level of Tree ...')
                    Cmax, lengthCmax = maxsubpatternhitset.formCmax(F1, period_j)
                    maxsubpatterntreeinstance = maxsubpatterntree.createRootnode(Cmax, None)
                    nonstar = (len(Cmax))
                    treedepth = maxsubpatterntreeinstance.depth()
                    # TODO change back if problems to getmaxPatterns
                    maxsubpatterntreeinstance, hitsetcomplete = maxsubpatternhitset.getmaxPatternsset(treedepth, Tslice, maxsubpatterntreeinstance,  counterindexpassed, Cmax, periods, period_j, nonstar, nodelist)
                    boolean = False
                    counter += 1
                    previoushitset = hitsetcomplete
            # returns after first iteration uncomment also boolean = False
            return maxsubpatterntreeinstance
