import MSfingerprinter.decoder as decoder
import MSfingerprinter.fingerprint as fingerprint
import MSfingerprinter.preprocessing as preprocessing
import MSfingerprinter.constructspectrograms as specmaker
import MSfingerprinter.clustering as clustering
import MSfingerprinter.pysax as SAX
import MSfingerprinter.periodicityfinder as periodicityfinder
import MSfingerprinter.maxsubpatterntree as maxsubpatterntree
import MSfingerprinter.apriori as apriori
import MSfingerprinter.datacube as datacube
import MSfingerprinter.miningperiodicpatterns as miningperiodicpatterns
import MSfingerprinter.maxsubpatternhitset as maxsubpatternhitset
import MSfingerprinter.reactiontree as reactiontree
import MSfingerprinter.reactiontreecompletespace as reactiontreecompletespace
import MSfingerprinter.postprocessing as postprocessing
import MSfingerprinter.comparetrees as comparetrees
import MSfingerprinter.findmasspatternssimple as findmasspatternssimple
import matplotlib.pyplot as plt
import multiprocessing
import os
import icoshift
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


        def doPostProcessing(self, resultfile1, resultfile2):
            cwd = os.getcwd()
            directory = cwd + '/MSfingerprinter/resultsPeriodicity/'

            # outputfiles postprocessing
            outfilefreqpatterns = os.path.join(directory,'postprocessedpatternsFREQ' + resultfile1.rstrip('.txt') + 'FREQPATTERNS.csv')
            outfilemasspatterns = os.path.join(directory, 'postprocessedpatternsMASS' + resultfile2.rstrip('.txt') + 'MASSPATTERNS.csv')

            # outfile massdifferences between and within patterns Mass and Freq
            outfilepatternsbetweenmass = os.path.join(directory,'postprocessed' + resultfile1.rstrip('.txt') + 'DIFFBETWEENPATTERNSMASS.txt')
            outfilepatternswithinmass = os.path.join('postprocessed' + resultfile1.rstrip('.txt') + 'DIFFWITHINPATTERNSMASS.txt')

            outfilepatternsbetweenfreq = os.path.join(directory,'postprocessed' + resultfile1.rstrip('.txt') + 'DIFFBETWEENPATTERNSFREQ.txt')
            outfilepatternswithinfreq = os.path.join(directory,'postprocessed' + resultfile1.rstrip('.txt') + 'DIFFWITHINPATTERNSFREQ.txt')

            # save similarity/dissimilarity of freq vs. mass space to file
            outfilesim = os.path.join(directory, 'postprocessed' + resultfile1.rstrip('.txt') + 'SIMILARITY.csv')
            # print(outfilesim)
            periodicfreqs, periodicmasses, similar, dissimilar = postprocessing.openresults(resultfile1, resultfile2, outfilesim)

            # gets from results periodicmasses and periodicfreqs
            postperiodicmasses, postperiodicfreqs = postprocessing.retrieveperiodicfreqsandmasses(periodicfreqs, periodicmasses)
            # returns massdifferencebetweenpatterns, massdifferencebetweenpatterns
            withinpatterndiffMass, betweenpatterndiffsMass = postprocessing.computeMassDiffsMASS(postperiodicmasses)

            return withinpatterndiffMass


        def doPostProcessingSingle(self, resultfile1, repeatunitarray):
            cwd = os.getcwd()
            directory = cwd + '/MSfingerprinter/resultsPeriodicity/'

            # outputfiles postprocessing
            outfilemasspatterns = os.path.join(directory, 'postprocessedpatternsMASS' + resultfile1.rstrip('.txt') + 'MASSPATTERNS.csv')

            # outfile massdifferences between and within patterns Mass and Freq
            outfilepatternswithinmass = os.path.join('postprocessed' + resultfile1.rstrip('.txt') + 'DIFFWITHINPATTERNSMASS.txt')

            # save similarity/dissimilarity of freq vs. mass space to file

            periodicmasses, periodicfreqs = postprocessing.openresultssingle(resultfile1)
            # print(periodicfreqs)
            # gets from results periodicmasses and periodicfreqs
            postperiodicmasses = postprocessing.retrieveperiodicmasses(periodicmasses)
            postperiodicfreqs = postprocessing.retrieveperiodicfreqs(periodicfreqs)

            # returns massdifferencebetweenpatterns, massdifferencebetweenpatterns
            withinpatterndiffMass = postprocessing.computeMassDiffsMASS(postperiodicmasses)
            withinpatterndiffFreq = postprocessing.computeMassDiffsFREQ(postperiodicfreqs)
            # withinpatterndiffMass = postprocessing.removefrequentitems(withinpatterndiffMass)
            meaningfuldiffpatternsmass, meaningfuldiffpatternsfreq, otherdiffpatterns, otherdiffpatternsfreq = postprocessing.filterMassMeaningfulMassdifferences(withinpatterndiffMass, withinpatterndiffFreq, repeatunitarray)

            meaningfulcount, notmeaningfulcount = postprocessing.countmeaninfuldiffpatternsmass(meaningfuldiffpatternsmass, otherdiffpatterns)
            postprocessing.tofile(resultfile1, meaningfuldiffpatternsmass, meaningfuldiffpatternsfreq, otherdiffpatterns, otherdiffpatternsfreq, meaningfulcount, notmeaningfulcount)
            return meaningfuldiffpatternsmass

        def searchinTrees(self, meaningfuldiffpatterns, path, extensions):
            treefilearray = list(comparetrees.find_files(path, extensions))
            for tree in treefilearray:
                nametree = ntpath.basename(tree).rstrip('.json')
                data = comparetrees.getdata(tree)
                counter = 0
                parent = None
                child = None
                reactiontreeinstance = treelib.Tree()
                reactiontreeinstance = comparetrees.retrievenodes(data, counter, parent, child, reactiontreeinstance)
                for i in range(len(meaningfuldiffpatterns)):
                    rootsubtree = round(meaningfuldiffpatterns[i][2], 4)
                    targetsubtree = round(meaningfuldiffpatterns[i][3],4)
                    pattern = meaningfuldiffpatterns[i][1]

                    try:

                        comparetrees.searchMasspattern(reactiontreeinstance, rootsubtree, targetsubtree, pattern, nametree)
                    except:
                        continue
            print('all trees in directory : ' + path.rstrip('completeInitiatortrees/') + 'resultsubtrees/')


        def searchinTreessimple(self, meaningfuldiffpatterns, path, extensions):
            treefilearray = list(comparetrees.find_files(path, extensions))
            counter = 0
            for tree in treefilearray:
                nametree = ntpath.basename(tree).rstrip('.json')
                data = comparetrees.getdata(tree)
                counter = 0
                parent = None
                child = None
                reactiontreeinstance = treelib.Tree()
                reactiontreeinstance = comparetrees.retrievenodes(data, counter, parent, child, reactiontreeinstance)
                for i in range(len(meaningfuldiffpatterns)):
                    counter += 1
                    rootsubtree = round(meaningfuldiffpatterns[i][0], 3)
                    targetsubtree = round(meaningfuldiffpatterns[i][2],3)

                    pattern = meaningfuldiffpatterns[i]

                    try:
                        comparetrees.searchMasspattern(reactiontreeinstance, rootsubtree, targetsubtree, pattern, nametree)
                    except:
                        continue
            print('all trees in directory : ' + path.rstrip('completeInitiatortrees/') + 'resultsubtrees/')

############################################################
#
#  reaction tree
#
###############################################################

        def constructInitiatorTrees(self, Initiator, nodelist, massrange):
            counter = 0
            previousname = None
            currenttrees = []
            boolean = True
            initiatorname = Initiator[0]
            maxlevel, maxleveltwo = reactiontree.getmaxtreelevel(nodelist, massrange)
            # maxlevel = 4
            treelen = len(Initiator)
            for i in Initiator:
                trees = reactiontree.createRootsInitiators(Initiator)
                counter += 1
                while boolean == True:
                    for i in trees:
                        if len(currenttrees) == treelen:
                            trees = currenttrees
                            currenttrees = []
                            continue
                        else:
                            reactiontreeinstance, counter, previousname = reactiontree.createnextlevelinitiator(nodelist, i, counter, maxlevel, previousname, initiatorname)
                            reactiontreeinstance.show()
                            if reactiontreeinstance != None and treelen > 0:
                                currenttrees.append(reactiontreeinstance)
                                reactiontreeinstance.show()
                            else:
                                treelen = treelen - 1
                                return

            return


        def constructInitiatorTreescompletespace(self, Initiator, nodelist, massrange):

            counter = 0
            previousname = None
            currenttrees = []
            boolean = True
            initiatorname = Initiator[0]
            maxlevel, maxleveltwo = reactiontreecompletespace.getmaxtreelevel(nodelist, massrange)
            # maxlevel = 4
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
                                reactiontreeinstance.show()
                            else:
                                treelen = treelen - 1
                                return

            return



        def constructTreesMassdiff(self, massdifferencewithinpatterns, rootlist, nodelist):
            reactiontreeshit = []
            treecollection = []
            trees = reactiontree.createRootsmassdiff(rootlist)
            counter = 0
            # for each value in a single massdifferencepattern
            # construct the trees
            for i in massdifferencewithinpatterns:
                # stop criterion for tree construction is when currentdifference equals
                # the sum in tree
                currentdifference = i
                # first level of tree
                if counter == 0:
                    for i in trees:
                        parent = 'root'
                        reactiontreeinstance = reactiontree.createnextlevelmassdiff(currentdifference, nodelist, i, parent, counter)
                        counter = reactiontreeinstance[1]
                        # list of trees for next iteration
                        tree = reactiontreeinstance[0]
                        treecollection.append(tree)
                        if type(reactiontreeinstance) != tuple:
                            reactiontreeshit.append((currentdifference, reactiontreeinstance))
                            return reactiontreeshit
                        if len(treecollection)==len(trees):
                            trees = treecollection
                            print('trees counter 0')
                            print(trees)
                            for i in trees:

                                parent = None
                                # feed trees single to create next level
                                reactiontreeinstance = reactiontree.createnextlevelmassdiff(currentdifference, nodelist, i, parent, counter)
                                counter = reactiontreeinstance[1]
                                treecollection.append(reactiontreeinstance[0])
                                # if this returns not a list its a feasible reactiontree
                                if type(reactiontreeinstance) != tuple:
                                    reactiontreeshit.append((currentdifference, reactiontreeinstance))
                                    return reactiontreeshit
                                if len(treecollection)==len(trees):
                                    trees = treecollection
                                    treecollection = []
                                    print('trees not counter 0')
                                    print(trees)
                                    continue
                return reactiontreeshit








######################################################
#
# Shazam fingerprinting function
#
###########################################################


        def fingerprint_directory(self, path, extensions, nprocesses=None):
            #  Try to use the maximum amount of processes if not given.
            filenames_to_fingerprint = []

            try:
                nprocesses = nprocesses or multiprocessing.cpu_count()
            except NotImplementedError:
                nprocesses = 1
            else:
                nprocesses = 1 if nprocesses <= 0 else nprocesses
            # filenames_to_fingerprint = []
            for filename, _ in decoder.find_files(path, extensions):

                 # don't refingerprint already fingerprinted files
                filenames_to_fingerprint.append(filename)
                try:
                    filename = filename

                except ValueError:
                    pass
                print('globally fingerprinting massspectrogram ...'+ filename)
                # Massspectrogramname, extension = os.path.splitext(os.path.basename(filename))
                # Massspectrogram_name = Massspectrogram_name or Massspectrogramname
                arr2D = decoder.read(filename)
                hashfreqtimepointlist = []
                hashfreqtimepoints = fingerprint.storefingerprintfreqs(arr2D, filename)
                # so now i get a list of lists with 2 lists (1 for each channel of freq1, freq2, time lists)
                # print(hashfreqtimepoints)
                #  this should successively be written to a csv file for later use including the Massspectrogram_name
                flatlist = list(itertools.chain(*hashfreqtimepoints))
                print('number of peaks included in fingerprint: ' + str(len(flatlist)))
                if len(flatlist) < 30:
                    flatlist.extend([0] * (30 - len(flatlist)))
                # print(flatlist)
                #  open file in append mode
                # outpath = 'home/trummelbummel/week11/resultfolder/'
                # outfile = os.path.join(outpath, 'globalfingerprints.csv')
                cwd = os.getcwd()
                directory = cwd + '/MSfingerprinter/results/'
                outpath = os.path.join(directory, 'globalfingerprintsnew.csv')
                name = ntpath.basename(filename).rstrip('.csv')
                samplename = re.search(r'[a-zA-Z]+\d+|[a-zA-Z]+|\-[a-zA-Z]+\_[A-Za-z]+', name)
                flatlist.insert(0, samplename.group(0))
                with open(outpath, "a") as output_file:
                    writer = csv.writer(output_file)
                    writer.writerow(flatlist)
                    output_file.close()
            return

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


###############################################################################
#
# Preprocessing
#
#################################################################################3


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


        def preprocessMSSpectracsv(self, filename):
            print('preprocessing raw data of ...'+ filename)
            originaldata, sampleID = decoder.readdataframecsv(filename)
            function = originaldata['intensity'].tolist()
            masspoints = originaldata['mass'].tolist()
            return function, masspoints


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
                    maxmass = len(originaldata)
                    if maxmass > maximum:
                        maximum = maxmass
                    if maxmass < minimum:
                        minimum = maxmass

                return maximum, minimum

    # preprocessing directory for clustering project

        def preprocess_directory(self, path, extensions, nprocesses=None):
                filenames_to_preprocess = []
                #  Try to use the maximum amount of processes if not given.
                try:
                    nprocesses = nprocesses or multiprocessing.cpu_count()
                except NotImplementedError:
                    nprocesses = 1
                else:
                    nprocesses = 1 if nprocesses <= 0 else nprocesses
                for filename, _ in decoder.find_files(path, extensions):
                     # don't refingerprint already fingerprinted files
                    filenames_to_preprocess.append(filename)
                    try:
                        filename = filename
                    except ValueError:
                        pass
                    # print('preprocessing raw data of ...'+ filename)
                    originaldata, sampleID = decoder.readdataframe(filename)
                    # excludes negative mass defect peaks
                    noiselessdata = preprocessing.excludenoise(originaldata)
                    # split up data into nominal masses such that sumsines get created within one nominal mass
                    seperatednominalmassdata = preprocessing.nominalmasses(pd.DataFrame(noiselessdata))
                    # normalizing the mass is still important as it reduces computational complexity
                    # if one had to use the original frequencies sampling rate to capture one period would
                    # imply huge arrays
                    normalizedmassdata = preprocessing.normalizefreq(seperatednominalmassdata)

                    # TODO: this currently overwrites rawdatafiles
                    cwd = os.getcwd()
                    directory = cwd + '/relativefreqresults/'
                    # outputfiles postprocessing
                    outfile = os.path.join(directory,'relativefreqs' + sampleID +'.csv')
                    # preprocessing.resultstofile(filename, normalizedmassdata)
                    # create and save the sumsine files
                    # preprocessing.createsines(normalizedmassdata[0], sampleID)
                return normalizedmassdata


        def preprocess_directoryASC(self, path, extensions, nprocesses=None):
                filenames_to_preprocess = []
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
                    originaldata, sampleID = decoder.readdataframeasc(filename)
                    # three digits to make join possible
                    # originaldata['mass'] = originaldata['mass'].round(3)
                    if counter == 0:
                        originaldataframe = originaldata
                        counter+= 1
                    else:
                        originaldata = originaldata
                        merged = pd.concat([originaldataframe, originaldata], axis=1)
                        originaldataframe = merged
                originaldataframe = icoshift(xt='max', xp=oiginaldataframe)
                # originaldataframe.fillna(value=999, inplace=True)
                # originaldataframe.sort(columns='mass', inplace=True)
                return originaldataframe

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


        def preprocess_directoryCSVfreq(self, path, extensions,  nprocesses=None):
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
                    originaldata, sampleID = decoder.readdataframecsvfreq(filename)
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
                # originaldataframe.columns = colnames
                # originaldataframe.sort(columns='freq', inplace=True)
                # originaldataframe = originaldataframe.apply(np.log)


                originaldataframe = originaldataframe.transpose()

                # originaldataframe = map(np.array, originaldataframe.values)
                originaldataframe = originaldataframe.apply(lambda x: x.values, axis=1)
                originaldataframe = originaldataframe.as_matrix()
                mask = np.isnan(originaldataframe)
                originaldataframe[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), originaldataframe[~mask])

                return originaldataframe, colnames

        def alignMS(self, X):
            X = decoder.alignMassspectra(X)
            return X

        def preprocess_directoryCSVmass(self, path, extensions,  nprocesses=None):
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
                    originaldata, sampleID = decoder.readdataframecsvmass(filename)
                    colnames.append(sampleID)
                    # three digits to make join possible
                    # originaldata['mass'] = originaldata['mass'].round(3)
                    if counter == 0:
                        originaldataframe = originaldata
                        counter+= 1
                    else:
                        originaldata = originaldata
                        merged = pd.concat([originaldataframe, originaldata], axis=1)
                        # merged = pd.merge(originaldataframe, originaldata, on='mass', how='outer')
                        originaldataframe = merged
                print(originaldataframe)
                originaldataframe.columns = colnames
                originaldataframe = originaldataframe.transpose()

                originaldataframe.fillna(value=np.nan, inplace=True)

                # originaldataframe.sort(columns='mass', inplace=True)
                # originaldataframe = originaldataframe.apply(np.log)
                # originaldataframe = originaldataframe.apply(lambda x: x.values, axis=1)
                originaldataframe = originaldataframe.as_matrix()
                # interpolate closest non nan value
                mask = np.isnan(originaldataframe)
                originaldataframe[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), originaldataframe[~mask])
                return originaldataframe, colnames




        def preprocess_directoryCSVfreqminmax(self, path, extensions, lowerbound, upperbound, nprocesses=None):
                filenames_to_preprocess = []
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
                    originaldata, sampleID = decoder.readdataframecsvfreqminmax(filename, lowerbound, upperbound)
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
                originaldataframe.fillna(value=999, inplace=True)
                # originaldataframe.sort(columns='freq', inplace=True)
                # print('freqdataframe')
                # print(originaldataframe)
                originaldataframe = originaldataframe.apply(np.log)

                return originaldataframe

        def preprocess_directoryCSVrelativefreq(self, path, extensions, lowerbound, upperbound, nprocesses=None):
                filenames_to_preprocess = []
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
                    originaldata, sampleID = decoder.readdataframecsvrelativefreq(filename, lowerbound, upperbound)
                    # three digits to make join possible
                    # drop the 20 Hz frequency cause they are masses with little information i.e. only one peak before relative freq calculated
                    # originaldata['freq'] = originaldata['freq'].round(0)
                #     if counter == 0:
                #         originaldataframe = originaldata
                #         originaldataframe = originaldataframe[originaldataframe.freq != 20.0]
                #         counter+= 1
                #     else:
                #         originaldata = originaldata
                #
                #         # inner join intersection not union like in others
                #         merged = pd.merge(originaldataframe, originaldata, on='freq', how='outer')
                #         originaldataframe = merged
                #         originaldataframe = originaldataframe[originaldataframe.freq != 20.0]
                # # drop rows that do not have at least 7 non nan values
                # # originaldataframe.astype(bool).sum(axis=0)
                # print(len(originaldataframe))
                # # originaldataframe = originaldataframe.dropna(axis=0, thresh=1)
                # print('length dataframe after dropping')
                # print(len(originaldataframe))
                # originaldataframe.fillna(value=999, inplace=True)
                # # originaldataframe.drop_duplicates(subset='freq',inplace=True)
                # originaldataframe.sort(columns='freq', inplace=True)
                # originaldataframe.to_csv('/home/trummelbummel/Desktop/resultingrelfreqs.csv', sep=',')
                # originaldataframe = originaldataframe.apply(np.log)
                return originaldata


        def preprocess_directoryCSVmassminmax(self, path, extensions, lowerbound, upperbound, nprocesses=None):
                filenames_to_preprocess = []
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
                    originaldata, sampleID = decoder.readdataframecsvmassminmax(filename, lowerbound, upperbound)
                    # three digits to make join possible

                    originaldata['mass'] = originaldata['mass'].round(3)
                    if counter == 0:
                        originaldataframe = originaldata
                        counter+= 1
                    else:
                        originaldata = originaldata
                        merged = pd.concat([originaldataframe, originaldata], axis=1)
                        # merged = pd.merge(originaldataframe, originaldata, on='mass', how='outer')
                        originaldataframe = merged
                # print('massdataframe')
                # print(originaldataframe)
                originaldataframe.fillna(value=999, inplace=True)
                originaldataframe.sort(columns='mass', inplace=True)
                originaldataframe.to_csv('/home/trummelbummel/Desktop/resultingmasses.csv', sep=',')
                originaldataframe = originaldataframe.apply(np.log)
                return originaldata




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
            # print(symbolicTS)
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
                for i in range(len(periodicityhintarray)):
                    period_j = i
                    print('periodicitysearch for period with length : ' + str(periodicityhintarray[period_j]))
                # workingcube = fingerprinter.constructWorkingcubediscretized(normalizedTS,cube, dictcube,  periodicityhintarray, period_j, mode, numbins)
                    workingcube = self.constructWorkingcubeComplete(normalizedTS, cube, dictcube, binarycube, periodicityhintarray, period_j, mode, numbins)
                    Tslice = self.constructTslice(workingcube, periodicityhintarray, period_j)
                    print(Tslice.keys())
                # # apriori mining of patterns for maxsubpatternstree construction
                # # minconfidencepattern here can be low because we are mining partial periodicities i.e. can also only occure rarely
                # TODO this threshold here in mineF1 results in empty set but should not be an issue mayebe minconfidence here is not even relevant
                    F1, F1counts = self.mineF1(Tslice, minconfidencepattern, periodicityhintarray, period_j)
                    print('Frequent one cycles')
                    print(F1)
                    print(len(F1))
                    if len(F1) != 0:
                        # TODO: remove break from here
                        maxsubpatterntree = self.mineMaxSubpatternTree(Tslice, F1, F1counts, periodicityhintarray, period_j, numbins, nodelist)
                        if space_type == 'mass':
                            print('doing massrun')
                            frequentperiods = self.mineFrequentPeriodicitiesMass(Tslice, maxsubpatterntree, minconf, periodicityhintarray, period_j, normalizedTS)
                            # break
                        elif space_type == 'freq':
                            print('doing FREQ run')
                            frequentperiods = self.mineFrequentPeriodicitiesFreqs(Tslice, maxsubpatterntree, minconf, periodicityhintarray, period_j, normalizedTS)
                            # break
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
# Mine frequent patterns from Tree
#
###########################################################

        def mineFrequentPeriodicitiesFreqs(self, Tslice, maxsubpatterntreeinstance, minconfidence, periods, period_j, normalizedTS):
            numperiods = float(len(normalizedTS))/float(periods[period_j])
            frequencies = maxsubpatterntree.getFrequencies(maxsubpatterntreeinstance, minconfidence, numperiods)
            masses = Tslice['m/z']
            if len(frequencies) == 0:
                print('no frequent patterns found for periodlength ' + str(periods[period_j]) + ' try different parametersettings.' )
            else:
                for i in frequencies:
                    # exclude maxpattern Cmax incidences from results
                    if type(i[2]) != dict:
                        # print(' i am here ending')
                        print('pattern is : ' + str(i[2]))
                        print('confidence of pattern is: ' + str(i[0]))
                        print('frequency of pattern is: ' + str(i[1]))
                        print('first period occurence starts at index : ' + str(i[3]))
                        print('occurences at indizes :' + str(i[5]))
                        # retrieve masses where the pattern occurs i.e. where there are non* letters
                        print('periodic freqs :' + str(i[4]))
                        print('periodlength is : ' + str(len(i[2])))

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
                        # print(' i am here ending')
                        print('periodic frequencies : ' + str(i[2]))
                        print('confidence of pattern is: ' + str(i[0]))
                        print('frequency of pattern is: ' + str(i[1]))
                        print('first period occurence starts at index : ' + str(i[3]))
                        print('occurences at indizes :' + str(i[5]))
                        # retrieve masses where the pattern occurs i.e. where there are non* letters
                        print('periodic masses :' + str(i[4]))
                        print('periodlength is : ' + str(periods[period_j]))



############################################################
#
# Construct max subpattern tree
#
###########################################################


        def mineMaxSubpatternTree(self, Tslice, F1, F1counts, periods, period_j, numbins, nodelist):
            print('constructing Tree...')
            counter = 1
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
                    # gets the maximum subpattern hitset in second scan of TS
                    maxsubpatternreps, transactions, hitsetcomplete, maxsubpatterntreeinstance = maxsubpatternhitset.getmaxPatterns(treedepth, Tslice, maxsubpatterntreeinstance,  counterindexpassed, Cmax, periods, period_j, nonstar, nodelist)
                    transactions = None
                    # returns nonstar elements count
                    nonstar = len(hitsetcomplete[0][0].nonzero()[0])

                    # maxsubpatterns, maxsubpatterntreeinstance, CmaxHitset = maxsubpatternhitset.findCiHit(Tslice, maxsubpatternreps, Cmax, maxsubpatterntreeinstance, periods, period_j, transactions, counterindexpassed=counterindexpassed)
                    counter += 1
                    counterindexpassed += 1

                    # hitset = CmaxHitset
                    if len(hitsetcomplete) == 0:
                        print('length hitset is 0')
                        boolean = False
                        # onecyclepatterns = maxsubpatternhitset.form1cyclepatterns(F1counts, periods, period_j)
                        # maxsubpatterntreeinstance = maxsubpatternhitset.addLeafnodes(onecyclepatterns, maxsubpatterntreeinstance)
                        filename = 'maxsubpatterntreeperiodlength' + str(periods[period_j]) + '.txt'
                        maxsubpatterntreeinstance.save2file(filename)
                        return maxsubpatterntreeinstance
                        #TODO: counter == 2 case fix it analogous to counter==1
                if counter == 2:
                        nonstar = nonstar
                        if nonstar > 2:
                            maxsubpatterntreeinstance = maxsubpatternhitset.Minesubpatterns(maxsubpatterntreeinstance)
                            nonstar = len(hitsetcomplete[0][0].nonzero()[0])
                        else:
                            print('nonstar smaller 2')
                            boolean = False
                            # onecyclepatterns = maxsubpatternhitset.form1cyclepatterns(F1counts, periods, period_j)
                            # maxsubpatterntreeinstance = maxsubpatternhitset.addLeafnodes(onecyclepatterns, maxsubpatterntreeinstance)
                            filename = 'maxsubpatterntreeperiodlength' + str(periods[period_j]) + '.txt'
                            maxsubpatterntreeinstance.save2file(filename)
                            return maxsubpatterntreeinstance



                    # else:
                    #     print('nonstar smaller equals 2')
                    #     boolean = False
                    #     # onecyclepatterns = maxsubpatternhitset.form1cyclepatterns(F1counts, periods, period_j)
                    #     # maxsubpatterntreeinstance = maxsubpatternhitset.addLeafnodes(onecyclepatterns, maxsubpatterntreeinstance)
                    #     filename = 'maxsubpatterntreeperiodlength' + str(periods[period_j]) + '.txt'
                    #     # maxsubpatterntreeinstance.save2file(filename)
                    #     return maxsubpatterntreeinstance

                # if counter == 2:
                #     print('constructing second level of tree...')
                #
                #     leaves = maxsubpatterntreeinstance.leaves('root')
                #     previousnodes = []
                #     for i in leaves:
                #         previousnodes.append(i.data.pattern)
                #     Ci, lengthCi = maxsubpatternhitset.formCi(hitsetcomplete, periods, period_j, counter)
                #
                #
                #     # Cmax, lengthCmax = maxsubpatternhitset.formCmax(F1, period_j)
                #     nonstar = (len(Ci)) - 1
                #     # parent is none for root
                #     treedepth = maxsubpatterntreeinstance.depth()
                #     if nonstar > 2:
                #
                #         maxsubpatternreps,  transactions, hitsetcomplete = maxsubpatternhitset.getPatterns(treedepth, previousnodes, Tslice, maxsubpatterntreeinstance,  counterindexpassed,Cmax, periods, period_j, nonstar, numbins, nodelist)
                #         # maxsubpatterns, maxsubpatterntreeinstance, CmaxHitset = maxsubpatternhitset.findCiHit(Tslice, maxsubpatternreps, Cmax, maxsubpatterntreeinstance, periods, period_j, transactions, counterindexpassed=counterindexpassed)
                #         counterindexpassed += 1
                #         # treedict = maxsubpatterntree.todict(maxsubpatterntreeinstance)
                #         # counter += 1
                #         # input for next iteration
                #         if len(hitsetcomplete) == 0:
                #             print('hitset empty')
                #             boolean = False
                #             # onecyclepatterns = maxsubpatternhitset.form1cyclepatterns(F1counts, periods, period_j)
                #             # maxsubpatterntreeinstance = maxsubpatternhitset.addLeafnodes(onecyclepatterns, maxsubpatterntreeinstance)
                #             filename = 'maxsubpatterntreeperiodlength' + str(periods[period_j]) + '.txt'
                #             # maxsubpatterntreeinstance.save2file(filename)
                #             return maxsubpatterntreeinstance
                #         else:
                #             counter+= 1
                #             continue
                #     else:
                #         print('nonstar samller equals 2')
                #         boolean = False
                #         # onecyclepatterns = maxsubpatternhitset.form1cyclepatterns(F1counts, periods, period_j)
                #         # maxsubpatterntreeinstance = maxsubpatternhitset.addLeafnodes(onecyclepatterns, maxsubpatterntreeinstance)
                #         filename = 'maxsubpatterntreeperiodlength' + str(periods[period_j]) + '.txt'
                #         # maxsubpatterntreeinstance.save2file(filename)
                #         return maxsubpatterntreeinstance
                #
                # # from second iteration onwards if hitset of previous run not empty and nonstar characters more than 2 in pattern
                # after second scan of time series we need to derive the counts of subpatterns from the nodes which can
                # generate the subpatterns
                # if counter > 1 :
                #     print('Constructing  level '  + str(counter) + ' of tree...')
                #     # TODO: bring hitset into original format again and then feed it to here
                #     # pass it from
                #     # previoushitsetrepresentation = miningperiodicpatterns.formhitset(hitset, period_j, counter)
                #     Ci, lengthCi = maxsubpatternhitset.formCi(hitsetcomplete, periods, period_j, counter)
                #
                #     nonstar = len(Ci)-(counter-1)
                #     treedepth = treedepth + 1
                #     leaves = maxsubpatterntreeinstance.leaves('root')
                #     previousnodes = []
                #     for i in leaves:
                #         previousnodes.append(i.data.pattern)
                #     # compute number of nonzero entries i.e. nonstar
                #     pattern = leaves[0].data.pattern
                #
                #
                #     if nonstar > 2 and len(hitsetcomplete) != 0:
                #         maxsubpatternreps,  transactions, hitsetcomplete = maxsubpatternhitset.getPatterns(treedepth, previousnodes, Tslice, maxsubpatterntreeinstance,  counterindexpassed,Cmax, periods, period_j, nonstar, numbins, nodelist)
                #         nonstar -= 1
                #         # # counter is also indicating at which i-pattern we are with i being the integer value denoting
                #         # # how many non* letters there are in a subpattern
                #         # counter += 1
                #         if len(hitsetcomplete) == 0:
                #             boolean = False
                #             # one cycle patterns not needed for chemical applications because the distance will usually be too large
                #             # onecyclepatterns = maxsubpatternhitset.form1cyclepatterns(F1counts, periods, period_j)
                #             # maxsubpatterntreeinstance = maxsubpatternhitset.addLeafnodes(onecyclepatterns, maxsubpatterntreeinstance)
                #             filename = 'maxsubpatterntreeperiodlength' + str(periods[period_j]) + '.txt'
                #             # maxsubpatterntreeinstance.save2file(filename)
                #             return maxsubpatterntreeinstance
                #         else:
                #             counter += 1
                #             continue
                #     else:
                #         print('nonstar smaller equals 2')
                #         boolean = False
                #         # onecyclepatterns = maxsubpatternhitset.form1cyclepatterns(F1counts, periods, period_j)
                #         # maxsubpatterntreeinstance = maxsubpatternhitset.addLeafnodes(onecyclepatterns, maxsubpatterntreeinstance)
                #         filename = 'maxsubpatterntreeperiodlength' + str(periods[period_j]) + '.txt'
                #         # maxsubpatterntreeinstance.save2file(filename)
                #         return maxsubpatterntreeinstance
                #






############################################################
#
# Construction of Spectrogram
#
###############################################################

        def make_spectrograms(self, path, extensions, nprocesses=None):
                filenames_to_specmaker = []
                foldernames = []
                totalfiles = []
                #  Try to use the maximum amount of processes if not given.
                try:
                    nprocesses = nprocesses or multiprocessing.cpu_count()
                except NotImplementedError:
                    nprocesses = 1
                else:
                    nprocesses = 1 if nprocesses <= 0 else nprocesses
                dirnames = decoder.crawl_folders(path, extensions)
                for i in dirnames:
                    directory =  i + '/'
                    specmaker.specmaker_worker(directory)
                return

###################################################################
#
# do Clustering
#
##############################################################

        # def doClustering(self, preprocessingmode, method, mode, numclusters):
        #     cwd = os.getcwd()
        #     X, labels = decoder.readinputclustering(cwd + '/MSfingerprinter/results/globalfingerprintsnew.csv', preprocessingmode=preprocessingmode)
        #     print(X)
        #     if method == 'Kmeans':
        #         clustering.runKmeans(X, preprocessingmode, numclusters=numclusters)
        #     elif method == 'affinity':
        #         clustering.Affinitypropagation(X, labels)
        #     else:
        #         clustering.runHierachicalclustering(X, labels, preprocessingmode=preprocessingmode, mode=mode)

        def storeClusteringresults(self, outfile, listofclusters):

            with open(outfile, 'a') as f:
                writer = csv.writer(f)
                writer.writerow(listofclusters)
                f.close()



        def doClusteringASC(self, Data, preprocessingmode, method, mode, numclusters):
            # print(Data)
            # cwd = os.getcwd()
            # X, labels = decoder.readinputclustering(cwd + '/MSfingerprinter/results/globalfingerprintsnew.csv', preprocessingmode=preprocessingmode)
            X = Data

            # print('Data')
            # print(X)
            labels = list(X.columns)

            # print(labels)

            if method == 'Kmeans':
                clustering.runKmeans(X, preprocessingmode, numclusters=numclusters)
            elif method == 'affinity':
                clustering.Affinitypropagation(X, labels)
            else:
                clustering.runHierachicalclustering(X, labels, preprocessingmode=preprocessingmode, mode=mode)


        def doClusteringCSVfreq(self,k,  Data, colnames, preprocessingmode, method, mode, numclusters):
            # print(Data)
            # cwd = os.getcwd()
            # X, labels = decoder.readinputclustering(cwd + '/MSfingerprinter/results/globalfingerprintsnew.csv', preprocessingmode=preprocessingmode)
            X = Data.T
            # del X['freq']
            # print('Data')
            # print(X)
            labels = colnames
            # print(labels)
            if method == 'Kmeans':
                clustering.runKmeans(X, preprocessingmode, numclusters=numclusters)
            elif method == 'affinity':
                clustering.Affinitypropagation(X, labels)
            else:
                clusters = clustering.runHierachicalclustering(k, X, labels, preprocessingmode=preprocessingmode, mode=mode)
            return clusters

        def doClusteringCSVmass(self, k, Data,colnames, preprocessingmode, method, mode, numclusters):
            # print(Data)
            # cwd = os.getcwd()
            # X, labels = decoder.readinputclustering(cwd + '/MSfingerprinter/results/globalfingerprintsnew.csv', preprocessingmode=preprocessingmode)
            X = Data.T
            # del X['mass']
            # print('Data')
            # print(X)
            labels = colnames
            # print(labels)
            if method == 'Kmeans':
                clustering.runKmeans(X, preprocessingmode, numclusters=numclusters)
            elif method == 'affinity':
                clustering.Affinitypropagation(X, labels)
            else:
                clusters = clustering.runHierachicalclustering(k, X, labels, preprocessingmode=preprocessingmode, mode=mode)
            return clusters
