import os
import pandas as pd
import numpy as np
import fnmatch
import matplotlib.pyplot as plt
# from MSfingerprinter import *
import MSfingerprinter as MSF
import MSfingerprinter.decoder as decoder
import MSfingerprinter.preprocessing as preprocessing
import MSfingerprinter.pysax as SAX
import MSfingerprinter.periodicityfinder as periodicityfinder
import MSfingerprinter.maxsubpatterntree as maxsubpatterntree
import MSfingerprinter.miningperiodicpatterns as miningperiodicpatterns
import MSfingerprinter.maxsubpatternhitset as maxsubpatternhitset
import MSfingerprinter.postprocessing as postprocessing


'''
postprocessing of patterns found, saves relevant mass differences in groundtruth.csv file
'''

# path input file is output of periodicity mining algorithm saved in text files, by default stored in current working directory
pathpatternresults = os.getcwd()
extensionspatterns = ['.txt']
extensions = ['.json']

# instantiate MSfingerprinter
fingerprinter = MSF.MSfingerprinter()

#######################################3
#
#
# mass patterns simple finding in CHO space
#
############################################3

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

# list of repeating units for pruning patterns that are chemically irrelevant
nodelist = np.array([ H2O,   CH2,  H2,  CO,  CO2, O])

extensionspatterns = ['.txt']



def findMassdifferences(MSarraymass, nodelist):
    '''
    finds mass differences that are relevant in CHO space, nodelist can be adapted for other mass differences if needed
    '''
    MSdifferences = []
    nodelistrounded = [round(elem, 6) for elem in nodelist]
    n = len(MSarraymass)
    MSarraymass = list(np.array(MSarraymass).reshape(1,n)[0])
    # removes elements for differences between peaks further apart
    for i in range(len(MSarraymass)):
        massi = MSarraymass[i]
        for j in range(len(MSarraymass)):
            massj = MSarraymass[j]
            diff = massi - massj
            if j != i:
                for x in nodelistrounded:
                    # check if the difference is valid i.e. if repeatunit division yields no remainder
                    # allow for error of 9 ppm
                    if round(diff,6)%x < 0.000009:
                        #print(massi, massj, x)
                        # appends mass i and position i, massj and position j and reactiontype as mass
                        MSdifferences.append([massi, i, massj, j, x])
    return MSdifferences

def find_files(path, extensions):
    # Allow json
    extensions = [e.replace(".", "") for e in extensions]
    for dirpath, dirnames, files in os.walk(path):
        for extension in extensions:
            for f in fnmatch.filter(files, "*.%s" % extension):
                p = os.path.join(dirpath, f)
                print(p)
                yield (p)

counter = 0
Dataframecomplete= pd.DataFrame()
Dataframecompleteirrelevant = pd.DataFrame()
# instantiate MSfingerprinter
fingerprinter = MSF.MSfingerprinter()
patternfilearray = list(find_files(pathpatternresults, extensionspatterns))
print('files to analyse : ')
print(patternfilearray)
print('\n')
allothers = 0
counternomassdifferencesfound = 0
countermassdifferencesfound = 0
completeDataframe = pd.DataFrame()
massdifferencedataframes = []

for patternfile in patternfilearray:
    allpatterns = fingerprinter.doPostProcessingSingleAlldiffs(patternfile, nodelist)
    counter += 1

for i in allpatterns:
    masspoints = i[0][0]
    massdifferences = findMassdifferences(masspoints, nodelist)
    massdifferenceDataframe = pd.DataFrame(massdifferences)
    massdifferencedataframes.append(massdifferenceDataframe)
    if len(massdifferences) == 0:
        counternomassdifferencesfound += 1
    else:
        countermassdifferencesfound += 1
dataframeMD = pd.concat(massdifferencedataframes, axis = 0)
dataframeMD.to_csv('groundtruth.csv')

print('total number of patterns ' + str(len(allpatterns)))
print('no massdifferences found in ' + str(counternomassdifferencesfound) + ' patterns')
print('massdifferences found in ' + str(countermassdifferencesfound) + ' patterns')
