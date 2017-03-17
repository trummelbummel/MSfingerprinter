import os
import pandas as pd
import numpy as np
import fnmatch
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
import MSfingerprinter.comparetrees as comparetrees
import json
import sys
import matplotlib.pyplot as plt

'''
Writes masses and sumformula found in tree to output file groundtruthsumformulaend.csv and groundtruthsumformulastart.csv
results can be found in subdirectory /results/nodesfound/
'''


# logs both to terminal and to outputfile LOGGER
filenameresults = 'resultsgroundtruthtreefinder.txt'


class Logger(object):
    def __init__(self, filename="Default.log"):
        self.terminal = sys.stdout
        self.log = open(filenameresults, "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
sys.stdout = Logger(filenameresults)

pathpatternresults = os.getcwd() +'/results/nodesfortreesearch/uniquetreenodescorrectconfidence.csv'
extensionspatterns = ['.txt']
extensions = ['.json']


path = os.getcwd() + '/results/completeInitiatortrees/'
# instantiate MSfingerprinter
fingerprinter = MSF.MSfingerprinter()
groundtruthpatterns = pd.read_csv(pathpatternresults, sep=',')
meaningfulMassdifferences = groundtruthpatterns.values.tolist()
patternsfoundintreesend = fingerprinter.searchinTreesgroundtruth(meaningfulMassdifferences, path, extensions)
patternsfoundintreesstart = fingerprinter.searchinTreesgroundtruthstart(meaningfulMassdifferences, path, extensions)

# drop column that contains mass that is not annotated in this case
headers = ['initiator', 'mass1', 'mass2', 'stoichiometricformula']
resultsincludingformula = pd.DataFrame(patternsfoundintreesend, columns=headers)
resultsincludingformula.to_csv(os.getcwd() + '/results/nodesfound/groundtruthsumformulaend.csv', sep=',', encoding='ascii')

headers = ['initiator', 'mass1', 'mass2', 'stoichiometricformula']
resultsincludingformula = pd.DataFrame(patternsfoundintreesstart, columns=headers)
resultsincludingformula.to_csv(os.getcwd() + '/results/nodesfound/groundtruthsumformulastart.csv', sep=',', encoding='ascii')
