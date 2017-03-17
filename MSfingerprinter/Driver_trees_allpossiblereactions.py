import os
import pandas as pd
import numpy as np
import json
import sys
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
import MSfingerprinter.reactiontreecompletespace as reactiontreecompletespace
import MSfingerprinter.postprocessing as postprocessing


'''
This script allows the construction of a tree covering parts of the search space in NOM considering only CHO
repeat units to be fragmented. The start at a root which is a reference substance used. Examples can be found in
results/completeInitiatortrees.
This approach however requires large amount of memory as it is a k-ary tree. Hence it
is only feasible to construct such trees without duplicate nodes and for a certain mass range.
'''

# instantiate MSfingerprinter
fingerprinter = MSF.MSfingerprinter()

# possible reactions and masses NOM
H = 1.007825
C = 12.0
O = 15.994915
# possible reactions

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


# SRFA calabration ESI negative mode

Initiators =    [[('C9H5O7',	225.0040764)]	,
                [('C11H7O7',	251.0197264)],
                [('C14H11O6',	275.0561117)],
                [('C15H9O7',	301.0353764)],
                [('C15H17O8',	325.092891)],
                [('C17H19O8',	351.108541)],
                [('C18H15O9',	375.0721556)],
                [('C19H13O10',	401.0514202)],
                [('C19H21O11',	425.1089348)],
                [('C21H23O11',	451.1245848)],
                [('C22H19O12',	475.0881995)],
                [('C23H17O13',	501.0674641)],
                [('C25H17O13',	525.0674641)],
                [('C27H19O13',	551.0831141)],
                [('C26H23O15',	575.1042433)],
                [('C27H21O16',	601.083508)],
                [('C29H21O16',	625.083508)],
                [('C31H23O16',	651.099158)],
                [('C33H23O16',	675.099158)],
                [('C31H25O19',	701.0995518)]]


# build all trees (computationally expensive)
nodelist = [('subH2O', subH2O), ('addH2O', H2O), ('subCH2', subCH2), ('addCH2', CH2), ('subH2', subH2), ('addH2', H2), ('subCO', subCO), ('addCO', CO), ('subCO2', subCO2),('addCO2', CO2), ('subO', subO), ('addO', O)]


# massrange gives maximum level of tree for each type of repeating units (computed by function maxlevel)
maxmass = 205
minmass = 300
massrange = abs(maxmass - minmass)

# create all trees possible combinations of Initiator and repeating unit
for i in Initiators:
    Initiator = i
    print('constructing tree for ' + Initiator[0][0] + 'allreactions....')
        # this creates initiatortrees
    Initiatortree = fingerprinter.constructInitiatorTreescompletespace(Initiator, nodelist, massrange)
