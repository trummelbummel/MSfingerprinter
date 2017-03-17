import re
import ast
import csv
import itertools
import numpy as np
import os
import ntpath

postprocessedperiodicfreqs = []
postprocessedmasses = []
postpostfreqs = []
postpostmasses = []


def openresultssingle(resultfile1):
    periodicfreqs = []
    periodicmasses = []
    with open(resultfile1, 'r') as f1:
        lines1 = f1.readlines()
    for i in range(len(lines1)):
        if lines1[i].startswith('masses :'):
            periodicmasses.append(lines1[i].lstrip('masses :').rstrip('\n'))
        if lines1[i].startswith('periodic frequencies :'):
            periodicfreqs.append([lines1[i].lstrip('periodic frequencies :').rstrip('\n')])
    return periodicmasses, periodicfreqs

def writesimtofile(outfilesimilarity, similar, dissimilar):
    with open(outfilesimilarity, 'w') as outf:
        outf.write('similar\n\n')
        outf.write(str(similar)+ '\n\n')
        outf.write('dissimilar\n\n')
        outf.write(str(dissimilar) + '\n\n')


def retrieveperiodicmasses(periodicmasses):
    for j in periodicmasses:
        if len(j) > 2:
            j.lstrip('[').rstrip(']')
            try:
                listj = ast.literal_eval(j)
            except:
                continue
            if len(listj) != 0:
                postprocessedmasses.append(listj)
    return postprocessedmasses

def retrieveperiodicfreqs(periodicfreqs):
    for i in periodicfreqs:
        i = i[0]
        if len(i) > 2:

            listi = i
            if len(listi) != 0:
                postprocessedperiodicfreqs.append(listi)
    return postprocessedperiodicfreqs


def countmeaninfuldiffpatternsmass(meaningfuldiffpatternsmass, otherdiffpatterns):
    meaningfulcount = len(meaningfuldiffpatternsmass)
    notmeaningfulcount = len(otherdiffpatterns)
    return meaningfulcount, notmeaningfulcount

def tofile(resultfile1, meaningfuldiffpatternsmass, meaningfuldiffpatternsfreq, otherdiffpatterns, otherdiffpatternsfreq, meaningfulcount, notmeaningfulcount):
    name = ntpath.basename(resultfile1).rstrip('.txt')
    name = os.getcwd() + '/results/Postprocessedresultspatternmining/' + name + 'RESULT.txt'
    filename = open(name, 'w')
    filename.write("count of meaningfullpatterns :" + str(meaningfulcount))
    filename.write('\n')
    for i in range(len(meaningfuldiffpatternsmass)):
        filename.write("%s\n" % list(str(element) for element in meaningfuldiffpatternsmass[i]))
        filename.write("%s\n" % list(str(element) for element in meaningfuldiffpatternsfreq[i]))
        filename.write('\n')
    filename.write("count of NOT meaningfullpatterns :" + str(notmeaningfulcount))
    filename.write('\n')
    for j in range(len(otherdiffpatterns)):
        filename.write("%s\n" % list(str(element) for element in otherdiffpatterns[j]))
        filename.write("%s\n" % list(str(element) for element in otherdiffpatternsfreq[j]))
        filename.write('\n')
    return
