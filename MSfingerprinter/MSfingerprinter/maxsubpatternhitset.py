import numpy as np
from treelib import Node, Tree
import MSfingerprinter.maxsubpatterntree as maxsubpatterntree
import miningperiodicpatterns
import operator
import itertools as it
import re
from collections import defaultdict, OrderedDict
import itertools

'''algorithm 3.2 in Efficient mining of partial periodic patterns in time series database Han et al.
form candidate frequent maxpattern maximal pattern generated from F1 uses regex to handle optional events
regex implementation of python can only handle 100 capturing group such that the maximum periodicity that can
be investigated with this implementation is 100 '''

def formCmax(F1, period_j):
    Cmax = []
    previous = 0
    # sort such that second element of (position, valuebin) tuple is ascending
    # so sort according to binvalues
    F1 = sorted(list(F1), key=lambda x: x[2][1])
    for i in F1:
         listcandidatepattern = list(i)
         # position in the period
         position = listcandidatepattern[2][0]
         # bin the value is in
         value = listcandidatepattern[2][1]
         Cmax.append((position, value))
    lengthCmax = len(Cmax)
    # Cmaxdict contains set of values at certain positions in the periodic pattern
    Cmaxdict = {}
    for key, value in Cmax:
        Cmaxdict.setdefault(key, set()).add(value)

    return Cmaxdict, lengthCmax


#test for exact equality in hitset
def arreq_in_list(maskedsubpattern, hitset):
    return next((True for elem in hitset if np.ma.allequal(elem, maskedsubpattern)), False)

# forms the 1 cycle patterns for tree base case fill the leaf level of the tree with 1-cycle patterns
def form1cyclepatterns(F1, periods, period_j):
    F1 = list(F1)
    periodlen = periods[period_j]
    onecyclepatterns = []
    periodlength = F1[0][0]
    for i in F1:
        pattern = np.zeros((periodlen,))
        np.put(pattern, [i[3][0]], [i[3][1]])
        onecyclepatterns.append((i[0], pattern))
    return onecyclepatterns


def createRegexfromCmax(Cmax):

    noncapturinggroups = []
    # one has to build regex as a string
    pattern = ''
    insert = []
    counter = 0
    for k,v in Cmax.items():
        subinsert = []
        if v == set([0]):
            counter += 1
            ins = '(?:' + '[0-9]'  + '{0,1}\,' + ')'
            insert.append(ins)

            noncapturinggroups.append(counter)
        elif v != set([0]):
            counter += 1
            for i in range(len(v)):

                subinsert.append(str(i) + '{0,1}')
                groupinsert = '|'.join(subinsert)
                ins = '(' + groupinsert  + ')\,'
            insert.append(ins)
    insert.append('?')
    regexpattern = ''.join(insert)

    return regexpattern, noncapturinggroups

def findhitregex(regexpattern, Segment_i):
    # replace anything that is not matched by pattern with 0
    subpattern = re.sub(regexpattern, '0', Segment_i)
    return subpattern




def findCiHitRegex(treedepth, previousnodes, Tslice, maxsubpattern, regexpattern, noncapturinggroups, maxsubpatterntreeinstance, periods, period_j,    counterindexpassed, nodelist):
    # leaves are last frequent patterns from which to generate new frequent subpatterns
    # occurences gets emptied out each pattern iteration...problem!
    occurences = []
    equalitybool = False
    previoustimes = np.zeros((1, period_j))
    hitset = []
    timepointsperiodicities = []
    globcounter = 0
    boolenarray = False
    onset = None
    periodlength = periods[period_j]
    testtimeseries = Tslice['symbolicTS']
    masspoints = Tslice['m/zPeriods']
    freqpoints = Tslice['freqPeriods']
    timepoints = Tslice['timepointsPeriods']

    # creates a regular expression object from regexpattern that represents Cmax
    regularexpressionobject = re.compile(regexpattern)
    counter = 0

    # compare each subpattern with all periodsegments
    if regularexpressionobject:
        if counterindexpassed == 0:
            for k,v in testtimeseries.items():
                masses = masspoints[k]
                times = timepoints[k]
                freqs = freqpoints[k]
                discreteTSperiod = testtimeseries[k]
                vint = v.astype('int')
                v = vint.astype('str')
                listarray = v.tolist()
                # TODO: save the discrete TS already as list and not as np.array
                stringpattern = ','.join(listarray)
                allpatterns = regularexpressionobject.findall(stringpattern)
                if len(allpatterns) != 0:
                        newnodepattern = list(allpatterns[0])
                        # insert star at all noncapturinggroups (-1 cause counter is 1 greater)
                        for i in noncapturinggroups:
                            index = i-1
                            newnodepattern.insert(index,'*')
                        # if hitset contains elements then check if it is already there, if not create node
                        if len(hitset) != 0:
                            for i in hitset:
                                name = ''.join(newnodepattern)
                                # if not in hitset then create node
                                if not maxsubpatterntreeinstance.contains(name):
                                    newnodepattern = np.array(newnodepattern)
                                    hitset.append([newnodepattern, 1])
                                    nonstarindizes = np.where(newnodepattern != '*')
                                    counter = 1
                                    parent = 'root'
                                    timepointsperiod = [np.take(times, nonstarindizes).tolist()]
                                    occuringat = [np.take(masses, nonstarindizes).tolist()]
                                    # this takes the freq from periodic pattern indizes
                                    frequencyoccuringat = np.take(freqs,nonstarindizes).tolist()
                                    print('creating node')
                                    maxsubpatterntreeinstance = maxsubpatterntree.createChildnode(frequencyoccuringat, name, maxsubpatterntreeinstance,  parent, counter, onset, occuringat, timepointsperiod)
                                # else increment count at node
                                else:
                                    index = -1
                                    for j in hitset:
                                        index += 1

                                        name = ''.join(newnodepattern)
                                        oldpatternname = ''.join(j[0])

                                        if name == oldpatternname:
                                            hitset[index][1] == hitset[index][1] + 1
                                            node = maxsubpatterntreeinstance.get_node(name)
                                            node.data.count = hitset[index][1]

                        else:

                            newnodepattern = np.array(newnodepattern)
                            hitset.append([newnodepattern, 1])
                            nonstarindizes = np.where(newnodepattern != '*')
                            counter = 1
                            parent = 'root'
                            name = ''.join(newnodepattern)
                            timepointsperiod = [np.take(times, nonstarindizes).tolist()]
                            occuringat = [np.take(masses, nonstarindizes).tolist()]
                            # this takes the freq from periodic pattern indizes
                            frequencyoccuringat = np.take(freqs,nonstarindizes).tolist()
                            print('creating node')
                            maxsubpatterntreeinstance = maxsubpatterntree.createChildnode(frequencyoccuringat, name, maxsubpatterntreeinstance,  parent, counter, onset, occuringat, timepointsperiod)

                else:
                    return maxsubpatterntreeinstance, hitset
            # maxsubpatterntreeinstance.show()
            return maxsubpatterntreeinstance, hitset



# now set approach intersect Cmax with each period segment
def  getmaxPatternsset(treedepth, Tslice, maxsubpatterntreeinstance, counterindexpassed, Cmax, periods, period_j, nonstarcount, nodelist):
    hitsetcomplete = []
    parent = None
    transactions = []
    counter = 0
    subpatterndict = OrderedDict()
    # function call to generate all possible patterns from Ci
    subpatterns = []
    subpatternstrings = []
    # keys are running number of index
    possiblekeys = range(periods[period_j])
    possiblekeys = set([float(i) for i in possiblekeys])
    # check in dictionary if one of the running number keys is missing and set the corresponding position to 0 always
    presentkeys = set(Cmax.keys())
    # setdifference to find nonmatching elements
    unmatchedkeys = possiblekeys.symmetric_difference(presentkeys)
    unmatchedkeyslist = list(unmatchedkeys)
    # insert set([0]) wherever key is missing (will always be don't care character)
    if unmatchedkeyslist:
        for k,v in Cmax.items():
            for i in unmatchedkeys:
                Cmax[i] = set([0])

    regexpattern, noncapturinggroups = createRegexfromCmax(Cmax)
    previousnodes = []
    maxsubpatterntreeinstance, hitset = findCiHitRegex(treedepth, previousnodes, Tslice, Cmax, regexpattern, noncapturinggroups, maxsubpatterntreeinstance, periods, period_j, counterindexpassed, nodelist)

    return maxsubpatterntreeinstance, hitset
