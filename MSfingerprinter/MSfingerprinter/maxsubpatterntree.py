from __future__ import unicode_literals
from treelib import Node, Tree
import numpy as np
import collections
import json


# dataobject passed to node
# dataobject passed to node
class subpattern(object):
    def __init__(self, pattern, count,   onset, occurences, timepointsperiodicities):
        self.pattern = pattern
        self.count = count
        # self.transaction = transaction
        self.onset = onset
        self.occurences = occurences
        self.timepointsperiodicities = timepointsperiodicities

def createChildnode(pattern, name, maxsubpatterntreeinstance,  parent, counter, onset, occurences, timepointsperiodicities):
    maxsubpatterntreeinstance.create_node(name, name, parent=parent,data=subpattern(pattern, counter,   onset, occurences, timepointsperiodicities))
    return maxsubpatterntreeinstance

def createChildnode1cycle(pattern, name, maxsubpatterntreeinstance,  parent, counter):
    onset = None
    occurence = []
    timepointsperiodicities = []
    maxsubpatterntreeinstance.create_node(name, name, parent=parent,data=subpattern(pattern, counter,   onset, occurence, timepointsperiodicities))
    return maxsubpatterntreeinstance

def traverseandgetData(maxsubpatterntreeinstance):
    allnodes = maxsubpatterntreeinstance.expand_tree('root')
    listallnodes = list(allnodes)
    for i in listallnodes:
        currentnode = maxsubpatterntreeinstance.get_node(i)
        print('pattern', currentnode.data.pattern)
        print('name', currentnode.tag)
        print('count', currentnode.data.count)
        print('onset', currentnode.data.onset)


def getFrequencies(maxsubpatterntreeinstance, confidencethreshold, numperiods):
    allpaths = []
    wrongpaths = []
    frequencylist = []
    pathlist = maxsubpatterntreeinstance.paths_to_leaves()
    paths = list(pathlist)
    for path in paths:
        freq = 0
        previouscount = 0
        wrongpath = []
        for j in path:
            node = maxsubpatterntreeinstance.get_node(j)
            count = node.data.count
            pattern = node.data.pattern
            onset = node.data.onset
            occurences = node.data.occurences
            timepoints = node.data.timepointsperiodicities
            # only count if it is not the same pattern (subpattern that has already been counted)
            if count != previouscount:
                freq = count
                confidence = float(freq)/float(numperiods)
            if confidence  > confidencethreshold:
                # freq -1 to deduct root node count
                frequencylist.append((confidence, freq, pattern, onset, occurences, timepoints))

    return frequencylist


def createRootnode(Cmax, parent):
    maxsubpatterntreeinstance = Tree()
    # parameters nodename, nodetype, data contained in node nodeID = 0
    onset = None
    occurence = None
    timepointsperiodicities = None
    maxsubpatterntreeinstance.create_node('Cmax', 'root', parent=parent, data=subpattern(Cmax, 1, onset, occurence, timepointsperiodicities))
    return maxsubpatterntreeinstance

def showPattern(maxsubpatterntreeinstance):
    print('pattern')
    print(maxsubpatterntreeinstance.show(data_property='pattern'))


def showCount(maxsubpatterntreeinstance):
    print('count')
    print(maxsubpatterntreeinstance.show(data_property="count"))
