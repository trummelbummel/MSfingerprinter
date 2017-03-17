# implementation based on binary tree structure module in python https://pypi.python.org/pypi/bintrees/2.0.2

# from bintrees.cython_trees import FastBinaryTree
from collections import defaultdict
from treelib import Node, Tree
import numpy as np
import collections
import json
import re


class ReactionNode(dict):
    def __init__(self, reactionpath, cost, reactiontype, sumcost, treename):
        self.path = [] # contains reactionpath so far
        self.reactiontype = reactiontype # contains reaction type as string
        self.cost = cost # contains cost of current Reaction
        self.sum = sumcost # contains sum resulting from all reactions
        self.treename = treename
        self._dict_= self


def createRootsmassdiff(rootlist):
    trees = []
    for i in rootlist:
        path = []
        parent=None
        treename = 'Tree' + i[0]
        reactiontreeinstance = Tree()
        # reactiontree takes itself, reactionpathsofar, cost, reactiontype, and current sum
        key = 0
        path.append(i[0])
        cost = i[1]
        sumcost = i[1]
        reactiontreeinstance.create_node(i[0], 'root', parent=None,data=ReactionNode(cost, path, i[0], sumcost, treename))
        trees.append(reactiontreeinstance)
    return trees

# cerates next level in a single reaction tree
def createnextlevelmassdiff(massdifference, nodelist, reactiontreeinstance, parent, counter):
    reactiontrees = []
    leaves = reactiontreeinstance.leaves('root')
    for leaf in leaves:
        parent = leaf.identifier
        parentname = leaf.tag
        for i in nodelist:
            counter += 1
            cost = i[1]
            # total cost of reactions so far
            # round to 6 digits
            sumcost = round(float(i[1]) + float(leaf.data.sum), 6)
            massdifference = round(massdifference, 6)
            # reactions so far
            path = leaf.data.path
            path.append(i[0])
            treename = leaf.data.treename
            key = i[0] + '_'+ str(counter)
            name = i[0] + '_'+ str(counter)
            # create node only if reaction won't result in neg total cost
            # if we have found a pathway to a massdifference return the treeinstance
            if sumcost == massdifference:
                rootnode = reactiontreeinstance.get_node('root')
                path = rootnode.data.path
                firstelementpath = path[0]
                # print('sumcost == massdiff')
                path.insert(0, firstelementpath)
                return path
            elif sumcost > 0 and num1 > 0 and num2 > 0 and num3 > 0:
                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))
            elif sumcost != massdifference and counter > 500:
                return 'no tree found for this pattern after 20 reactionsteps'

    # every level return reactiontreeinstances created for next iteration
    return (reactiontreeinstance, counter)

def getmaxtreelevel(nodelist, massrange):
    for i in nodelist:
        print(i)
        if i[0] == 'addH2':
            cost = i[1]
            print(cost)
            maxlevel = int(massrange/abs(int(cost)))
        if i[0] == 'addCH2':
            cost = i[1]
            print(cost)
            maxleveltwo = int(massrange/abs(int(cost)))

    return maxlevel, maxleveltwo


def createRootsInitiators(Initiators):
    trees = []
    for i in Initiators:
        path = []
        parent=None
        treename = 'Tree' + i[0]
        reactiontreeinstance = Tree()
        # reactiontree takes itself, reactionpathsofar, cost, reactiontype, and current sum
        key = 0
        path.append(i[0])
        cost = i[1]
        sumcost = i[1]
        reactiontreeinstance.create_node(i[0], 'root', parent=None,data=ReactionNode(cost, path, i[0], sumcost, treename))
        trees.append(reactiontreeinstance)
    return trees

# regex retrieves nodes and adds/subtracts appropriate number of CHO atoms
def createnextlevelinitiator(nodelist, reactiontreeinstance,  counter, maxlevel, maxleveltwo, previousmaxdepth, initiatorname):
    leaves = reactiontreeinstance.leaves('root')
    maxdepth = reactiontreeinstance.depth()
    if maxdepth < maxlevel and maxdepth < maxleveltwo:
        if counter == 1:
            for leaf in leaves:
                parent = leaf.identifier
                previousname = leaf.data.reactiontype

                for i in nodelist:
                    counter += 1
                    name = ''
                    nameordered = ['0','0','0']

                    if i[0] == 'addCO2':
                        m = re.match(r'([A-Z])([0-9]+)([A-Z])([0-9]+)([A-Z])([0-9]+)', previousname)
                        let1 = m.group(1)
                        let2 = m.group(3)
                        let3 = m.group(5)
                        # TODO: use enumeration to check which letter is which
                        if let1 == 'C':
                            name =''
                            num1 = int(m.group(2)) + 1
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'C':
                            name =''
                            num2 = int(m.group(4)) + 1
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'C':
                            name =''
                            num3 = int(m.group(6)) + 1
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'O':
                            name =''
                            num1 = int(m.group(2)) + 2
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'O':
                            name =''
                            num2 = int(m.group(4)) + 2
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'O':
                            name =''
                            num3 = int(m.group(6)) + 2
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'H':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'H':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'H':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name

                        cost = i[1]
                        sumcost = round((float(i[1])) + float(leaf.data.sum), 6)
                        path = leaf.data.path
                        path.append(i[0])
                        treename = leaf.data.treename
                        key = "".join(nameordered)
                        name = key

                        if sumcost > 0 and num1 > 0 and num2 > 0 and num3 > 0:
                            if not reactiontreeinstance.contains(name) :
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))

                            else:
                                name = name + '_' + str(counter)
                                key = name
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))



                        else:
                            continue
                    if i[0] == 'addH2O':

                        m = re.match(r'([A-Z])([0-9]+)([A-Z])([0-9]+)([A-Z])([0-9]+)', previousname)

                        let1 = m.group(1)
                        let2 = m.group(3)
                        let3 = m.group(5)
                        if let1 == 'C':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'C':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'C':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'O':
                            name =''
                            num1 = int(m.group(2)) + 1
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'O':
                            name =''
                            num2 = int(m.group(4)) + 1
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'O':
                            name =''
                            num3 = int(m.group(6)) + 1
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'H':
                            name =''
                            num1 = int(m.group(2)) + 2
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'H':
                            name =''
                            num2 = int(m.group(4)) + 2
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'H':
                            name =''
                            num3 = int(m.group(6)) + 2
                            name += let3
                            name += str(num3)
                            nameordered[2] = name

                        cost = i[1]
                        sumcost = round((float(i[1])) + float(leaf.data.sum), 6)
                        path = leaf.data.path
                        path.append(i[0])
                        treename = leaf.data.treename
                        key = "".join(nameordered)
                        name = key

                        if sumcost > 0 and num1 > 0 and num2 > 0 and num3 > 0:
                            if not reactiontreeinstance.contains(name) :
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))

                            else:
                                name = name + '_' + str(counter)
                                key = name
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))

                        else:
                            continue
                    if i[0] == 'addO':

                        m = re.match(r'([A-Z])([0-9]+)([A-Z])([0-9]+)([A-Z])([0-9]+)', previousname)

                        let1 = m.group(1)
                        let2 = m.group(3)
                        let3 = m.group(5)
                        if let1 == 'C':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'C':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'C':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'O':
                            name =''
                            num1 = int(m.group(2)) + 1
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'O':
                            name =''
                            num2 = int(m.group(4)) + 1
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'O':
                            name =''
                            num3 = int(m.group(6)) + 1
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'H':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'H':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'H':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name

                        cost = i[1]
                        sumcost = round((float(i[1])) + float(leaf.data.sum), 6)
                        path = leaf.data.path
                        path.append(i[0])
                        treename = leaf.data.treename
                        key = "".join(nameordered)
                        name = key

                        if sumcost > 0 and num1 > 0 and num2 > 0 and num3 > 0:
                            if not reactiontreeinstance.contains(name) :
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))
                            else:
                                name = name + '_' + str(counter)
                                key = name
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))

                        else:
                            continue
                    if i[0] == 'addH2':

                        m = re.match(r'([A-Z])([0-9]+)([A-Z])([0-9]+)([A-Z])([0-9]+)', previousname)

                        let1 = m.group(1)
                        let2 = m.group(3)
                        let3 = m.group(5)
                        if let1 == 'C':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'C':

                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'C':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'O':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'O':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'O':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'H':
                            name =''
                            num1 = int(m.group(2)) + 2
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'H':
                            name =''
                            num2 = int(m.group(4)) + 2
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'H':
                            name =''
                            num3 = int(m.group(6)) + 2
                            name += let3
                            name += str(num3)
                            nameordered[2] = name

                        cost = i[1]
                        sumcost = round((float(i[1])) + float(leaf.data.sum), 6)
                        path = leaf.data.path
                        path.append(i[0])
                        treename = leaf.data.treename
                        key = "".join(nameordered)
                        name = key
                        if sumcost > 0 and num1 > 0 and num2 > 0 and num3 > 0:
                            if not reactiontreeinstance.contains(name):
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))
                            else:
                                name = name + '_' + str(counter)
                                key = name
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))

                        else:
                            continue
                    if i[0] == 'addCO':

                        m = re.match(r'([A-Z])([0-9]+)([A-Z])([0-9]+)([A-Z])([0-9]+)', previousname)

                        let1 = m.group(1)
                        let2 = m.group(3)
                        let3 = m.group(5)
                        if let1 == 'C':
                            name =''
                            num1 = int(m.group(2)) + 1
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        if let2 == 'C':
                            name =''
                            num2 = int(m.group(4))  + 1
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'C':
                            name =''
                            num3 = int(m.group(6)) + 1
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'O':
                            name =''
                            num1 = int(m.group(2)) + 1
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'O':
                            name =''
                            num2 = int(m.group(4)) + 1
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'O':
                            name =''
                            num3 = int(m.group(6)) + 1
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'H':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'H':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'H':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name


                        cost = i[1]
                        sumcost = round((float(i[1])) + float(leaf.data.sum), 6)
                        path = leaf.data.path
                        path.append(i[0])
                        treename = leaf.data.treename
                        key = "".join(nameordered)
                        name = key
                        if sumcost > 0 and num1 > 0 and num2 > 0 and num3 > 0:
                            if not reactiontreeinstance.contains(name) :
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))
                                #reactiontreeinstance.show()
                            else:
                                name = name + '_' + str(counter)
                                key = name
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))
                                #reactiontreeinstance.show()

                        else:
                            continue
                    if i[0] == 'addCH2':
                        m = re.match(r'([A-Z])([0-9]+)([A-Z])([0-9]+)([A-Z])([0-9]+)', previousname)
                        let1 = m.group(1)
                        let2 = m.group(3)
                        let3 = m.group(5)
                        if let1 == 'C':
                            name =''
                            num1 = int(m.group(2)) + 1
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'C':
                            name =''
                            num2 = int(m.group(4)) + 1
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'C':
                            name =''
                            num3 = int(m.group(6)) + 1
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'O':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'O':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'O':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'H':
                            name =''
                            num1 = int(m.group(2)) + 2
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'H':
                            name =''
                            num2 = int(m.group(4)) + 2
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'H':
                            name =''
                            num3 = int(m.group(6)) + 2
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        cost = i[1]
                        sumcost = round((float(i[1])) + float(leaf.data.sum), 6)
                        path = leaf.data.path
                        path.append(i[0])
                        treename = leaf.data.treename
                        key = "".join(nameordered)
                        name = key

                        if sumcost > 0 and num1 > 0 and num2 > 0 and num3 > 0:
                            if not reactiontreeinstance.contains(name) :
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))
                            else:
                                name = name + '_' + str(counter)
                                key = name
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))

                        else:
                            continue
                    if i[0] == 'subCO2':

                        m = re.match(r'([A-Z])([0-9]+)([A-Z])([0-9]+)([A-Z])([0-9]+)', previousname)

                        let1 = m.group(1)
                        let2 = m.group(3)
                        let3 = m.group(5)
                        if let1 == 'C':
                            name =''
                            num1 = int(m.group(2)) - 1
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'C':
                            name =''
                            num2 = int(m.group(4)) - 1
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'C':
                            name =''
                            num3 = int(m.group(6)) - 1
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'O':
                            name =''
                            num1 = int(m.group(2)) - 2
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'O':
                            name =''
                            num2 = int(m.group(4)) - 2
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'O':
                            name =''
                            num3 = int(m.group(6)) - 2
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'H':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'H':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'H':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        cost = i[1]
                        sumcost = round((float(i[1])) + float(leaf.data.sum), 6)
                        path = leaf.data.path
                        path.append(i[0])
                        treename = leaf.data.treename
                        key = "".join(nameordered)
                        name = key

                        if sumcost > 0 and num1 > 0 and num2 > 0 and num3 > 0:
                            if not reactiontreeinstance.contains(name) :
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))
                            else:
                                name = name + '_' + str(counter)
                                key = name
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))


                        else:
                            continue
                    if i[0] == 'subH2O':

                        m = re.match(r'([A-Z])([0-9]+)([A-Z])([0-9]+)([A-Z])([0-9]+)', previousname)

                        let1 = m.group(1)

                        let2 = m.group(3)
                        let3 = m.group(5)
                        if let1 == 'C':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'C':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'C':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'O':
                            name =''
                            num1 = int(m.group(2)) - 1
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'O':
                            name =''
                            num2 = int(m.group(4)) - 1
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'O':

                            name =''
                            num3 = int(m.group(6)) - 1
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'H':
                            name =''
                            num1 = int(m.group(2)) - 2
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'H':
                            name =''
                            num2 = int(m.group(4)) - 2
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'H':
                            name =''
                            num3 = int(m.group(6)) - 2
                            name += let3
                            name += str(num3)
                            nameordered[2] = name


                        cost = i[1]
                        sumcost = round((float(i[1])) + float(leaf.data.sum), 6)
                        path = leaf.data.path
                        path.append(i[0])
                        treename = leaf.data.treename
                        key = "".join(nameordered)
                        name = key

                        if sumcost > 0 and num1 > 0 and num2 > 0 and num3 > 0:
                            if not reactiontreeinstance.contains(name) :
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))
                            else:
                                name = name + '_' + str(counter)
                                key = name
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))
                                # reactiontreeinstance.show()



                        else:
                            continue
                    if i[0] == 'subO':

                        m = re.match(r'([A-Z])([0-9]+)([A-Z])([0-9]+)([A-Z])([0-9]+)', previousname)

                        let1 = m.group(1)
                        let2 = m.group(3)
                        let3 = m.group(5)
                        if let1 == 'C':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'C':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'C':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'O':
                            name =''
                            num1 = int(m.group(2)) - 1
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'O':
                            name =''
                            num2 = int(m.group(4)) - 1
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'O':
                            name =''
                            num3 = int(m.group(6)) - 1
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'H':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'H':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'H':

                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name

                        cost = i[1]
                        sumcost = round((float(i[1])) + float(leaf.data.sum), 6)
                        path = leaf.data.path
                        path.append(i[0])
                        treename = leaf.data.treename
                        key = "".join(nameordered)
                        name = key

                        if sumcost > 0 and num1 > 0 and num2 > 0 and num3 > 0:
                            if not reactiontreeinstance.contains(name) :
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))
                            else:
                                name = name + '_' + str(counter)
                                key = name
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))

                        else:
                            continue
                    if i[0] == 'subH2':
                        m = re.match(r'([A-Z])([0-9]+)([A-Z])([0-9]+)([A-Z])([0-9]+)', previousname)

                        let1 = m.group(1)
                        let2 = m.group(3)
                        let3 = m.group(5)
                        if let1 == 'C':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'C':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'C':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'O':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'O':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'O':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'H':
                            name =''
                            num1 = int(m.group(2)) -2
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'H':
                            name =''
                            num2 = int(m.group(4)) - 2
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'H':
                            name =''
                            num3 = int(m.group(6)) - 2
                            name += let3
                            name += str(num3)
                            nameordered[2] = name

                        cost = i[1]
                        sumcost = round((float(i[1])) + float(leaf.data.sum), 6)
                        path = leaf.data.path
                        path.append(i[0])
                        treename = leaf.data.treename
                        key = "".join(nameordered)
                        name = key

                        if sumcost > 0 and num1 > 0 and num2 > 0 and num3 > 0:
                            if not reactiontreeinstance.contains(name) :
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))
                            else:
                                name = name + '_' + str(counter)
                                key = name
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))

                        else:
                            continue
                    if i[0] == 'subCO':

                        m = re.match(r'([A-Z])([0-9]+)([A-Z])([0-9]+)([A-Z])([0-9]+)', previousname)

                        let1 = m.group(1)
                        let2 = m.group(3)
                        let3 = m.group(5)
                        if let1 == 'C':
                            name =''
                            num1 = int(m.group(2))  - 1
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'C':
                            name =''
                            num2 = int(m.group(4)) - 1
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'C':
                            name =''
                            num3 = int(m.group(6)) - 1
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'O':
                            name =''
                            num1 = int(m.group(2)) - 1
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'O':
                            name =''
                            num2 = int(m.group(4)) - 1
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'O':
                            name =''
                            num3 = int(m.group(6)) - 1
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'H':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'H':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'H':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name


                        cost = i[1]
                        sumcost = round((float(i[1])) + float(leaf.data.sum), 6)
                        path = leaf.data.path
                        path.append(i[0])
                        treename = leaf.data.treename
                        key = "".join(nameordered)
                        name = key
                        if sumcost > 0 and num1 > 0 and num2 > 0 and num3 > 0:
                            if not reactiontreeinstance.contains(name) :
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))
                            else:
                                name = name + '_' + str(counter)
                                key = name
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))

                        else:
                            continue
                    if i[0] == 'subCH2':

                        m = re.match(r'([A-Z])([0-9]+)([A-Z])([0-9]+)([A-Z])([0-9]+)', previousname)

                        let1 = m.group(1)
                        let2 = m.group(3)
                        let3 = m.group(5)
                        if let1 == 'C':
                            name =''
                            num1 = int(m.group(2)) - 1
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'C':
                            name =''
                            num2 = int(m.group(4)) - 1
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'C':
                            name =''
                            num3 = int(m.group(6)) - 1
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'O':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'O':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'O':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'H':
                            name =''
                            num1 = int(m.group(2)) - 2
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'H':
                            name =''
                            num2 = int(m.group(4)) - 2
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'H':
                            name =''
                            num3 = int(m.group(6)) - 2
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        cost = i[1]
                        sumcost = round((float(i[1])) + float(leaf.data.sum), 6)
                        path = leaf.data.path
                        path.append(i[0])
                        treename = leaf.data.treename
                        key = "".join(nameordered)
                        name = key

                        if sumcost > 0 and num1 > 0 and num2 > 0 and num3 > 0:
                            if not reactiontreeinstance.contains(name) :
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))
                            else:
                                name = name + '_' + str(counter)
                                key = name
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))

                        else:
                            continue
            previousmaxdepth = maxdepth
            return reactiontreeinstance, counter, previousmaxdepth
        else:
            leaves = reactiontreeinstance.leaves()
            for leaf in leaves:
                parent = leaf.identifier
                maxdepth = reactiontreeinstance.depth()
                previousname = leaf.tag

                for i in nodelist:
                    counter += 1

                    nameordered = ['0','0','0']
                    name =''
                    if i[0] == 'addCO2':

                        m = re.match(r'([A-Z])([0-9]+)([A-Z])([0-9]+)([A-Z])([0-9]+)', previousname)

                        let1 = m.group(1)
                        let2 = m.group(3)
                        let3 = m.group(5)
                        if let1 == 'C':
                            name =''
                            num1 = int(m.group(2)) + 1
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'C':
                            name =''
                            num2 = int(m.group(4)) + 1
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'C':
                            name =''
                            num3 = int(m.group(6)) + 1
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'O':
                            name =''
                            num1 = int(m.group(2)) + 2
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'O':
                            name =''
                            num2 = int(m.group(4)) + 2
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'O':
                            name =''
                            num3 = int(m.group(6)) + 2
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'H':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'H':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'H':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name

                        cost = i[1]
                        sumcost = round((float(i[1])) + float(leaf.data.sum), 6)
                        path = leaf.data.path
                        path.append(i[0])
                        treename = leaf.data.treename
                        key = "".join(nameordered)
                        name = key

                        if sumcost > 0 and num1 > 0 and num2 > 0 and num3 > 0 :
                            if reactiontreeinstance.contains(name):
                                name = name + '_' + str(counter)
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))


                                continue
                            else:
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))

                        else:
                            continue
                    if i[0] == 'addH2O':


                        m = re.match(r'([A-Z])([0-9]+)([A-Z])([0-9]+)([A-Z])([0-9]+)', previousname)

                        let1 = m.group(1)
                        let2 = m.group(3)
                        let3 = m.group(5)
                        if let1 == 'C':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'C':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'C':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'O':
                            name =''
                            num1 = int(m.group(2)) + 1
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'O':
                            name =''
                            num2 = int(m.group(4)) + 1
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'O':
                            name =''
                            num3 = int(m.group(6)) + 1
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'H':
                            name =''
                            num1 = int(m.group(2)) + 2
                            name += let1
                            name += str(num1)

                            nameordered[0] = name
                        elif let2 == 'H':
                            name =''
                            num2 = int(m.group(4)) + 2
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'H':
                            name =''
                            num3 = int(m.group(6)) + 2
                            name += let3
                            name += str(num3)
                            nameordered[2] = name

                        cost = i[1]
                        sumcost = round((float(i[1])) + float(leaf.data.sum), 6)
                        path = leaf.data.path
                        path.append(i[0])
                        treename = leaf.data.treename
                        key = "".join(nameordered)
                        name = key
                        if sumcost > 0 and num1 > 0 and num2 > 0 and num3 > 0:
                            if reactiontreeinstance.contains(name):
                                name = name + '_' + str(counter)
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))


                                continue
                            else:
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))

                        else:
                            continue
                    if i[0] == 'addO':

                        m = re.match(r'([A-Z])([0-9]+)([A-Z])([0-9]+)([A-Z])([0-9]+)', previousname)

                        let1 = m.group(1)
                        let2 = m.group(3)
                        let3 = m.group(5)
                        if let1 == 'C':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'C':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'C':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'O':
                            name =''
                            num1 = int(m.group(2)) + 1
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'O':

                            name =''
                            num2 = int(m.group(4)) + 1
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'O':
                            name =''
                            num3 = int(m.group(6)) + 1
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'H':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'H':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'H':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name

                        cost = i[1]
                        sumcost = round((float(i[1])) + float(leaf.data.sum), 6)
                        path = leaf.data.path
                        path.append(i[0])
                        treename = leaf.data.treename
                        key = "".join(nameordered)
                        name = key
                        if sumcost > 0 and num1 > 0 and num2 > 0 and num3 > 0:
                            if reactiontreeinstance.contains(name):
                                name = name + '_' + str(counter)
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))


                                continue
                            else:
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))

                        else:
                            continue
                    if i[0] == 'addH2':

                        m = re.match(r'([A-Z])([0-9]+)([A-Z])([0-9]+)([A-Z])([0-9]+)', previousname)
                        let1 = m.group(1)
                        let2 = m.group(3)
                        let3 = m.group(5)
                        if let1 == 'C':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'C':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'C':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'O':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'O':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'O':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'H':
                            name =''
                            num1 = int(m.group(2)) + 2
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'H':
                            name =''
                            num2 = int(m.group(4)) + 2
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'H':
                            name =''
                            num3 = int(m.group(6)) + 2
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        cost = i[1]
                        sumcost = round((float(i[1])) + float(leaf.data.sum), 6)
                        path = leaf.data.path
                        path.append(i[0])
                        treename = leaf.data.treename
                        key = "".join(nameordered)
                        name = key
                        if sumcost > 0 and num1 > 0 and num2 > 0 and num3 > 0:
                            if reactiontreeinstance.contains(name):
                                name = name + '_' + str(counter)
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))


                                continue
                            else:

                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))

                        else:
                            continue
                    if i[0] == 'addCO':


                        m = re.match(r'([A-Z])([0-9]+)([A-Z])([0-9]+)([A-Z])([0-9]+)', previousname)

                        let1 = m.group(1)
                        let2 = m.group(3)
                        let3 = m.group(5)
                        if let1 == 'C':
                            name =''
                            num1 = int(m.group(2)) + 1
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'C':
                            name =''
                            num2 = int(m.group(4))  + 1
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'C':
                            name =''
                            num3 = int(m.group(6)) + 1
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'O':
                            name =''
                            num1 = int(m.group(2)) + 1
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'O':
                            name =''
                            num2 = int(m.group(4)) + 1
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'O':
                            name =''
                            num3 = int(m.group(6)) + 1
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'H':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'H':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'H':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        name = let1 +  str(num1) + let2 + str(num2) + let3 + str(num3)
                        cost = i[1]
                        sumcost = round((float(i[1])) + float(leaf.data.sum), 6)
                        path = leaf.data.path
                        path.append(i[0])
                        treename = leaf.data.treename
                        key = "".join(nameordered)
                        name = key
                        if sumcost > 0 and num1 > 0 and num2 > 0 and num3 > 0:
                            if reactiontreeinstance.contains(name):
                                name = name + '_' + str(counter)
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))


                                continue
                            else:
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))

                        else:
                            continue
                    if i[0] == 'addCH2':

                        m = re.match(r'([A-Z])([0-9]+)([A-Z])([0-9]+)([A-Z])([0-9]+)', previousname)

                        let1 = m.group(1)
                        let2 = m.group(3)
                        let3 = m.group(5)
                        if let1 == 'C':
                            name =''
                            num1 = int(m.group(2)) + 1
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'C':
                            name =''
                            num2 = int(m.group(4)) + 1
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'C':
                            name =''
                            num3 = int(m.group(6)) + 1
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'O':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'O':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'O':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'H':
                            name =''
                            num1 = int(m.group(2)) + 2
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'H':
                            name =''
                            num2 = int(m.group(4)) + 2
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'H':
                            name =''
                            num3 = int(m.group(6)) + 2
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        name = let1 +  str(num1) + let2 + str(num2) + let3 + str(num3)
                        cost = i[1]
                        sumcost = round((float(i[1])) + float(leaf.data.sum), 6)
                        path = leaf.data.path
                        path.append(i[0])
                        treename = leaf.data.treename
                        key = "".join(nameordered)
                        name = key
                        if sumcost > 0 and num1 > 0 and num2 > 0 and num3 > 0:
                            if reactiontreeinstance.contains(name):
                                name = name + '_' + str(counter)
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))


                                continue
                            else:
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))

                        else:
                            continue
                    if i[0] == 'subCO2':


                        m = re.match(r'([A-Z])([0-9]+)([A-Z])([0-9]+)([A-Z])([0-9]+)', previousname)

                        let1 = m.group(1)
                        let2 = m.group(3)
                        let3 = m.group(5)
                        if let1 == 'C':
                            name =''
                            num1 = int(m.group(2)) - 1
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'C':
                            name =''
                            num2 = int(m.group(4)) - 1
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'C':
                            name =''
                            num3 = int(m.group(6)) - 1
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'O':
                            name =''
                            num1 = int(m.group(2)) - 2
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'O':
                            name =''
                            num2 = int(m.group(4)) - 2
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'O':
                            name =''
                            num3 = int(m.group(6)) - 2
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'H':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'H':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'H':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name

                        cost = i[1]
                        sumcost = round((float(i[1])) + float(leaf.data.sum), 6)
                        path = leaf.data.path
                        path.append(i[0])
                        treename = leaf.data.treename
                        key = "".join(nameordered)
                        name = key
                        if sumcost > 0 and num1 > 0 and num2 > 0 and num3 > 0:
                            if reactiontreeinstance.contains(name):
                                name = name + '_' + str(counter)
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))


                                continue
                            else:
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))

                        else:
                            continue
                    if i[0] == 'subH2O':


                        m = re.match(r'([A-Z])([0-9]+)([A-Z])([0-9]+)([A-Z])([0-9]+)', previousname)

                        let1 = m.group(1)
                        let2 = m.group(3)
                        let3 = m.group(5)
                        if let1 == 'C':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'C':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'C':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'O':
                            name =''
                            num1 = int(m.group(2)) - 1
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'O':
                            name =''
                            num2 = int(m.group(4)) - 1
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'O':
                            name =''
                            num3 = int(m.group(6)) - 1
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'H':
                            name =''
                            num1 = int(m.group(2)) - 2
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'H':
                            name =''
                            num2 = int(m.group(4)) - 2
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'H':
                            name =''
                            num3 = int(m.group(6)) - 2
                            name += let3
                            name += str(num3)
                            nameordered[2] = name

                        cost = i[1]
                        sumcost = round((float(i[1])) + float(leaf.data.sum), 6)
                        path = leaf.data.path
                        path.append(i[0])
                        treename = leaf.data.treename
                        key = "".join(nameordered)
                        name = key
                        if sumcost > 0 and num1 > 0 and num2 > 0 and num3 > 0:
                            if reactiontreeinstance.contains(name):
                                name = name + '_' + str(counter)
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))


                                continue
                            else:
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))

                        else:
                            continue
                    if i[0] == 'subO':

                        m = re.match(r'([A-Z])([0-9]+)([A-Z])([0-9]+)([A-Z])([0-9]+)', previousname)
                        let1 = m.group(1)
                        let2 = m.group(3)
                        let3 = m.group(5)
                        if let1 == 'C':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'C':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'C':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'O':
                            name =''
                            num1 = int(m.group(2)) - 1
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'O':
                            name =''
                            num2 = int(m.group(4)) - 1
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'O':
                            name =''
                            num3 = int(m.group(6)) - 1
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'H':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'H':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'H':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name

                        cost = i[1]
                        sumcost = round((float(i[1])) + float(leaf.data.sum), 6)
                        path = leaf.data.path
                        path.append(i[0])
                        treename = leaf.data.treename
                        key = "".join(nameordered)
                        name = key
                        if sumcost > 0 and num1 > 0 and num2 > 0 and num3 > 0:

                            if reactiontreeinstance.contains(name):
                                name = name + '_' + str(counter)
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))

                                continue

                                # try:
                                #     node = reactiontreeinstance.get_node(name)
                                #     node.update_fpointer(previousname)
                                # except:
                                #     node = reactiontreeinstance.get_node(name)
                                #     node.update_fpointer('root')

                            else:
                                    reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))


                        else:
                            continue
                    if i[0] == 'subH2':
                        m = re.match(r'([A-Z])([0-9]+)([A-Z])([0-9]+)([A-Z])([0-9]+)', previousname)

                        let1 = m.group(1)
                        let2 = m.group(3)
                        let3 = m.group(5)
                        if let1 == 'C':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'C':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'C':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'O':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'O':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'O':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'H':
                            name =''
                            num1 = int(m.group(2)) -2
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'H':
                            name =''
                            num2 = int(m.group(4)) - 2
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'H':
                            name =''
                            num3 = int(m.group(6)) - 2
                            name += let3
                            name += str(num3)
                            nameordered[2] = name

                        cost = i[1]
                        sumcost = round((float(i[1])) + float(leaf.data.sum), 6)
                        path = leaf.data.path
                        path.append(i[0])
                        treename = leaf.data.treename
                        key = "".join(nameordered)
                        name = key
                        if sumcost > 0 and num1 > 0 and num2 > 0 and num3 > 0:
                            if reactiontreeinstance.contains(name):
                                name = name + '_' + str(counter)
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))

                                continue
                            else:
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))

                        else:
                            continue
                    if i[0] == 'subCO':
                            # previousname = leaf.data.reactiontype

                            m = re.match(r'([A-Z])([0-9]+)([A-Z])([0-9]+)([A-Z])([0-9]+)', previousname)

                            let1 = m.group(1)
                            let2 = m.group(3)
                            let3 = m.group(5)
                            if let1 == 'C':
                                name =''
                                num1 = int(m.group(2))  - 1
                                name += let1
                                name += str(num1)
                                nameordered[0] = name
                            elif let2 == 'C':
                                name =''
                                num2 = int(m.group(4)) - 1
                                name += let2
                                name += str(num2)
                                nameordered[1] = name
                            elif let3 == 'C':
                                name =''
                                num3 = int(m.group(6)) - 1
                                name += let3
                                name += str(num3)
                                nameordered[2] = name
                            if let1 == 'O':
                                name =''
                                num1 = int(m.group(2)) - 1
                                name += let1
                                name += str(num1)
                                nameordered[0] = name
                            elif let2 == 'O':
                                name =''
                                num2 = int(m.group(4)) - 1
                                name += let2
                                name += str(num2)
                                nameordered[1] = name
                            elif let3 == 'O':
                                name =''
                                num3 = int(m.group(6)) - 1
                                name += let3
                                name += str(num3)
                                nameordered[2] = name
                            if let1 == 'H':
                                name =''
                                num1 = int(m.group(2))
                                name += let1
                                name += str(num1)
                                nameordered[0] = name
                            elif let2 == 'H':
                                name =''
                                num2 = int(m.group(4))
                                name += let2
                                name += str(num2)
                                nameordered[1] = name
                            elif let3 == 'H':
                                name =''
                                num3 = int(m.group(6))
                                name += let3
                                name += str(num3)
                                nameordered[2] = name
                            cost = i[1]
                            sumcost = round((float(i[1])) + float(leaf.data.sum), 6)
                            path = leaf.data.path
                            path.append(i[0])
                            treename = leaf.data.treename
                            key = "".join(nameordered)
                            name = key
                            if sumcost > 0 and num1 > 0 and num2 > 0 and num3 > 0:
                                if reactiontreeinstance.contains(name):
                                    name = name + '_' + str(counter)
                                    reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))

                                    # try:
                                    #     node = reactiontreeinstance.get_node(name)
                                    #     node.update_fpointer(previousname)
                                    # except:
                                    #     node = reactiontreeinstance.get_node(name)
                                    #     node.update_fpointer('root')
                                    continue
                                else:
                                    reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))

                            else:
                                continue
                    if i[0] == 'subCH2':
                        # previousname = leaf.data.reactiontype
                        #
                        m = re.match(r'([A-Z])([0-9]+)([A-Z])([0-9]+)([A-Z])([0-9]+)', previousname)

                        let1 = m.group(1)
                        let2 = m.group(3)
                        let3 = m.group(5)
                        if let1 == 'C':
                            name =''
                            num1 = int(m.group(2)) - 1
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'C':
                            name =''
                            num2 = int(m.group(4)) - 1
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'C':
                            name =''
                            num3 = int(m.group(6)) - 1
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'O':
                            name =''
                            num1 = int(m.group(2))
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'O':
                            name =''
                            num2 = int(m.group(4))
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'O':
                            name =''
                            num3 = int(m.group(6))
                            name += let3
                            name += str(num3)
                            nameordered[2] = name
                        if let1 == 'H':
                            name =''
                            num1 = int(m.group(2)) - 2
                            name += let1
                            name += str(num1)
                            nameordered[0] = name
                        elif let2 == 'H':
                            name =''
                            num2 = int(m.group(4)) - 2
                            name += let2
                            name += str(num2)
                            nameordered[1] = name
                        elif let3 == 'H':
                            name =''
                            num3 = int(m.group(6)) - 2
                            name += let3
                            name += str(num3)
                            nameordered[2] = name

                        cost = i[1]
                        sumcost = round((float(i[1])) + float(leaf.data.sum), 6)
                        path = leaf.data.path
                        path.append(i[0])
                        treename = leaf.data.treename
                        key = "".join(nameordered)
                        name = key

                        if sumcost > 0 and num1 > 0 and num2 > 0 and num3 > 0:
                            if reactiontreeinstance.contains(name):
                                name = name + '_' + str(counter)
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))

                                # try:
                                #     node = reactiontreeinstance.get_node(name)
                                #     node.update_fpointer(previousname)
                                # except:
                                #     node = reactiontreeinstance.get_node(name)
                                #     node.update_fpointer('root')
                                continue
                            else:
                                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))

                        else:
                            continue
        previousmaxdepth = maxdepth

        return reactiontreeinstance, counter, previousmaxdepth

    elif maxdepth < maxlevel and maxdepth > maxleveltwo:
        nodelist = [('subH2', subH2), ('addH2', H2)]
        print('now adding only h')
        for leaf in leaves:
            parent = leaf.identifier
            previousname = leaf.data.reactiontype

            for i in nodelist:
                counter += 1
                name = ''
                nameordered = ['0','0','0']
                if i[0] == 'subH2':

                    m = re.match(r'([A-Z])([0-9]+)([A-Z])([0-9]+)([A-Z])([0-9]+)', previousname)

                    let1 = m.group(1)

                    let2 = m.group(3)
                    let3 = m.group(5)
                    if let1 == 'C':
                        name =''
                        num1 = int(m.group(2))
                        name += let1
                        name += str(num1)
                        nameordered[0] = name
                    elif let2 == 'C':
                        name =''
                        num2 = int(m.group(4))
                        name += let2
                        name += str(num2)
                        nameordered[1] = name
                    elif let3 == 'C':
                        name =''
                        num3 = int(m.group(6))
                        name += let3
                        name += str(num3)
                        nameordered[2] = name
                    if let1 == 'O':
                        name =''
                        num1 = int(m.group(2))
                        name += let1
                        name += str(num1)
                        nameordered[0] = name
                    elif let2 == 'O':
                        name =''
                        num2 = int(m.group(4))
                        name += let2
                        name += str(num2)
                        nameordered[1] = name
                    elif let3 == 'O':

                        name =''
                        num3 = int(m.group(6))
                        name += let3
                        name += str(num3)
                        nameordered[2] = name
                    if let1 == 'H':
                        name =''
                        num1 = int(m.group(2)) - 2
                        name += let1
                        name += str(num1)
                        nameordered[0] = name
                    elif let2 == 'H':
                        name =''
                        num2 = int(m.group(4)) - 2
                        name += let2
                        name += str(num2)
                        nameordered[1] = name
                    elif let3 == 'H':
                        name =''
                        num3 = int(m.group(6)) - 2
                        name += let3
                        name += str(num3)
                        nameordered[2] = name


                    cost = i[1]
                    sumcost = round((float(i[1])) + float(leaf.data.sum), 6)
                    path = leaf.data.path
                    path.append(i[0])
                    treename = leaf.data.treename
                    key = "".join(nameordered)
                    name = key

                    if sumcost > 0 and num1 > 0 and num2 > 0 and num3 > 0:
                        if not reactiontreeinstance.contains(name) :
                            reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))
                        else:
                            name = name + '_' + str(counter)
                            key = name
                            reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))


                    else:
                        continue
            if i[0] == 'addH2':

                m = re.match(r'([A-Z])([0-9]+)([A-Z])([0-9]+)([A-Z])([0-9]+)', previousname)

                let1 = m.group(1)
                let2 = m.group(3)
                let3 = m.group(5)
                if let1 == 'C':
                    name =''
                    num1 = int(m.group(2))
                    name += let1
                    name += str(num1)
                    nameordered[0] = name
                elif let2 == 'C':

                    name =''
                    num2 = int(m.group(4))
                    name += let2
                    name += str(num2)
                    nameordered[1] = name
                elif let3 == 'C':
                    name =''
                    num3 = int(m.group(6))
                    name += let3
                    name += str(num3)
                    nameordered[2] = name
                if let1 == 'O':
                    name =''
                    num1 = int(m.group(2))
                    name += let1
                    name += str(num1)
                    nameordered[0] = name
                elif let2 == 'O':
                    name =''
                    num2 = int(m.group(4))
                    name += let2
                    name += str(num2)
                    nameordered[1] = name
                elif let3 == 'O':
                    name =''
                    num3 = int(m.group(6))
                    name += let3
                    name += str(num3)
                    nameordered[2] = name
                if let1 == 'H':
                    name =''
                    num1 = int(m.group(2)) + 2
                    name += let1
                    name += str(num1)
                    nameordered[0] = name
                elif let2 == 'H':
                    name =''
                    num2 = int(m.group(4)) + 2
                    name += let2
                    name += str(num2)
                    nameordered[1] = name
                elif let3 == 'H':
                    name =''
                    num3 = int(m.group(6)) + 2
                    name += let3
                    name += str(num3)
                    nameordered[2] = name

                cost = i[1]
                sumcost = round((float(i[1])) + float(leaf.data.sum), 6)
                path = leaf.data.path
                path.append(i[0])
                treename = leaf.data.treename
                key = "".join(nameordered)
                name = key
                print('added only H')
                break
                if sumcost > 0 and num1 > 0 and num2 > 0 and num3 > 0:
                    if not reactiontreeinstance.contains(name):
                        reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))
                    else:

                        name = name + '_' + str(counter)
                        key = name
                        reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))

                else:
                    continue
        previousmaxdepth = maxdepth

        return reactiontreeinstance, counter, previousmaxdepth

    else:
        print('saving tree to working directory...')
        # saves file with initiatorname that served as root and
        filename = 'startInitiatorname' + initiatorname[0] + 'maxlevel' + str(maxlevel) + '.txt'
        dicttree = reactiontreeinstance.to_dict()
        print('converted tree to dict...')
        visualdict = toVisualDict(reactiontreeinstance)
        print('converted tree to visualizationdict...')
        dictfilename = 'DICTIOARNYTREE' + filename.rstrip('.txt') + '.json'
        json.dump(visualdict, open(dictfilename,'w'))
        print('converted to json...')
        reactiontreeinstance.save2file(filename)
        print('saved to file...')
        return None, counter, previousmaxdepth




def toVisualDict(reactiontreeinstance, nid=None, key=None, reverse=False):
        """transform self into a json serializable object for visualization"""


        nid = reactiontreeinstance.root if (nid is None) else nid
        # , "data":[reactiontreeinstance[nid].data.sum, reactiontreeinstance[nid].data.path]
        tree_dict = {"name": reactiontreeinstance[nid].tag , "mass": reactiontreeinstance[nid].data.sum, "reactiontype": reactiontreeinstance[nid].data.reactiontype, "children": []}
        if reactiontreeinstance[nid].expanded:
            queue = [reactiontreeinstance[i] for i in reactiontreeinstance[nid].fpointer]
            key = (lambda x: x) if (key is None) else key
            queue.sort(key=key, reverse=reverse)
            # iterate through all nodes
            for elem in queue:
                identifier = elem.identifier
                # recursively append children
                tree_dict["children"].append(toVisualDict(reactiontreeinstance, identifier))
            if tree_dict["children"] == []:
                tree_dict = {"name": reactiontreeinstance[nid].tag, "mass": reactiontreeinstance[nid].data.sum}


            return tree_dict


# cerates next level in a single reaction tree
def createnextlevelmassdiff(massdifference, nodelist, reactiontreeinstance, parent, counter):
    reactiontrees = []
    leaves = reactiontreeinstance.leaves('root')
    for leaf in leaves:
        parent = leaf.identifier
        counter += 1
        for i in nodelist:

            cost = i[1]
            # print(cost)
            # total cost of reactions so far
            # round to 6 digits
            sumcost = round(float(i[1]) + float(leaf.data.sum), 6)
            massdifference = round(massdifference, 6)
            # answer = str(round(answer, 2))
            # reactions so far
            path = leaf.data.path
            path.append(i[0])
            # print(path)
            treename = leaf.data.treename
            key = i[0] + '_'+ str(counter)
            name = i[0] + '_'+ str(counter)
            # create node only if reaction won't result in neg total cost
            # if we have found a pathway to a massdifference return the treeinstance
            if sumcost == massdifference:
                #TODO get rootnode value for path first element
                rootnode = reactiontreeinstance.get_node('root')
                path = rootnode.data.path
                firstelementpath = path[0]
                # print('sumcost == massdiff')
                path.insert(0, firstelementpath)
                return path
            elif sumcost > 0 and num1 > 0 and num2 > 0 and num3 > 0:
                reactiontreeinstance.create_node(key, name, parent=parent,data=ReactionNode(cost, path, i[0], sumcost, treename))

            elif sumcost != massdifference and counter > 500:
                return 'no tree found for this pattern after 20 reactionsteps'

    # every level return reactiontreeinstances created for next iteration
    # print(reactiontreeinstance)
    # print(counter)
    return (reactiontreeinstance, counter)
