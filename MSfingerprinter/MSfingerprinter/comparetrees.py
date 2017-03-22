import json
from treelib import Node, Tree
import os
import fnmatch

class ReactionNode(dict):
    def __init__(self, mass, reactiontype, name):
        self.mass = mass # contains reactionpath so far
        self.reactiontype = reactiontype # contains reaction type as string
        self.name = name
        self.children = []

def getdata(filename):
    with open(filename) as f:
        data = json.load(f)
    return data

# recursively build tree from file, since we want to search for masses we now have the mass as nodeID and chemical substance
# as human readable tag
def retrievenodes(data, counterpassed,  parent, child, reactiontreeinstance):
    counter = counterpassed
    name = round(data['mass'],6)
    mass = data['name']
    # cheating to get original mass in the mix
    reactiontype = data['mass']

    # nodeID, human readable tag
    reactiontreeinstance.create_node(mass, name, parent = parent, data=ReactionNode(mass, reactiontype, name))
    parent = name
    for k, v in data.iteritems():
        counter += 1
        if counter >= 1:
            if isinstance(v, list):
                for i in v:
                    data = i
                    try:
                        children = data['children']
                        if children:
                            for j in children:
                                mass = data['mass']
                                name = data['name']
                                reactiontype = data['reactiontype']
                                retrievenodes(data, counter, parent, j, reactiontreeinstance)
                    except:
                        continue
    return reactiontreeinstance

# search subtree and return json file to visualize sub
def searchMasspattern(reactiontreeinstance, rootsubtree, endsubtree, pattern, filename):
    nodes = []
    paths = []
    patternandtreeroots = []
    rootnodeid = reactiontreeinstance.root
    # gets initiatorname
    rootnode = reactiontreeinstance.get_node(rootnodeid).data.mass
    start = rootsubtree
    target = endsubtree
    node = getNode(reactiontreeinstance, start)
    if node != None:
        # skip firstnode as it is start node of subtree
        # return soft copy of the subtree rooted at start to search for paths in this tree
        # means all nodes are shared between original and soft copy
        subtree = reactiontreeinstance.subtree(start)
        # check if second chemical entity involved in pattern is in subtree
        if  subtree.contains(target):
            print('contains path')
            print(start,target)
            visualdict = toVisualDict(subtree)
            resultpath = os.getcwd + '/results/resultsubtrees/'
            dictfilename = 'DICTIOARNYTREEPATH' + 'from' + str(start) + 'to' + str(target) + str(filename)+ '.json'
            savepath = os.path.join(resultpath, dictfilename)
            json.dump(visualdict, open(savepath,'w'))
            return rootnode, pattern

        else:

            print('node found in Initiatortree' + str(filename))
            print(start)
            print(target)
            visualdict = toVisualDict(subtree)
            resultpath = os.getcwd + '/results/resultsubtrees/'
            dictfilename = 'DICTIOARNYTREENODE' + str(start) + 'from' + str(target) + 'to' + str(filename)+  '.json'
            savepath = os.path.join(resultpath, dictfilename)
            json.dump(visualdict, open(savepath,'w'))
            return rootnode, pattern
    else:
        return  None, None

def searchMasspatterngroundtruth(reactiontreeinstance, rootsubtree, endsubtree, pattern, filename):
    nodes = []
    paths = []
    patternandtreeroots = []
    rootnodeid = reactiontreeinstance.root
    # gets initiatorname
    rootnode = reactiontreeinstance.get_node(rootnodeid).data.mass

    start = rootsubtree
    target = endsubtree
    # if reactiontreeinstance.
    node = getNode(reactiontreeinstance, target)
    if node != None:

        # skip firstnode as it is start node of subtree
        # return soft copy of the subtree rooted at start to search for paths in this tree
        # means all nodes are shared between original and soft copy
        subtree = reactiontreeinstance.subtree(target)
        # check if second chemical entity involved in pattern is in subtree
        if  subtree.contains(start):
            print('contains path')
            print(start,target)
            print('stoichiometformula')
            stoichiometformula = subtree.get_node(target).data.mass
            visualdict = toVisualDict(subtree)
            resultpath = os.getcwd + '/results/resultssubtrees/'
            dictfilename = 'DICTIOARNYTREEPATH' + 'from' + str(start) + 'to' + str(target) + str(filename)+ '.json'
            savepath = os.path.join(resultpath, dictfilename)
            json.dump(visualdict, open(savepath,'w'))
            return rootnode, stoichiometformula

        else:

            print('node found in Initiatortree' + str(filename))
            print(start)
            print(target)
            stoichiometformula = subtree.get_node(target).data.mass
            visualdict = toVisualDict(subtree)
            resultpath = os.getcwd + '/results/resultssubtrees/'
            dictfilename = 'DICTIOARNYTREENODE' + str(start) + 'from' + str(target) + 'to' + str(filename)+  '.json'
            savepath = os.path.join(resultpath, dictfilename)
            json.dump(visualdict, open(savepath,'w'))
            return rootnode, stoichiometformula
    else:
        return  None, None




def toVisualDict(reactiontreeinstance, nid=None, key=None, reverse=False):
        """transform self into a json serializable object for visualization"""
        nid = reactiontreeinstance.root if (nid is None) else nid
        # , "data":[reactiontreeinstance[nid].data.sum, reactiontreeinstance[nid].data.path]
        tree_dict = {"name": reactiontreeinstance[nid].data.mass , "mass": reactiontreeinstance[nid].data.name, "reactiontype": reactiontreeinstance[nid].data.reactiontype, "children": []}
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
                tree_dict = {"name": reactiontreeinstance[nid].data.mass, "mass": reactiontreeinstance[nid].data.name}

            return tree_dict



# find files path, reads csv files only unless specified differently in extensions
def find_files(path, extensions):
    # Allow json
    extensions = [e.replace(".", "") for e in extensions]
    for dirpath, dirnames, files in os.walk(path):
        for extension in extensions:
            for f in fnmatch.filter(files, "*.%s" % extension):
                p = os.path.join(dirpath, f)
                yield (p)

def getNode(reactiontreeinstance, mass):

    node = reactiontreeinstance.get_node(mass)
    return node
