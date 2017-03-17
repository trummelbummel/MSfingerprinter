import pandas as pd
import numpy as np
import os
from sklearn.cluster import AffinityPropagation
from sklearn import metrics
import re
import math
# scikitlearn
from sklearn.metrics.pairwise import euclidean_distances
import scipy.spatial.distance
from scipy import sparse
# plotting
from matplotlib import pyplot as plt
from itertools import cycle



#########################################3
#
# Entropy based feature selection
#
######################################3

def getAlpha(X):
    numpoints = len(X)
    Xtrans = X.T
    completedist = euclidean_distances(Xtrans)
    distmat = np.matrix(completedist)
    sumdist = distmat.sum()
    averagedistance = sumdist/(numpoints*numpoints)
    alpha = (-math.log(0.5)/averagedistance)
    return alpha, distmat

# results in mostly 0 similarity for higher values
def getSimilarities(alpha, completedistances):
    similarities = []
    lineararraydist = completedistances.ravel()
    for i in np.nditer(lineararraydist):
        distance = i
        Sij = math.exp(-alpha*distance)
        similarities.append(Sij)
    return similarities

def computetotalEntropyData(X):
    alpha, completedistances = getAlpha(X)
    similarities = getSimilarities(alpha, completedistances)
    arraysimilarities = np.array(similarities)
    # mask because of minus infinity
    logsimilarities = np.ma.log(arraysimilarities)
    ones = np.ones(arraysimilarities.size)
    oneminussim = ones-arraysimilarities
    # mask because of minus infinity
    logoneminussim = np.ma.log(oneminussim)
    Entropy = - np.sum(arraysimilarities * logsimilarities + oneminussim * logoneminussim)
    return Entropy

# Citation: feature selection for clustering, Dash
def doRankfeatures(X, featurenames, rangestart, rangeend):
    totalentropy = computetotalEntropyData(X)
    featureimportancearray = []
    indifferencefeatures = []
    originalX = X
    # set minimum to infinity
    minimumentropy = np.inf
    for i in range(len(featurenames)):
        Xtoremove = originalX
        Xnew = np.delete(Xtoremove, i, axis=0)
        X = Xnew.T
        totalentropynew = computetotalEntropyData(X)
        # if new entropy without feature is new minimum etropy
        # then feature is least important
        if totalentropynew < minimumentropy:
            featureimportancearray.append((featurenames[i], totalentropynew))
            minimumentropy = totalentropynew
            X = originalX
    return featureimportancearray


def statsentropy(rankedentropyvector, type):
    rankedentropyvectorentropies = rankedentropyvector
    statistics = []
    histogram = np.histogram(rankedentropyvectorentropies)
    plt.hist(rankedentropyvectorentropies, bins=20)  # plt.hist passes it's arguments to np.histogram
    plt.title("Histogram of entropyvector in "+ type +" space")
    plt.show()
    mean = np.mean(rankedentropyvectorentropies)
    statistics.append(['mean', mean])
    variance = np.var(rankedentropyvectorentropies)
    statistics.append(['variance', variance])
