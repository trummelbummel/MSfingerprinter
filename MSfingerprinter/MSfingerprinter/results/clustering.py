import pandas as pd
import numpy as np
from sklearn.cluster import AffinityPropagation
from sklearn import metrics
import re
# Clustering
from sklearn.cluster import DBSCAN
import sklearn.metrics.pairwise
import scipy.spatial.distance
from sklearn.cluster import MiniBatchKMeans, KMeans
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.neighbors import NearestNeighbors
from scipy import sparse
# plotting
from matplotlib import pyplot as plt
# hierachical Clustering
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist
from sklearn import preprocessing

numberofsamples = 9
# creates similarity matrix from DTW output

def readinputclustering(filename):
    df = pd.read_csv(filename, header=None)
    X = df.ix[:, 1::]
    print('nan value count')
    print(X.isnull().sum())
    X.fillna(0, inplace=True)
    min_max_scaler = preprocessing.MinMaxScaler()
    x_scaled = min_max_scaler.fit_transform(X)
    X = pd.DataFrame(x_scaled)
    # maxlength = len(X.columns)
    # print(maxlength)
    labels = df.ix[:, 0]
    # print(X)
    # print(labels)
    return X, labels



# # log transform the dataframe cause values differ by orders of magnitude
# X = np.log(X)
# X[~np.isfinite(X)] = 0
# if dataframe rows (fingerprint freqs) are unequal length fill with zeros



#
# # ##################################################3
# # #
# # # Parameter estimation
# # #
# # ##################################################
# #
# # DBSCAN parameters
# # look for at the distance distribution of dataset which only contains the
# # parameters we actually cluster over i.e. the penalties and set eps
# # to be a sensible value based on this analysis. For eps to be appropriate
# # we should pick a value for the distance of all points to its nearest neighbors
# # that includes most of the dataset such that one has a small number of outliers.
# #
# neighbors = NearestNeighbors(n_neighbors=2, algorithm='ball_tree').fit(X)
# distances, indices = neighbors.kneighbors(X)
# stepsize = 0.1
# histo = np.histogram(distances, bins=np.arange(0, 10, stepsize))
# # cutoff in percent
# cutoff = 1
# # look if change greater some percent value cutoff
# def percentagechange(histogram):
#     percentchange = []
#     cumsum = np.cumsum(histo[0])
#     sumtotal = cumsum[-1]
#     for i in range(len(histogram)):
#         change = abs((float(histogram[i-1] - histogram[i]) /float(sumtotal)) * 100)
#         percentchange.append(change)
#     return percentchange
#
# # looks if cutoff smaller
# def geteps(percentchange, cutoff):
#     for i in range(len(percentchange)):
#         if percentchange[i] < cutoff:
#             return i
#
# percentchange = percentagechange(histo[0])
#
# eps = stepsize*geteps(percentchange, cutoff)
#
# numberofdatapoints = len(df)
# numtry = 10
# epsparam = geteps(percentchange, 0.01)
# min_samplesparam = 3
# # creates an array for clustersim
# minptsarray = np.arange(min_samplesparam, 30, 2)
#
#
# # as an example plot the histogram
# print(epsparam)
#
# # visualize heuristic for choosing eps
# plt.hist(distances, bins=np.arange(0,3, stepsize))
# plt.xlabel('Distances')
# plt.ylabel('Frequency')
# plt.show()

############################################33
#
# DBScan
#
######################################################33
#
# min_samplesparam = 3
#
# db = DBSCAN(eps=epsparam, min_samples=min_samplesparam).fit(X)
# # extract clusterlabels and outliers
# core_samples = db.core_sample_indices_
# # define a boolean mask so we can plot the core, boundary, and noise points separately
# core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
# core_samples_mask[db.core_sample_indices_] = True
# labels = db.labels_
# n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
#
#
# print("Silhouette Coefficient: %0.3f"
#       % metrics.silhouette_score(X, labels))
#
# # 2D plot using overlaps and stretches on x and y axis
# plt.scatter(X.iloc[:, 0], X.iloc[:, 1], c=labels)
# plt.title("DBScan Clustering")
# plt.show()
#
# # returns clusters and outliers as lists
# clusters = [X[labels == i] for i in xrange(n_clusters_)]
# outliers = [X[labels == -1]]
#
##############################################
#
# K-MEANS
###########################################

# returns all datapoints in a cluster when given the clusternumber and the labels
def ClusterIndicesNumpy(clustNum, labels_array): #numpy
    return np.where(labels_array == clustNum)[0]


# PCA reduced K means
def runKmeans(X, numcluster):
    reduced_data = PCA(n_components=numcluster).fit_transform(X)
    kmeans = KMeans(init='k-means++', n_clusters=numcluster, n_init=20)
    kmeans.fit(reduced_data)
    labels = kmeans.labels_
    n_clusters_ = len(set(labels))
    for i in range(len(n_clusters)):
        'clusters' + str(i) +': ' + str(ClusterIndicesNumpy(i, labels))
    # returns clusters and outliers as lists
    print('Estimated number of clusters: %d' % n_clusters_)
    # plot kmeans results
    plt.scatter(reduced_data[:, 0], reduced_data[:, 1], c=labels)
    plt.title("Kmeans with K" + str(numcluster))
    plt.show()


#
#
# kmeans = KMeans(init='k-means++', n_clusters=2, n_init=20)
# kmeans.fit(X)
# labels = kmeans.labels_
# n_clusters_ = len(set(labels))
#
# # returns all datapoints in a cluster when given the clusternumber and the labels
# def ClusterIndicesNumpy(clustNum, labels_array): #numpy
#     return np.where(labels_array == clustNum)[0]
#
# clusters1 = ClusterIndicesNumpy(1, labels)
# clusters2 = ClusterIndicesNumpy(2, labels)
#
#
# print('cluster1')
# print(clusters1)
# print('cluster2')
# print(clusters2)
#
#
# # returns clusters and outliers as lists
# print('Estimated number of clusters: %d' % n_clusters_)
# #
# # # plot kmeans results
# # plt.scatter(X[:, 0], X[:, 1], c=labels)
# # plt.title("K = 2")
# # plt.show()


######################################################33
#
# Hierachical Clustering
#
######################################################


def plotDendrogram(linkage_mat, labels):
    titles = labels.tolist()
    print(labels)
    print(titles)
    # calculate full dendrogram
    plt.figure(figsize=(25, 10))
    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('sample index')
    plt.ylabel('distance')
    dendrogram(
        linkage_mat,
        p = 12, # show only last p merges
        # labels=titles,
        show_leaf_counts = False,
        leaf_rotation=90.,  # rotates the x axis labels
        leaf_font_size=8.,  # font size for the x axis labels
    plt.show()





def hierachicalclustering(X, labels, mode='ward'):
    X = np.array(X)
    if mode == 'single':
        linkage_mat = linkage(X, 'single')
        c, coph_dists = cophenet(linkage_mat, pdist(X))
        plotDendrogram(linkage_mat, labels)
    elif mode == 'average':
        linkage_mat = linkage(X, 'average')
        c, coph_dists = cophenet(linkage_mat, pdist(X))
        plotDendrogram(linkage_mat, labels)
    # defaults to ward merge strategy
    else:
        # generate linkage matrix
        linkage_mat = linkage(X, 'ward')
        # The closer the value of c is to 1, the better the clustering preserves the original distance
        # c is the copheny correlation coefficient
        c, coph_dists = cophenet(linkage_mat, pdist(X))
        plotDendrogram(linkage_mat, labels)
    return linkage_mat
