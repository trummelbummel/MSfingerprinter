# MSfingerprinter
Finds periodicities in mass spectral data based on the intensity domain and construction of reactiontrees describing mass differences in mass spectra

Implements:

Efficient mining of partial periodic patterns in time series database Han et al. 1999 (maximum subpattern hitset algorithm)

Multiple and partial periodicity mining in time series databases Berberedis et al. 2002

Software used:

peak_detection (MIT license 1.0.4) https://github.com/demotu/BMC

pysax https://github.com/dolaameng/pysax

d3js.org (BSD 3)

treelib (Apache License) http://treelib.readthedocs.io/en/latest/

Dependencies using Anaconda run:

pip install joblib


pip install treelib


Explanation of main functions:

Driver_entropy: Given inputfiles in rawdata/freatureranking with the format ['index', 'mass', 'freq', 'intensity']
this main function can identify the entropy in both mass and frequency space.

Driver_trees_allpossiblereactions: This lets one construct Initiatortrees in CHO space. Given Initiators i.e.
known sum formulas and masses and a massrange defined by minmass and maxmass one can cover part of the
search space adding and subtracting known repeat units.

Driver_periodicitybothspaces: Finds periodicities in Mass spectral data provided as a file of the format
['index', 'mass', 'freq', 'intensity'] alternatively ['index', 'mass', 'mass', 'intensity'] is also a possible input
if no frequency data is available for the mass spectrum. Space type lets one define if one wants to mine
patterns in either 'mass' (mass space) or 'freq' (frequency space). Output will be saved in files that can
then be processed using Driver_postprocessingpatterns.

Driver_postprocessingpatterns: Retrieves and saves mass differences that are chemically relevant in CHO space from periodic patterns
mined.

Driver_treefinder : finds masses that have been retrieved using Driver_postprocessing patterns
in Initiatortrees that are already constructed. Results saved in folder nodesfound.

Visualization of Initiatortrees possible replacing filename to respective tree
in index_initiatorTrees.html saved in subdirectory completeInitiatortrees
