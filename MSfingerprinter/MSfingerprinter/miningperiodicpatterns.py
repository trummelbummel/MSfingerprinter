from collections import defaultdict
import numpy as np
from bitarray import bitarray
import re

'''
script which mines frequent one cycles
'''

# F1 is also generally kept and only further iterations have join and prune steps involved
def findOneCyclePatterns(Tslice, minconfidence, periodicityhintarray, period_j):
    '''finds frequent 1-cycles in the working cube W[Time, Period, Value] in Tslice
    time index corresponds to offset, size of time index corresponds to period p,
    size of period index is total number of periods occuring in the time series'''
    # event is frequent if it occurs no less frequent than minconf * max(periodindex)
    # scanning once the periodindexall summary slice finds patterns that are more frequent
    # than the acceptable minfrequency
    # F1 set
    frequentonecyclepatterns = []
    frequentonecyclepatternscounts = []
    F1 = []
    F1counts = []
    periodlength = periodicityhintarray[period_j]
    # get maximum number of periods from ordered dict
    maxperiods = int(Tslice['period_indizes'].keys()[-1])
    valuelist = list(Tslice['periodindexall'].values())
    keylist = list(Tslice['periodindexall'].keys())
    for k,v in Tslice['periodindexall'].items():

        for i in range(len(v)):
            support = v[i]
            # maxperiods is: length of TS / periodlength
            minfrequency = minconfidence * float(maxperiods)
            binval = re.match(r'^[0-9]', k).group(0)
            if support >= minfrequency and support != 0:
                # C= (every nth position, starting at offset o, repeats value i) gets appended
                # C =(p,o.V) p = length or period of the cycle, o is offset indicating the first time at which cycle occurs, and V values that form cycle
                # is a concept category thus output of algorithm will be: a set containing patterns like this
                # P0 = (3,1,{(0,1)}) first position every nth position, 1-cyclepattern, dict key of cube where value is found
                # binposition indicates where in the cube value is located
                binposition = k
                binmatch1 = re.search(r'(^[0-9])', binposition).group(0)
                binmatch = (float(i) , float(binmatch1))
                # i = periodlength, 1 indicates it's a 1 cycle pattern
                # frozen set because set can only contain immutable datatype
                frequentonecyclepatterns.append((periodlength, 1, binmatch))
                frequentonecyclepatternscounts.append((support, periodlength, 1, binmatch))
        F1 = set(frequentonecyclepatterns)
        F1counts = set(frequentonecyclepatternscounts)
    return F1, F1counts
