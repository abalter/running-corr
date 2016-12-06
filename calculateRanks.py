import numpy as np
from collections import defaultdict as dd
from collections import Counter as ctr
from collections import OrderedDict as od
from math import ceil

#a = (np.random.geometric(p=0.01, size=10**7)-1)*10
#np.savetxt("test_depths.txt", a, fmt="%d")

filename = "short_depths.txt"

depths = dd(lambda:0)

with open(filename, 'r') as f:
    for line in f:
        depth = int(line.split("\t")[0])
        depths[depth] += 1
        
a = [5,1,2,2,9,345,7,100,7,100,7]
depths = ctr()
for el in a:
    depths[el] += 1

# frequencies has:
# Keys: Depth values
# Values: Frequency of those depths
# Sorted by Values
frequencies = od(sorted(depths.items(), reverse=True))

rankdata = dict()
ranks = []
current_rank = 1
for depth, frequency in frequencies.items():
    #print("depth: {}, frequency: {}".format(frequency, count))
    if frequency > 1:
        N = current_rank - 1
        M = N + frequency
        sum_of_raw_ranks = M*(M+1)/2.0 - N*(N+1)/2.0
        average_rank = sum_of_raw_ranks*1.0/frequency
        rankdata[depth]['frequency'] = frequency
        rankdata[depth]['rank'] = average_rank
        current_rank = current_rank + frequency
    else:
        current_rank += 1
        ranks.append(current_rank)

        
    print("current rank: {}".format(current_rank))

    