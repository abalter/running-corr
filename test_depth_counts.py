#!/usr/bin/python

import numpy as np
import running_corr as rc

# rc.getDepthCounts

MAX_DEPTH = 100

rows = 100
cols = 5

fake_depths = np.random.uniform(0, MAX_DEPTH, size=(rows, cols))

depths = np.zeros(size=(rows, cols))

for i in range(rows):
    
