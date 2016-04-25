#!/usr/bin/python

import numpy as np
import running_corr as rc
from bisect import insort

# rc.getDepthCounts

MAX_DEPTH = 100

rows = 100
cols = 5

fake_depths = np.random.uniform(0, MAX_DEPTH, size=(rows, cols))

depths = np.zeros(size=(rows, cols))

for i in range(rows):
    



def getDepthCounts(data_filename):
    
    # get the number of sequences in the multipilup 
    # in order to initialize the correlation matrices
    num_columns = getNumColumns(data_filename)
            
    # initialize the values from a single line of the multipileup file
    values = zeros(num_columns)

    # initialize depth_counts with empty lists
    depth_counts = [{} for dummy in range(num_columns) ]
    max_depth = 0
    
    ### Read in the data
    
    # Read through the file and create a dictionary of reads.
    # Each column is a dictionary with keys that are the integer values
    # of reads, and the values are dictionaries with keys 'count' and 'rank'
    with open(data_filename, 'r') as data_file:
        for line in data_file: # read through data file line by line
            
            values = line.split() # split the depths in the line on whitespace]
            for i in range(values.length):
                read = int(values[i])
                # if that depth exists for that column then increment. Else set to 1.
                if read in depth_counts[i]:
                    depth_counts[i][read]['count'] += 1  # increment the number of times that depth appeared
                else:
                    depth_counts[i][read] = {'count':1, 'rank':0}
                
                # Track the maximum depth so that at the end can make a vector
                # with max_depth elements to hold the count of each depth
                max_depth = max(max_depth, max(values)) 
                
    ### Calculate the ranks
    
    # Fill in the 'rank' attribute of each read object.
    # For each column:
    #    1) Create an array of dictionaries with keys 'depth', 'count' and 'rank'
    #       that are populated from the column. Called 'ranked_depths'
    #    2) Sort the array on depths in ascending order
    #    3) Loop through sorted ranked_depths and assign consecutive ranks
    #       according to the fractional method that assigns ranks for tied
    #       counts equal to the average of their ordinal ranks.
    #    4) Finally, using the ranked_depths table, update the 'rank' field
    #       of the depth_count objects.
    
    for column in depth_counts:
        ranked_depths = [{'depth':key, 'count':val['count'], 'rank':0} for key, val in column]
        ranked_depths.sort(key=lambda x: x['depth'])
        
        count = 0
        next_rank = 1
        for i in range(ranked_depths.length):
            count = item['count']
            if count == 0:
                break
            else if count == 1:
                item['rank'] = next_rank
            else:
                rank = calculateTiedRank(current_rank, count)
                item['rank'] = rank
            
            next_rank += count
            
        for read in column:
            column[read]['rank'] = ranked_depths[read]
    
    
     
    
    return depth_counts


def writeTempRankFile(depth_counts):
    temp_file = open("temp", "w")
    array_of_ranks = []
    array_of_reads = []
    with open(data_filename, 'r') as data_file:
        for line in data_file: # read through data file line by line
            # Put the line of reads into an array
            array_of_reads = [int(read) for read in line.split()]
            # Translate to an array of ranks using depth_counts
            array_of_ranks = [depth_counts[read]['rank'] for read in array_of_reads]
            # Write a line of ranks
            temp_file.writeln(array_of_ranks)

            



            
def calculateTiedRank(current_rank, count):

"""
sum 1...N = N*(N-1)/2
sum N...M = M(M-1)/2 - (N-1)*(N-2)/2
"""

N = current_rank + 1
M = N + depth
sum_of_raw_ranks = M*(M-1)/2 - (N-1)*(N-2)/2
average_rank = sum_of_raw_ranks/depth
return average_rank



