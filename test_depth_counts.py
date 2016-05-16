#!/usr/bin/python

import numpy as np
import running_corr as rc
from bisect import insort
import pprint
pp = pprint.PrettyPrinter()

# rc.getDepthsAndRanks


MAX_DEPTH = 100

COUNTS_MEAN = 50
COUNTS_SIGMA = 20
NUM_RANKS = 100
COVARIANCES_RANGE = 3
MEANS_RANGE = 3

def main(rows, cols):
    rank_counts = np.random.poisson(1,NUM_RANKS).sort();


    ### Create a list of fake ranks
    #fake_counts = np.random.lognormvariate(COUNTS_MEAN, COUNT_SIGMA, )

    covariances = np.zeros((cols, cols))
    corr_coeffs = np.ones((cols, cols))
    correlators = np.random.uniform(-1, 1, size=(cols, cols))*COVARIANCES_RANGE
    print("correlators", correlators)
    covariances = correlators.dot(correlators.T)
    means = np.random.uniform(-1, 1, cols)*MEANS_RANGE

    for i in range(cols):
        for j in range(cols):
            corr_coeffs[i][j] = covariances[i][j]/np.sqrt(covariances[i][i]*covariances[j][j])
                                
    print("means")
    pp.pprint(means)

    print("covariances")
    pp.pprint(covariances)

    print("corr_coeffs")
    pp.pprint(corr_coeffs)

    writeTempReadFile(rows, cols, means, correlators);
    
    depths_counts_ranks = getDepthsAndRanks("temp_read_file")
    
    print("in main depths_counts_ranks")
    pp.pprint(depths_counts_ranks)
    
    writeTempRankFile("temp_read_file", depths_counts_ranks)



def writeTempReadFile(rows, cols, means, covariances):
    
    
    temp_file = open("temp_read_file", 'w')
    
    for i in range(rows):
        correlated_normals = np.random.multivariate_normal(means, covariances)
        #print("correlated normals")
        pp.pprint(correlated_normals)
        correlated_lognormals = [int(np.exp(x)) for x in correlated_normals]
        #print("correlated lognormals")
        pp.pprint(correlated_lognormals)
        # raw_input("")
        dummy = np.random.randint(1,10,cols)
        # line = ",".join( str(x) for x in correlated_lognormals)
        line = ",".join(str(x) for x in dummy)
        temp_file.write(line)
        temp_file.write("\n")
    
    temp_file.close()
    

def getCorrelatedLognormals(means, covariances, num_columns):
    #print("num columns", num_columns)
    normals = np.random.multivariate_normal(means, covariances)
    lognormals = np.exp(normals)
    correlated_lognormals = lognormals.dot(correlators) + means
    return correlated_lognormals   
    
    
### Calculate spearmans

### Create a list of data based on fake ranks

### Write data to file

### Test algorithm



def getDepthsAndRanks(data_filename):
    
    # get the number of sequences in the multipilup 
    # in order to initialize the correlation matrices
    num_columns = getNumColumns(data_filename)
            
    # initialize the values from a single line of the multipileup file
    values = np.zeros(shape=(num_columns,1))

    # initialize depth_counts with empty lists
    depth_counts = [{} for dummy in range(num_columns)]
    max_depth = 0
    
    ### Read in the data
    
    # Read through the file and create a dictionary of reads.
    # Each column is a dictionary with keys that are the integer values
    # of reads, and the values are dictionaries with keys 'count' and 'rank'
    with open(data_filename, 'r') as data_file:
        for line in data_file: # read through data file line by line
            
            values = line.split(",") # split the depths in the line on whitespace]
            #print("len(values)", len(values))
            #print("len(depth_counts", len(depth_counts))
            for i in range(len(values)):
                read = int(values[i])
                #print(i, "read", read)
                # if that depth exists for that column then increment. Else set to 1.
                if read in depth_counts[i]:
                    #print("depth exits, incrementing")
                    depth_counts[i][read]['count'] += 1  # increment the number of times that depth appeared
                    #print(depth_counts[i][read])
                    #print
                else:
                    #print("depth doesn't exist, initializing")
                    depth_counts[i][read] = {'count':1, 'rank':0}
                    #print
                
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
    
    print("depth counts")
    pp.pprint(depth_counts)
    raw_input("")
    
    for column in depth_counts:
        # Create a list of objects that have depth, count, rank for that depth.
        depths_counts_ranks = [{'depth':key, 'count':val['count'], 'rank':0} for key, val in column.iteritems()]
        # Sort the list by depth
        depths_counts_ranks.sort(key=lambda x: x['depth'])
        print("sorted depths_counts_ranks no ranks")
        pp.pprint(depths_counts_ranks)
        
        # Calculate the ranks using the fractional method for ties
        count = 0
        current_rank = 1
        for i in range(len(depths_counts_ranks)):
        #for item in depths_counts_ranks:
            item = depths_counts_ranks[i]
            #print("item[rank]", item['rank'])
            count = item['count']
            #print("count", count, "next_rank", current_rank)
            
            if count == 0:
                break
            elif count == 1:
                item['rank'] = current_rank
            else:
                tied_rank = calculateTiedRank(current_rank, count)
                #print("tied_rank", tied_rank)
                item['rank'] = tied_rank
                
            depths_counts_ranks[i] = item
            
            current_rank += count
            
            #print("ending current_rank", current_rank)
            #print("item", item)
            
        #print("column")
        #pp.pprint(column)
        
        print("depths_counts_ranks with ranks")
        pp.pprint(depths_counts_ranks)
        
        #print("ranked_depths.shape")
        #print(len(depths_counts_ranks))
        
        print("filling depth_counts array for current column")
        print("column")
        pp.pprint(column)
        for item in depths_counts_ranks:
            print("item")
            pp.pprint(item)
            column[item['depth']]['rank'] = item['rank']
            
    print("depth_counts")
    pp.pprint(depth_counts)
    
    
    return depth_counts


def writeTempRankFile(data_filename, depths_counts_ranks):
    print("writing temp rank file")
    print("depths_counts_ranks")
    pp.pprint(depths_counts_ranks)
    raw_input("")
    
    #temp_read_file = opent(data_filename, "r")
    temp_rank_file = open("temp_rank_file", "w")
    with open(data_filename, 'r') as data_file:
        for line in data_file: # read through data file line by line
            array_of_ranks = []
            print("line")
            print(line)
            # Put the line of reads into an array
            reads = line.split(",")
            int_reads = [int(read) for read in reads]
            print("int_reads")
            pp.pprint(int_reads)
            for i in range(len(int_reads)):
                rank = depths_counts_ranks[i][int_reads[i]]['rank']
                array_of_ranks.append(str(rank))
            print("array of ranks")
            pp.pprint(array_of_ranks)
            line_of_ranks = ",".join(array_of_ranks)
            print("line of ranks")
            print(line_of_ranks)
            # Write a line of ranks
            temp_rank_file.write(line_of_ranks + "\n")

            
def calculateTiedRank(current_rank, count):

    """
    sum 1...N = N*(N+1)/2
    sum N...M = M(M+1)/2 - N*(N+1)/2
    """

    N = current_rank - 1
    M = N + count
    #print("current rank", current_rank, "count", count)
    #print("N", N, "M", M)
    sum_of_raw_ranks = M*(M+1)/2.0 - N*(N+1)/2.0
    average_rank = sum_of_raw_ranks*1.0/count
    return average_rank
    

def getNumColumns(data_filename):

    data_file = open(data_filename, 'r')
    line = data_file.readline()
    data_file.close()

    return len(line.split(","))
    
    

if __name__ == "__main__":
    rows = 20
    cols = 2
    main(rows, cols)
