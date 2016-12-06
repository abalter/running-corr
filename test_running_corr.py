#!/bin/python

import numpy as np
import pprint
import running_corr as c


MEANS_RANGE = 100
COVARIANCES_RANGE = MEANS_RANGE**(3/2)

pp = pprint.PrettyPrinter()

def testRunningCorr(rows, cols):
    
    corr_coeffs, covariances, means = getCorrelationMatrix(cols)
    writeTempFileOfRVs(rows, cols, means, covariances)
    pcorr, pcov, means = c.pearsons("temp.txt")
    pcorr = np.matrix(pcorr)
    print("actual")
    pp.pprint(corr_coeffs.round(3))
    print("read")
    pp.pprint(pcorr.round(3))
    
    difference = pcorr - corr_coeffs
    print("difference")
    pp.pprint((pcorr - corr_coeffs).round(3))
    
    mean_square_test = 0
    
    test = difference*difference/rows**2
    
    print("mean square test: ", mean_square_test)
    

def getCorrelationMatrix(cols):
    
    # print("cols", cols)
    
    # Create a symmetric matrix with values
    # uniformly distributed in [-1,1]

    covariances = np.zeros((cols, cols))
    corr_coeffs = np.ones((cols, cols))
    A = np.random.uniform(-1, 1, size=(cols, cols))*COVARIANCES_RANGE
    covariances = A.dot(A.T)
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
    
    return corr_coeffs, covariances, means
    

def generateMeans(cols):
    
    means = np.random.uniform(-1, 1, cols)*MEANS_RANGE
    return means
    

def writeTempFileOfRVs(rows, cols, means, covariances):
    
    
    temp_file = open("temp.txt", 'w')
    
    for i in range(rows):
        correlated_normals = np.random.multivariate_normal(means, covariances)
        line = " ".join( [str(x) for x in correlated_normals] )
        temp_file.write(line)
        temp_file.write("\n")
    
    temp_file.close()
    

if __name__ == "__main__":
    
    print("__name__ == __main__")

    import sys
    
    args = sys.argv;
    
    print(args)
    
    rows = int(args[1])
    cols = int(args[2])
    
    print("rows", rows, "cols", cols)
    
    testRunningCorr(rows, cols)
    
else:
    getCorrelationMatrix(5)

# http://comisef.wikidot.com/tutorial:correlateduniformvariates
# http://stackoverflow.com/questions/32718752/how-to-generate-correlated-uniform0-1-variables
# http://math.stackexchange.com/questions/446093/generate-correlated-normal-random-variables
# 
