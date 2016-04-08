#!/bin/python

import numpy as np
import pprint

MEANS_RANGE = 100

def generateCorrelationMatrix(cols):
    cov = np.random.uniform(-1, 1, size=(cols, cols))
    corr_coeffs = np.zeros((cols, cols))
    for i in range(cols):
        for j in range(cols):
            corr_coeffs[i][j] = cov[i][j]/np.sqrt(cov[i][i]*cov[j][j])

    pp = pprint.PrettyPrinter()
    print("corr_coeffs")
    pp.pprint(corr_coeffs)
    
    return corr_coeffs, cov
    

def generateMeans(cols):
    
    means = np.random.uniform(-1, 1, size(cols, 1))*MEANS_RANGE
    return means
    

def getCorrelatedRandomNormals(cols):
    corr_coeffs, cov = generateCorrelationMatrix(cols)
    means = generateMeans(cols)
    
    correlated_normals = np.random.multivariate_normal(means, cov)
    
    return correlated_normals
    

def writeTempFileOfRVs(rows, cols):
    
    correlated_normals = getCorrelatedRandomNormals(cols)
    
    temp_file = open("temp", 'w')
    
    for i in range(rows):
        line = " ".join( [str(x) for x in correlated_normals] )
        temp_file.write(line)
    
    temp_file.close()
    

if __name__ == "__main__":

    import sys
    
    args = sys.argv;
    
    rows = args[1]
    cols = args[2]
    
    testRunningCorr(rows, cols)


# http://comisef.wikidot.com/tutorial:correlateduniformvariates
# http://stackoverflow.com/questions/32718752/how-to-generate-correlated-uniform0-1-variables
# http://math.stackexchange.com/questions/446093/generate-correlated-normal-random-variables
# 
