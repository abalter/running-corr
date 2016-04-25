#!/usr/bin/python

import numpy as np
import sys

MEANS_RANGE = 100
COVARIANCES_RANGE = MEANS_RANGE**(3/2)

rows = 10000
cols = 5

covariances = np.zeros((cols, cols))
corr_coeffs = np.ones((cols, cols))
A = np.random.uniform(-1, 1, size=(cols, cols))*COVARIANCES_RANGE
covariances = A.dot(A.T)
means = np.random.uniform(-1, 1, cols)*MEANS_RANGE

for i in range(cols):
    for j in range(cols):
        corr_coeffs[i][j] = covariances[i][j]/np.sqrt(covariances[i][i]*covariances[j][j])
        

print("Corr coeffs")
print(corr_coeffs)
print("covariances")
print(covariances)
print()
print()



X = np.random.multivariate_normal(means, covariances, rows)

# print("X")
#print(X)
print()
print()

C = np.cov(X.T);
S = np.corrcoef(X.T)

print("cov")
print(C-covariances)
print("corr")
print(S-corr_coeffs)
