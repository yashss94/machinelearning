import numpy as np
from scipy.spatial import distance

B_temp = [[2.0, -2.0, -1.0, 3.0],
          [-1.0, 3.0, 3.0, -1.0],
          [0.0, 2.0, 3.0, 0.0],
          [1.0, 3.0, 1.0, 3.0],
          [1.0, 0.0, -1.0, 2.0],
          [-3.0, 2.0, 4.0, -1.0],
          [5.0, -1.0, 5.0, 3.0],
          [2.0, 1.0, 2.0, 0.0]]
A_temp = [[0.0, 0.0, 0.0, 0.0],
          [0.0, 0.0, 0.0, 0.0],
          [0.0, 0.0, 0.0, 0.0],
          [0.0, 0.0, 0.0, 0.0],
          [0.0, 0.0, 0.0, 0.0],
          [0.0, 0.0, 0.0, 0.0],
          [0.0, 0.0, 0.0, 0.0],
          [0.0, 0.0, 0.0, 0.0]]

scoringMatrix = [[0.0, 0.0, 0.0, 0.0],
                 [0.0, 0.0, 0.0, 0.0],
                 [0.0, 0.0, 0.0, 0.0],
                 [0.0, 0.0, 0.0, 0.0],
                 [0.0, 0.0, 0.0, 0.0],
                 [0.0, 0.0, 0.0, 0.0],
                 [0.0, 0.0, 0.0, 0.0],
                 [0.0, 0.0, 0.0, 0.0]]

C = np.transpose(A_temp)

B = np.array(B_temp)

sum_rows = [0.0] * 8
mean_rows = [0.0] * 8

for i in range(0, 8):
    for j in range(0, 4):
        sum_rows[i] += B[i][j]

for i in range(0, 8):
    mean_rows[i] = sum_rows[i] / 4.0

for i in range(0, 8):
    for j in range(0, 4):
        A_temp[i][j] = B[i][j] - mean_rows[i]

A = np.array(A_temp)

A_transpose = np.transpose(A)

C = 0.25 * (A.dot(A_transpose))

# U, s, Vh = np.linalg.svd(A, full_matrices=True)
mew_temp = np.array([mean_rows])

s, U = np.linalg.eig(C)
np.set_printoptions(formatter={'float_kind': '{:f}'.format})
np.set_printoptions(suppress=True)
# print ("U&&&&")
# print U.real
mew = np.transpose(mew_temp)

U_transpose = np.transpose(U)
scoringMatrix = U_transpose.dot(A)
# print ("scoring Matrix:")
# print scoringMatrix
scoringMatrix = np.delete(scoringMatrix, [1, 4, 5, 6, 7], 0)

print ("The Matrix B is:", B)
print ("\n")
print ("The Matrix A is:", A)
print ("\n")
print ("The co-variance matrix is: ", C.real)
print ("\n")
print ("The Eigen-Vector Matrix is:", U.real)
print ("\n")
print ("The Eigen Values are:", s.real)
print ("\n")
print ("The scoring Matrix is: ", scoringMatrix)
print ("\n")


def minEuclideanDist():
    low = float('inf')
    for i in range(0, 4):
        dist = distance.euclidean(W, scoringMatrix[:, i: i + 1])
        if dist < low:
            low = dist
        # print dist
    return low


u1 = U[0:8, 0:1]
u2 = U[0:8, 2:3]
u3 = U[0:8, 3:4]

Y1_temp = [[1],
           [5],
           [1],
           [5],
           [5],
           [1],
           [1],
           [3]]

Y_bar = Y1_temp - mew
Y1 = np.transpose(Y_bar)
W = np.array([Y1.dot(u1), Y1.dot(u2), Y1.dot(u3)])
result = minEuclideanDist()
print "Score for Y1 is:"
# print ("\n")
print result
print ("\n")
Y2_temp = [[-2],
           [3],
           [2],
           [3],
           [0],
           [2],
           [-1],
           [1]]

Y2_bar = Y2_temp - mew
Y2 = np.transpose(Y2_bar)
W = np.array([Y2.dot(u1), Y2.dot(u2), Y2.dot(u3)])
result = minEuclideanDist()
print ("Score for Y2 is:")
# print ("\n")
print result
print ("\n")
Y3_temp = [[2],
           [-3],
           [2],
           [3],
           [0],
           [0],
           [2],
           [-1]]

Y3_bar = Y3_temp - mew
Y3 = np.transpose(Y3_bar)

W = np.array([Y3.dot(u1), Y3.dot(u2), Y3.dot(u3)])
result = minEuclideanDist()
print ("Score for Y3 is:")
# print ("\n")
print result
print ("\n")

Y4_temp = [[2],
           [-2],
           [2],
           [2],
           [-1],
           [1],
           [2],
           [2]]

Y4_bar = Y4_temp - mew
Y4 = np.transpose(Y4_bar)

W = np.array([Y4.dot(u1), Y4.dot(u2), Y4.dot(u3)])
result = minEuclideanDist()
print ("Score for Y4 is:")
# print ("\n")
print result
print ("\n")
