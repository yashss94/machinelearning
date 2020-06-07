import numpy as np
from scipy.spatial import distance

delta = np.array([[-1.1069, 1.2794, -2.6800, 2.5076],
                  [1.5480, 0.5484, -1.2085, -0.8879]])

u1 = np.array([[0.1641],
               [0.6278],
               [-0.2604],
               [-0.5389],
               [0.4637],
               [0.0752]])

u2 = np.array([[0.2443],
               [0.1070],
               [-0.8017],
               [0.4277],
               [-0.1373],
               [-0.2904]])

mew = np.array([[1.75],
                [1.75],
                [1.25],
                [2],
                [2],
                [1]])


def minEuclideanDist():
    low = float('inf')
    for i in range(0, 4):
        dist = distance.euclidean(W, delta[:, i: i + 1])
        if dist < low:
            low = dist
        # print dist
    return low


M1_temp = np.array([[1],
                    [-1],
                    [1],
                    [-1],
                    [-1],
                    [1]])

M1_bar = M1_temp - mew
M1 = np.transpose(M1_bar)
W = np.array([M1.dot(u1), M1.dot(u2)])
result = minEuclideanDist()
print ("Score of M1 is: ")
print result

M2_temp = np.array([[-2],
                    [2],
                    [2],
                    [-1],
                    [-2],
                    [2]])

M2_bar = M2_temp - mew
M2 = np.transpose(M2_bar)
W = np.array([M2.dot(u1), M2.dot(u2)])
result = minEuclideanDist()
print ("Score of M2 is: ")
print result


M3_temp = np.array([[1],
                    [3],
                    [0],
                    [1],
                    [3],
                    [1]])

M3_bar = M3_temp - mew
M3 = np.transpose(M3_bar)
W = np.array([M3.dot(u1), M3.dot(u2)])
result = minEuclideanDist()
print ("Score of M3 is: ")
print result

M4_temp = np.array([[2],
                    [3],
                    [1],
                    [1],
                    [-2],
                    [0]])

M4_bar = M4_temp - mew
M4 = np.transpose(M4_bar)
W = np.array([M4.dot(u1), M4.dot(u2)])
result = minEuclideanDist()
print ("Score of M4 is: ")
print result
