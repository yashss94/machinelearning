import numpy as np
import math


def EM(tow1, tow2, mew1, mew2, S1, S2, X):
    iters = 0

    det_S1 = np.linalg.det(S1)
    det_S2 = np.linalg.det(S2)

    inv_S1 = np.linalg.inv(S1)
    inv_S2 = np.linalg.inv(S2)

    while iters < 1:
        temp_tow1 = 0.0
        temp_tow2 = 0.0

        temp_mew1 = 0.0
        temp_mew2 = 0.0

        temp_S1 = [[0.0, 0.0],
                   [0.0, 0.0]]
        temp_S2 = [[0.0, 0.0],
                   [0.0, 0.0]]

        for i in range(len(X)):
            temp_p_a = tow1 * ((1 / (2 * math.pi * math.sqrt(det_S1))) * (math.e ** (
                        -0.5 * np.dot(np.dot(np.subtract(X[i], mew1), inv_S1), np.transpose(np.subtract(X[i], mew1))))))
            temp_p_b = tow2 * ((1 / (2 * math.pi * math.sqrt(det_S2))) * (math.e ** (
                        -0.5 * np.dot(np.dot(np.subtract(X[i], mew2), inv_S2), np.transpose(np.subtract(X[i], mew2))))))

            p_a = (temp_p_a / (temp_p_a + temp_p_b))
            p_b = (temp_p_b / (temp_p_a + temp_p_b))

            temp_S1 = np.add(temp_S1, np.dot(np.dot(p_a, np.atleast_2d(np.subtract(X[i], mew1)).T),np.atleast_2d(np.subtract(X[i], mew1))))
            temp_S2 = np.add(temp_S2, np.dot(np.dot(p_b, np.atleast_2d(np.subtract(X[i], mew2)).T),np.atleast_2d(np.subtract(X[i], mew2))))

            temp_tow1 += p_a
            temp_tow2 += p_b

            temp_mew1 += np.dot(p_a, X[i])
            temp_mew2 += np.dot(p_b, X[i])

        tow1 = temp_tow1 / len(X)
        tow2 = temp_tow2 / len(X)

        mew1 = np.divide(temp_mew1, temp_tow1)
        mew2 = np.divide(temp_mew2, temp_tow2)

        S1 = np.divide(temp_S1, temp_tow1)
        S2 = np.divide(temp_S2, temp_tow2)

        iters += 1

    return tow1, tow2, mew1, mew2, S1, S2


X = [
    [3.6, 79],
    [1.8, 54],
    [2.283, 62],
    [3.333, 74],
    [2.883, 55],
    [4.533, 85],
    [1.950, 51],
    [1.833, 54],
    [4.7, 88],
    [3.6, 85],
    [1.600, 52],
    [4.350, 85],
    [3.917, 84],
    [4.2, 78],
    [1.750, 62],
    [1.8, 51],
    [4.7, 83],
    [2.167, 52],
    [4.800, 84],
    [1.750, 47]]

mew1 = [2.5, 65.0]

S1 = [[1.0, 5.0],
      [5.0, 100.0]]

mew2 = [3.5, 70.0]

S2 = [[2.0, 10.0],
      [10.0, 200.0]]

tow1 = 0.6
tow2 = 0.4

tow1, tow2, mew1, mew2, S1, S2 = EM(tow1, tow2, mew1, mew2, S1, S2, X)


print("tow1", tow1)
print("tow2",tow2)
print("mew1",np.atleast_2d(mew1).T)
print("mew2",np.atleast_2d(mew2).T)

print("S1")
print S1
print("S2")
print S2
