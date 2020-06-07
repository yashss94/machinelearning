import numpy as np
import math

def EM(tau1, tau2, tau3, mew1, mew2, mew3, S1, S2, S3, X):
    iters = 0
        
    while iters < 100:
        p_a = []
        p_b = []
        p_c = []
        
        temp_tau1 = 0.0
        temp_tau2 = 0.0
        temp_tau3 = 0.0

        temp_mew1 = 0.0
        temp_mew2 = 0.0
        temp_mew3 = 0.0
        
        temp_S1 = [[0.0, 0.0],
                   [0.0, 0.0]]
        temp_S2 = [[0.0, 0.0],
                   [0.0, 0.0]]
        temp_S3 = [[0.0, 0.0],
                   [0.0, 0.0]]
        
        det_S1 = np.linalg.det(S1)
        det_S2 = np.linalg.det(S2)
        det_S3 = np.linalg.det(S3)

        inv_S1 = np.linalg.inv(S1)
        inv_S2 = np.linalg.inv(S2)
        inv_S3 = np.linalg.inv(S3)

        for i in range(len(X)):
            temp_p_a = tau1 * ((1 / (2 * math.pi * math.sqrt(det_S1))) * (math.e ** (
                            -0.5 * np.dot(np.dot(np.subtract(X[i], mew1), inv_S1), np.transpose(np.subtract(X[i], mew1))))))
            temp_p_b = tau2 * ((1 / (2 * math.pi * math.sqrt(det_S2))) * (math.e ** (
                            -0.5 * np.dot(np.dot(np.subtract(X[i], mew2), inv_S2), np.transpose(np.subtract(X[i], mew2))))))
            temp_p_c = tau3 * ((1 / (2 * math.pi * math.sqrt(det_S3))) * (math.e ** (
                            -0.5 * np.dot(np.dot(np.subtract(X[i], mew3), inv_S3), np.transpose(np.subtract(X[i], mew3))))))
    
            p_a.append(temp_p_a / (temp_p_a + temp_p_b + temp_p_c))
            p_b.append(temp_p_b / (temp_p_a + temp_p_b + temp_p_c))
            p_c.append(temp_p_c / (temp_p_a + temp_p_b + temp_p_c))

            
            temp_tau1 += p_a[i]
            temp_tau2 += p_b[i]
            temp_tau3 += p_c[i]

            temp_mew1 += np.dot(p_a[i], X[i])
            temp_mew2 += np.dot(p_b[i], X[i])
            temp_mew3 += np.dot(p_c[i], X[i])
                
        tau1 = temp_tau1 / len(X)
        tau2 = temp_tau2 / len(X)
        tau3 = temp_tau3 / len(X)
    
        mew1 = np.divide(temp_mew1, temp_tau1)
        mew2 = np.divide(temp_mew2, temp_tau2)
        mew3 = np.divide(temp_mew3, temp_tau3)
        
        for i in range(len(X)):
            temp_S1 = np.add(temp_S1, np.dot(np.dot(p_a[i], np.atleast_2d(np.subtract(X[i], mew1)).T),np.atleast_2d(np.subtract(X[i], mew1))))
            temp_S2 = np.add(temp_S2, np.dot(np.dot(p_b[i], np.atleast_2d(np.subtract(X[i], mew2)).T),np.atleast_2d(np.subtract(X[i], mew2))))
            temp_S3 = np.add(temp_S3, np.dot(np.dot(p_c[i], np.atleast_2d(np.subtract(X[i], mew3)).T),np.atleast_2d(np.subtract(X[i], mew3))))

        S1 = np.divide(temp_S1, temp_tau1)
        S2 = np.divide(temp_S2, temp_tau2)
        S3 = np.divide(temp_S3, temp_tau3)
    
        iters += 1
    
    return tau1, tau2, tau3, mew1, mew2, mew3, S1, S2, S3



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

mew3 = [3.0, 67.5]
S3 = [[1.5, 7.5],
        [7.5, 150.0]]
    

tau1 = 0.3
tau2 = 0.5
tau3 = 0.2


p_a = []
p_b = []
p_c = []

tau1, tau2, tau3, mew1, mew2, mew3, S1, S2, S3 = EM(tau1, tau2, tau3, mew1, mew2, mew3, S1, S2, S3, X)

print("Converged values of Tau, Theta and S for (c) part for three clusters:")
print("Tau1",tau1)
print("Tau2",tau2)
print("Tau3",tau3)

print("Mew1",np.atleast_2d(mew1).T)
print("Mew2",np.atleast_2d(mew2).T)
print("Mew3",np.atleast_2d(mew3).T)


print("S1")
print(S1)
print("S2")
print(S2)
print("S3")
print(S3)


