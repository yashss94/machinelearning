import numpy as np
import math

def EM(tau1, tau2, mew1, mew2, S1, S2, X):
    iters = 0
        
    while iters < 100:
        p_a = []
        p_b = []
        
        temp_tau1 = 0.0
        temp_tau2 = 0.0
    
        temp_mew1 = 0.0
        temp_mew2 = 0.0
    
        temp_S1 = [[0.0, 0.0],
                   [0.0, 0.0]]
        temp_S2 = [[0.0, 0.0],
                   [0.0, 0.0]]
        
        det_S1 = np.linalg.det(S1)
        det_S2 = np.linalg.det(S2)
    
        inv_S1 = np.linalg.inv(S1)
        inv_S2 = np.linalg.inv(S2)
        
        for i in range(len(X)):
            temp_p_a = tau1 * ((1 / (2 * math.pi * math.sqrt(det_S1))) * (math.e ** (
                            -0.5 * np.dot(np.dot(np.subtract(X[i], mew1), inv_S1), np.transpose(np.subtract(X[i], mew1))))))
            temp_p_b = tau2 * ((1 / (2 * math.pi * math.sqrt(det_S2))) * (math.e ** (
                            -0.5 * np.dot(np.dot(np.subtract(X[i], mew2), inv_S2), np.transpose(np.subtract(X[i], mew2))))))
    
            p_a.append(temp_p_a / (temp_p_a + temp_p_b))
            p_b.append(temp_p_b / (temp_p_a + temp_p_b))
            
            temp_tau1 += p_a[i]
            temp_tau2 += p_b[i]
            
            temp_mew1 += np.dot(p_a[i], X[i])
            temp_mew2 += np.dot(p_b[i], X[i])
    
        tau1 = temp_tau1 / len(X)
        tau2 = temp_tau2 / len(X)
    
        mew1 = np.divide(temp_mew1, temp_tau1)
        mew2 = np.divide(temp_mew2, temp_tau2)
        
        for i in range(len(X)):
            temp_S1 = np.add(temp_S1, np.dot(np.dot(p_a[i], np.atleast_2d(np.subtract(X[i], mew1)).T),np.atleast_2d(np.subtract(X[i], mew1))))
            temp_S2 = np.add(temp_S2, np.dot(np.dot(p_b[i], np.atleast_2d(np.subtract(X[i], mew2)).T),np.atleast_2d(np.subtract(X[i], mew2))))
    
        S1 = np.divide(temp_S1, temp_tau1)
        S2 = np.divide(temp_S2, temp_tau2)
    
        iters += 1
    
    return tau1, tau2, mew1, mew2, S1, S2, p_a, p_b



def plotClusters(p_a, p_b, X, n, centroids):
    clusters=[]
    for i in range(n):
        if(p_a[i] > p_b[i]):
            clusters.append(1)
        else:
            clusters.append(2)
    colormap = {1:'blue' , 2:'orange'}
    colors = map(lambda x: colormap[x], clusters)
    colors1 = list(colors)
   
    plt.scatter(X[:,0], X[:,1], color=colors1, alpha=0.5,edgecolor='5')
    for idx, centroid in enumerate(centroids):
        plt.scatter(*centroid, color=colormap[idx+1])
    
    plt.xlim(0,5)
    plt.ylim(40,90)
    plt.show()

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

tau1 = 0.6
tau2 = 0.4

p_a = []
p_b = []

tau1, tau2, mew1, mew2, S1, S2, p_a, p_b = EM(tau1, tau2, mew1, mew2, S1, S2, X)

print("Converged values of Tau, Theta and S for (d) part using given initializations:")
print("Tau1",tau1)
print("Tau2",tau2)
print("Mew1",np.atleast_2d(mew1).T)
print("Mew2",np.atleast_2d(mew2).T)

print("S1")
print(S1)
print("S2")
print(S2)

centroids = np.array([mew1, mew2])

print("The centroids of the data will be the converged values of mew")
print("The centroids are: ",centroids)
plotClusters(p_a, p_b, np.array(X), len(X), centroids)