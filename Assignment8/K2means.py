# -*- coding: utf-8 -*-
"""
Created on Sat May  9 20:37:48 2020

@author: hardi
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

iterations = 100

def centroidCluster():
    cluster1duration = df[df['clusters']==1]['duration'].mean()
    cluster2duration = df[df['clusters']==2]['duration'].mean()
    cluster1wait = df[df['clusters']==1]['wait'].mean()
    cluster2wait = df[df['clusters']==2]['wait'].mean()
    
    centroids[0][0] = cluster1duration
    centroids[0][1] = cluster1wait
    centroids[1][0] = cluster2duration
    centroids[1][1] = cluster2wait
    
def clusterAssign(X):
    minDist = float("inf")
    cluster = 1
    for i in range(2):
        dist = 0
        for j in range(2):
            dist += (X[j] - centroids[i][j]) ** 2
        if dist < minDist:
            minDist = dist
            cluster = (i+1)
    
    return cluster


df = pd.DataFrame({
    'duration': [3.600, 1.800, 2.283, 3.333, 2.883, 4.533, 1.950, 1.833, 4.700, 3.600, 1.600, 4.350, 3.917, 4.200, 1.750,1.800,
          4.700, 2.167, 4.800, 1.750],
    'wait': [79, 54, 62, 74, 55, 85, 51, 54, 88, 85, 52, 85, 84, 78, 62, 51, 83, 52, 84, 47]
})

# Number of clusters = 2
K = 2 

minDuration = df['duration'].min()
maxDuration = df['duration'].max()

minWait = df['wait'].min()
maxWait = df['wait'].max()

dur_bar = (maxDuration - minDuration) / (K+1)
wait_bar = (maxWait - minWait) / (K+1)

centroids = []
for i in range(K):
    centr_dur = minDuration + (i+1)*dur_bar
    centr_wait = minWait + (i+1)*wait_bar
    centroids.append([centr_dur,centr_wait])
    
clusters = []
for index, rows in df.iterrows():
    
    # Create list for the current row 
    entry =[rows.duration, rows.wait] 
       
    clusters.append(clusterAssign(entry))
    
df['clusters'] = clusters

for iter in range(iterations+100):
    centroidCluster()
    
    clusters.clear()
    
    for index, rows in df.iterrows():
        # Create list for the current row 
        entry =[rows.duration, rows.wait]    
        clusters.append(clusterAssign(entry))
    df['clusters'] = clusters
    
colmap = {1:'r' , 2:'g'}
fig = plt.figure(figsize=(5,5))
colors = map(lambda x: colmap[x], clusters)
colors1 = list(colors)
plt.scatter(df['duration'],df['wait'], color=colors1, alpha=0.5,edgecolor='2')
for idx, centroid in enumerate(centroids):
    plt.scatter(*centroid, color=colmap[idx+1])
plt.xlim(0,5)
plt.ylim(0,90)
plt.show()
    
