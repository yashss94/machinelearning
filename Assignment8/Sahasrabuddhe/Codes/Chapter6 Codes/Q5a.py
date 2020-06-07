# Name: Yash Sahasrabuddhe
# Student ID: 014498887

import pandas as pd
import matplotlib.pyplot as plt

iters = 100

def assignCluster(X):
    minDist = float("inf")
    cluster = 1
    for i in range(K):
        dist = 0
        for j in range(2):
            dist += (X[j] - centroids[i][j]) ** 2
        if dist < minDist:
            minDist = dist
            cluster = (i+1)
    
    return cluster
   
dataset = pd.DataFrame({
    'time': [3.600, 1.800, 2.283, 3.333, 2.883, 4.533, 1.950, 1.833, 4.700, 3.600, 1.600, 4.350, 3.917, 4.200, 1.750,1.800,
          4.700, 2.167, 4.800, 1.750],
    'wait_time': [79, 54, 62, 74, 55, 85, 51, 54, 88, 85, 52, 85, 84, 78, 62, 51, 83, 52, 84, 47]
})

K = 2 

min_time = dataset['time'].min()
max_time = dataset['time'].max()

min_wait_time = dataset['wait_time'].min()
max_wait_time = dataset['wait_time'].max()

time_bar = (max_time - min_time) / (K+1)
wait_bar = (max_wait_time - min_wait_time) / (K+1)

centroids = []
for i in range(K):
    centr_time = min_time + (i+1)*time_bar
    centr_wait = min_wait_time + (i+1)*wait_bar
    centroids.append([centr_time,centr_wait])
    
clusters = []
for index, rows in dataset.iterrows():
    
    result =[rows.time, rows.wait_time] 
       
    clusters.append(assignCluster(result))
    
dataset['clusters'] = clusters

for iter in range(iters):
    cluster1time = dataset[dataset['clusters']==1]['time'].mean()
    cluster2time = dataset[dataset['clusters']==2]['time'].mean()
    cluster1wait = dataset[dataset['clusters']==1]['wait_time'].mean()
    cluster2wait = dataset[dataset['clusters']==2]['wait_time'].mean()
    
    centroids[0][0] = cluster1time
    centroids[0][1] = cluster1wait
    centroids[1][0] = cluster2time
    centroids[1][1] = cluster2wait

        
    clusters.clear()
    
    for index, rows in dataset.iterrows():
        result = [rows.time, rows.wait_time]    
        clusters.append(assignCluster(result))
    dataset['clusters'] = clusters
    
colmap = {1:'r' , 2:'b'}
fig = plt.figure(figsize=(7,7))
colors = map(lambda x: colmap[x], clusters)
colors1 = list(colors)
plt.scatter(dataset['time'],dataset['wait_time'], color=colors1, alpha=0.4,edgecolor='3')
for idx, centroid in enumerate(centroids):
    plt.scatter(*centroid, color=colmap[idx+1])
plt.xlim(0,5)
plt.ylim(0,100)
plt.show()
    
