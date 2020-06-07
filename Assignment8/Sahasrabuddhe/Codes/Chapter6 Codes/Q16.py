# Name: Yash Sahasrabuddhe
# Student ID: 014498887



import matplotlib.pyplot as plt
import math

def getneighbor(point, clusterIndex):
    global dictCore
    global neighbor
    global Clusters
    global X
    global visited

    point_neighbor = neighbor[point]
    for i in range(len(point_neighbor)):
        if visited[point_neighbor[i]]: continue
        visited[point_neighbor[i]] = True
        Clusters[clusterIndex].append(point_neighbor[i])
        if point_neighbor[i] in dictCore:
            getneighbor(point_neighbor[i], clusterIndex)


X = [[1.0, 5.0],
     [1.25, 5.35],
     [1.25, 5.75],
     [1.5, 6.25],
     [1.75, 6.75],
     [2.0, 6.5],
     [3.0, 7.75],
     [3.5, 8.25],
     [3.75, 8.75],
     [3.95, 9.1],
     [4.0, 8.5],
     [2.5, 7.25],
     [2.25, 7.75],
     [2.0, 6.5],
     [2.75, 8.25],
     [4.5, 8.9],
     [9.0, 5.0],
     [8.75, 5.85],
     [9.0, 6.25],
     [8.0, 7.0],
     [8.5, 6.25],
     [8.5, 6.75],
     [8.25, 7.65],
     [7.0, 8.25],
     [6.0, 8.75],
     [5.5, 8.25],
     [5.25, 8.75],
     [4.9, 8.75],
     [5.0, 8.5],
     [7.5, 7.75],
     [7.75, 8.25],
     [6.75, 8.0],
     [6.25, 8.25],
     [4.5, 8.9],
     [5.0, 1.0],
     [1.25, 4.65],
     [1.25, 4.25],
     [1.5, 3.75],
     [1.75, 3.25],
     [2.0, 3.5],
     [3.0, 2.25],
     [3.5, 1.75],
     [3.75, 8.75],
     [3.95, 0.9],
     [4.0, 1.5],
     [2.5, 2.75],
     [2.25, 2.25],
     [2.0, 3.5],
     [2.75, 1.75],
     [4.5, 1.1],
     [5.0, 9.0],
     [8.75, 5.15],
     [8.0, 2.25],
     [8.25, 3.0],
     [8.5, 4.75],
     [8.5, 4.25],
     [8.25, 3.35],
     [7.0, 1.75],
     [8.0, 3.5],
     [6.0, 1.25],
     [5.5, 1.75],
     [5.25, 1.25],
     [4.9, 1.25],
     [5.0, 1.5],
     [7.5, 2.25],
     [7.75, 2.75],
     [6.75, 2.0],
     [6.25, 1.75],
     [4.5, 1.1],
     [3.0, 4.5],
     [7.0, 4.5],
     [5.0, 3.0],
     [4.0, 3.35],
     [6.0, 3.35],
     [4.25, 3.25],
     [5.75, 3.25],
     [3.5, 3.75],
     [6.5, 3.75],
     [3.25, 4.0],
     [6.75, 4.0],
     [3.75, 3.55],
     [6.25, 3.55],
     [4.75, 3.05],
     [5.25, 3.05],
     [4.5, 3.15],
     [5.5, 3.15],
     [4.0, 6.5],
     [4.0, 6.75],
     [4.0, 6.25],
     [3.75, 6.5],
     [4.25, 6.5],
     [4.25, 6.75],
     [3.75, 6.25],
     [6.0, 6.5],
     [6.0, 6.75],
     [6.0, 6.25],
     [5.75, 6.75],
     [5.75, 6.25],
     [6.25, 6.75],
     [6.25, 6.25],
     [9.5, 9.5],
     [2.5, 9.5],
     [1.0, 8.0]]
dictCore = {}
neighbor = []
Clusters = []
visited = []
epsilon = 2
m = 10

for i in range(len(X)): neighbor.append([])

for i in range(0, len(X)):
    nCount = 0
    for j in range(i + 1, len(X)):
        dist = ((X[i][0] - X[j][0]) ** 2) + ((X[i][1] - X[j][1]) ** 2)
        dist = math.sqrt(dist)
        if dist <= epsilon:
            nCount += 1
            neighbor[i].append(j)
            neighbor[j].append(i)

    if nCount >= m:
        dictCore[i] = True

visited = []
for i in range(len(X)):
    visited.append(False)

clusterIndex = 0
for key in dictCore:
    if visited[key]: continue

    Clusters.append([])
    visited[key] = True
    getneighbor(key, clusterIndex)

    clusterIndex += 1

colors = ['red', 'green', 'black', 'cyan', 'pink', 'orange', 'violet', 'yellow', 'blue', 'grey', 'darkblue', 'purple']

clusteredpoints = 0
for i in range(len(Clusters)):
    x_coordinate = []
    y_coordinate = []
    currentCluster = Clusters[i]

    print("Elements in Cluster " + str(i+1) + ": " + str(len(currentCluster)))
    clusteredpoints += len(currentCluster)

    for j in range(len(currentCluster)):
        x_coordinate.append(X[currentCluster[j]][0])
        y_coordinate.append(X[currentCluster[j]][1])

        plt.scatter(x_coordinate, y_coordinate, color=colors[i])

print("Number of Outliers: " + str(len(X) - clusteredpoints))

print("m:",m)
print("epsilon:",epsilon)

plt.show()


