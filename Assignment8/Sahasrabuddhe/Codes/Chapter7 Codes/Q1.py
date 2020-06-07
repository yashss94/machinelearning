# Name: Yash Sahasrabuddhe
# Student ID: 014498887


from sklearn.neighbors import KNeighborsClassifier
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np 

X = np.array([
    [0.5, 3.0],
    [1.0, 4.25],
    [1.5, 2.0],
    [2.0, 2.75],
    [2.5, 1.65],
    [3.0, 2.7],
    [3.5, 1.0],
    [4.0, 2.5],
    [4.5, 2.1],
    [5.0, 2.75],
    [0.5, 1.75],
    [1.5, 1.5],
    [2.5, 4.0],
    [2.5, 2.1],
    [3.0, 1.5],
    [3.5, 1.85],
    [4.0, 3.5],
    [5.0, 1.45]
])

Y = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0])

k = 3   
h = 0.02   
map_color = ListedColormap(['#FFC4B4', '#C4FFB4', '#89FFFF'])
point_color = ListedColormap(['#E60000', '#00E600', '#0000CC'])

clf = KNeighborsClassifier(n_neighbors=k)
clf.fit(X, Y)

X_min, X_max = X[:, 0].min() - 1, X[:, 0].max() + 1
Y_min, Y_max = X[:, 0].min() - 1, X[:, 0].max() + 1
xx, yy = np.meshgrid(np.arange(X_min, X_max, h), np.arange(Y_min, Y_max, h))  
Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])


Z = Z.reshape(xx.shape)
plt.figure()
plt.pcolormesh(xx, yy, Z, cmap=map_color)
plt.scatter(X[:, 0], X[:, 1], c=Y, cmap=point_color)
plt.xlim(xx.min(), xx.max())
plt.ylim(yy.min(), yy.max())


plt.show()

