import numpy as np

A = np.array([[-0.5, 1, -0.5],
              [0.5, -1, 0.5]])

B = np.array([[-0.5, 0.5],
              [1, -1],
              [-0.5, 0.5]])

A_B = A.dot(B)
B_A = B.dot(A)

B_A = 0.5 * B_A

print B_A
