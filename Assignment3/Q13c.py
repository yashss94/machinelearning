import numpy as np
from scipy.spatial import distance

M1 = [[1],
      [-1],
      [1],
      [-1],
      [-1],
      [1]]

M2 = [[-2],
      [2],
      [2],
      [-1],
      [-2],
      [2]]
M3 = [[1],
      [3],
      [0],
      [1],
      [3],
      [1]]
M4 = [[2],
      [3],
      [1],
      [1],
      [-2],
      [0]]

M = np.column_stack((M1, M2, M3, M4))

M_A = np.array([[0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0]])

sum_rows = [0.0] * 6
mew_malign = [0.0] * 6

for i in range(0, 6):
    for j in range(0, 4):
        sum_rows[i] += M[i][j]

for i in range(0, 6):
    mew_malign[i] = sum_rows[i] / 4.0

for i in range(0, 6):
    for j in range(0, 4):
        M_A[i][j] = M[i][j] - mew_malign[i]

mew_malign = np.transpose([mew_malign])
M_A_transpose = np.transpose(M_A)

C_M = (1.0 / 4) * (np.dot(M_A, M_A_transpose))

s_M, U_M = np.linalg.eig(C_M)
np.set_printoptions(formatter={'float_kind': '{:f}'.format})
np.set_printoptions(suppress=True)
print("Eigen Vectors for Malign Samples are:")
print U_M
print("Eigen Values for Malign Samples are:")
print s_M
U_M_transpose = np.transpose(U_M)

scoreMatM = U_M_transpose.dot(M_A)
# print("$$$$$")
# print scoreMatM
delta = np.delete(scoreMatM, [1, 4, 5], 0)
print("Delta Matrix for Malign Samples is:")
print delta
u1_mal = U_M[0:6, 0:1]
u2_mal = U_M[0:6, 2:3]
u3_mal = U_M[0:6, 3:4]

B1 = np.array([[-1],
               [2],
               [1],
               [2],
               [-1],
               [0]])

B2 = np.array([[-2],
               [1],
               [2],
               [3],
               [2],
               [1]])

B3 = np.array([[-1],
               [3],
               [0],
               [1],
               [3],
               [-1]])

B4 = np.array([[0],
               [2],
               [3],
               [1],
               [1],
               [-2]])

B = np.column_stack((B1, B2, B3, B4))

A = np.array([[0.0, 0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0, 0.0]])

sum_rows = [0.0] * 6
mew = [0.0] * 6

for i in range(0, 6):
    for j in range(0, 4):
        sum_rows[i] += B[i][j]

for i in range(0, 6):
    mew[i] = sum_rows[i] / 4.0

for i in range(0, 6):
    for j in range(0, 4):
        A[i][j] = B[i][j] - mew[i]

A_transpose = np.transpose(A)

C = (1.0 / 4) * (np.dot(A, A_transpose))
# print C

mew_benign_temp = np.array([mew])
mew_benign = np.transpose(mew_benign_temp)

s, U = np.linalg.eig(C)
np.set_printoptions(formatter={'float_kind': '{:f}'.format})
np.set_printoptions(suppress=True)
# print s
# print U
U_transpose = np.transpose(U)
print("Eigen Vectors for Benign Samples are:")
print U
print("Eigen Values for Benign Samples are:")
print s
scoringMatrix = U_transpose.dot(A)
InverseDelta = np.delete(scoringMatrix, [3, 4, 5], 0)
print("The Inverse Delta Matrix for Benign Samples is:")
print InverseDelta


def minEuclideanDistMalign():
    low = float('inf')
    for i in range(0, 4):
        dist = distance.euclidean(W_malign, delta[:, i: i + 1])
        dist = float(('{0:1.2f}'.format(dist)))
        if dist < low:
            low = dist
        # print dist
    return low


def minEuclideanDistBenign():
    low = float('inf')
    for i in range(0, 4):
        dist = distance.euclidean(W_benign, InverseDelta[:, i: i + 1])
        dist = float(('{0:1.2f}'.format(dist)))
        if dist < low:
            low = dist
        # print dist
    return low


u1_ben = U[0:6, 0:1]
u2_ben = U[0:6, 1:2]
u3_ben = U[0:6, 2:3]


Y1 = [[1],
      [5],
      [1],
      [5],
      [5],
      [1]]

Y1_bar_benign = Y1 - mew_benign
Y1_bar_malign = Y1 - mew_malign

Y1_benign_transpose = np.transpose(Y1_bar_benign)
Y1_malign_transpose = np.transpose(Y1_bar_malign)

W_benign = np.array([Y1_benign_transpose.dot(u1_ben), Y1_benign_transpose.dot(u2_ben), Y1_benign_transpose.dot(u3_ben)])
W_malign = np.array([Y1_malign_transpose.dot(u1_mal), Y1_malign_transpose.dot(u2_mal), Y1_malign_transpose.dot(u3_mal)])
# print("W_Benign")
# print W_benign
# print("W_Malign")
# print W_malign

# print("Benign scores")
result1 = minEuclideanDistBenign()
# print("Malign scores")
result2 = minEuclideanDistMalign()

if result1 < result2:
    print("Y1 is benign")
else:
    print("Y1 is malign")

Y2 = [[-2],
      [3],
      [2],
      [3],
      [0],
      [2]]

Y2_bar_benign = Y2 - mew_benign
Y2_bar_malign = Y2 - mew_malign

Y2_benign_transpose = np.transpose(Y2_bar_benign)
Y2_malign_transpose = np.transpose(Y2_bar_malign)

W_benign = np.array([Y2_benign_transpose.dot(u1_ben), Y2_benign_transpose.dot(u2_ben), Y2_benign_transpose.dot(u3_ben)])
W_malign = np.array([Y2_malign_transpose.dot(u1_mal), Y2_malign_transpose.dot(u2_mal), Y2_malign_transpose.dot(u3_mal)])
# print("W_Benign")
# print W_benign
# print("W_Malign")
# print W_malign

# print("Benign scores")
result1 = minEuclideanDistBenign()
# print("Malign scores")
result2 = minEuclideanDistMalign()

if result1 < result2:
    print("Y2 is benign")
else:
    print("Y2 is malign")

Y3 = [[2],
      [-3],
      [2],
      [3],
      [0],
      [0]]

Y3_bar_benign = Y3 - mew_benign
Y3_bar_malign = Y3 - mew_malign

Y3_benign_transpose = np.transpose(Y3_bar_benign)
Y3_malign_transpose = np.transpose(Y3_bar_malign)

W_benign = np.array([Y3_benign_transpose.dot(u1_ben), Y3_benign_transpose.dot(u2_ben), Y3_benign_transpose.dot(u3_ben)])
W_malign = np.array([Y3_malign_transpose.dot(u1_mal), Y3_malign_transpose.dot(u2_mal), Y3_malign_transpose.dot(u3_mal)])
# print("W_Benign")
# print W_benign
# print("W_Malign")
# print W_malign

# print("Benign scores")
result1 = minEuclideanDistBenign()
# print("Malign scores")
result2 = minEuclideanDistMalign()

if result1 < result2:
    print("Y3 is benign")
else:
    print("Y3 is malign")

Y4 = [[2],
      [-2],
      [2],
      [2],
      [-1],
      [1]]

Y4_bar_benign = Y4 - mew_benign
Y4_bar_malign = Y4 - mew_malign

Y4_benign_transpose = np.transpose(Y4_bar_benign)
Y4_malign_transpose = np.transpose(Y4_bar_malign)

W_benign = np.array([Y4_benign_transpose.dot(u1_ben), Y4_benign_transpose.dot(u2_ben), Y4_benign_transpose.dot(u3_ben)])
W_malign = np.array([Y4_malign_transpose.dot(u1_mal), Y4_malign_transpose.dot(u2_mal), Y4_malign_transpose.dot(u3_mal)])
# print("W_Benign")
# print W_benign
# print("W_Malign")
# print W_malign

# print("Benign scores")
result1 = minEuclideanDistBenign()
# print("Malign scores")
result2 = minEuclideanDistMalign()

if result1 < result2:
    print("Y4 is benign")
else:
    print("Y4 is malign")
