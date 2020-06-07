import math as m

X0 = [0.6, 0.1, 0.8, 0.3, 0.7, 0.7, 0.2]
X1 = [0.4, 0.2, 0.6, 0.7, 0.3, 0.7, 0.9]
Z = [1, 0, 0, 1, 1, 0, 1]
W = [1, 2, -1, 1, -2, 1]

learn_rate = 0.1
epoch = 10000

count = 0

while count < epoch:
    for i in range(len(X0)):
        v0 = W[0]
        v1 = W[1]
        v2 = W[2]
        v3 = W[3]
        v4 = W[4]
        v5 = W[5]
        v6 = X0[i] * v0 + X1[i] * v2
        v7 = X0[i] * v1 + X1[i] * v3
        v8 = 1 + m.exp(-v6)
        v9 = 1 + m.exp(-v7)
        v10 = v4 / v8
        v11 = v5 / v9
        v12 = (v10 + v11 - Z[i]) ** 2 / 2
        z = v12
        dz = 1
        dv11 = v10 + v11 - Z[i]
        dv10 = v10 + v11 - Z[i]
        dv9 = -v5 / (v9 ** 2) * dv11
        dv8 = -v4 / (v8 ** 2) * dv10
        dv7 = -m.exp(-v7) * dv9
        dv6 = -m.exp(-v6) * dv8
        dv5 = dv11 / v9
        dv4 = dv10 / v8
        dv3 = X1[i] * dv7
        dv2 = X1[i] * dv6
        dv1 = X0[i] * dv7
        dv0 = X0[i] * dv6
        W[0] -= learn_rate * dv0
        W[1] -= learn_rate * dv1
        W[2] -= learn_rate * dv2
        W[3] -= learn_rate * dv3
        W[4] -= learn_rate * dv4
        W[5] -= learn_rate * dv5
    count += 1

print("The updated weights are:")
print(W)

prediction = []
X0_Test = [0.55, 0.32, 0.24, 0.86, 0.53, 0.46, 0.16, 0.52, 0.46, 0.96]
X1_Test = [0.11, 0.21, 0.64, 0.68, 0.79, 0.54, 0.51, 0.94, 0.87, 0.63]
Z_Test = [1, 0, 1, 0, 0, 1, 1, 0, 1, 0]

for i in range(len(X0_Test)):
    Y = (W[4] / (1 + m.exp(-(W[0] * X0_Test[i] + W[2] * X1_Test[i])))) + (W[5] / (1 + m.exp(-(W[1] * X0_Test[i] + W[3] * X1_Test[i]))))
    if Y > 0.5:
        prediction.append(1)
    else:
        prediction.append(0)
    print("Y for X0_test ", X0_Test[i], "and for X1_test ", X1_Test[i], " is: ", Y)

print("\nPredicted Values for Test Case:")
print(prediction)

print("\nActual Values for Test Case")
print(Z_Test)

count = 0
for i in range(len(prediction)):
    if prediction[i] == Z_Test[i]:
        count += 1

accuracy = count * 1.0 / len(prediction) * 100
print("Accuracy is:", accuracy)
