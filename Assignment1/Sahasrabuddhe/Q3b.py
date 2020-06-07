def alphaPass(A, B, Pi, ObsSeq):
    N=2
    T=4
    alphaProb = []
    for m in range(0,2):
        alphaProb.append([1,1,1,1])
    for i in range(0, 2):
        alphaProb[i][0] = Pi[i]*B[i][ObsSeq[0]]
    for j in range(1, 4):
        for k in range(0, 2):
            alphaProb[k][j] = 0
            for l in range(0, 2):
                alphaProb[k][j] += (alphaProb[l][j-1]*A[l][k])
            alphaProb[k][j] *= B[k][ObsSeq[j]]
    aprob = alphaProb[N-2][T-1] + alphaProb[N-1][T-1]
    return aprob


Pi = [0.6, 0.4]
A = [[0.7, 0.3], [0.4, 0.6]]
B = [[0.1, 0.4, 0.5], [0.7, 0.2, 0.1]]
FinalProb = 0
alphaprobs = 0
for O0 in range(0, 3):
    for O1 in range(0, 3):
        for O2 in range(0, 3):
            for O3 in range(0, 3):
                ObsSeq = [O0, O1, O2, O3]
                alphaprobs = alphaPass(A, B, Pi, ObsSeq)
                FinalProb += alphaprobs
                print ("Observation Sequence     Alpha Probability")
                print (str(ObsSeq) + "\t" + str(alphaprobs))
print ("Final Sum of Probabilities is: ", str(FinalProb))












