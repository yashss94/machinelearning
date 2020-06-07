Pi = [0.6, 0.4]
A = [[0.7, 0.3], [0.4, 0.6]]
B = [[0.1, 0.4, 0.5], [0.7, 0.2, 0.1]]
FinalProb = 0
for O0 in range(0, 3):
    for O1 in range(0, 3):
        for O2 in range(0, 3):
            for O3 in range(0, 3):
                SequenceProb = 0
                for t in range(0, 2):
                    for x in range(0, 2):
                        for y in range(0, 2):
                            for z in range(0, 2):
                                SequenceProb += (Pi[t]*B[t][O0]*A[t][x]*B[x][O1]*A[x][y]*B[y][O2]*A[y][z]*B[z][O3])
                                round(SequenceProb, 5)
                FinalProb += SequenceProb
                ObsSeq = [O0, O1, O2, O3]
                print ("Observation Sequence\tProbability")
                print (str(ObsSeq) + "\t\t\t" + str(SequenceProb))
print ("Final Sum of Probabilities is: ", str(FinalProb))












