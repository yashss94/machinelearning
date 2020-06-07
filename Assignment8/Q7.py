
data = [[8, 2],
        [5, 5],
        [9, 1],
        [4, 6],
        [7, 3]]

theta1 = 0.6
theta2 = 0.5

tow1 = 0.7
tow2 = 0.3

iters = 0

while iters < 2:
    temp_tow1 = 0
    temp_tow2 = 0

    mew1 = 0
    mew2 = 0
    if iters == 1:
        print("\nP-values: ")
    for i in range(len(data)):
        temp_p_a = tow1 * (theta1 ** data[i][0]) * ((1 - theta1) ** data[i][1])
        temp_p_b = tow2 * (theta2 ** data[i][0]) * ((1 - theta2) ** data[i][1])

        p_a = (temp_p_a / (temp_p_a + temp_p_b))
        p_b = (temp_p_b / (temp_p_a + temp_p_b))

        if iters == 1:
            print "p1",i+1,"    ",p_a
            print "p2",i+1,"    ",p_b

        temp_tow1 += p_a
        temp_tow2 += p_b

        mew1 += p_a * data[i][0]
        mew2 += p_b * data[i][0]

    tow1 = temp_tow1 / len(data)
    tow2 = temp_tow2 / len(data)

    theta1 = mew1 / (10 * len(data) * tow1)
    theta2 = mew2 / (10 * len(data) * tow2)

    iters += 1


