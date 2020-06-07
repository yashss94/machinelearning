
def EM(tau1, tau2, theta1, theta2, data):
    iters = 0

    while iters < 150:
        temp_tau1 = 0
        temp_tau2 = 0

        mew1 = 0
        mew2 = 0

        for i in range(len(data)):

            temp_p_a = tau1 * (theta1 ** data[i][0]) * ((1 - theta1) ** data[i][1])
            temp_p_b = tau2 * (theta2 ** data[i][0]) * ((1 - theta2) ** data[i][1])

            p_a = (temp_p_a / (temp_p_a + temp_p_b))
            p_b = (temp_p_b / (temp_p_a + temp_p_b))

            temp_tau1 += p_a
            temp_tau2 += p_b

            mew1 += p_a * data[i][0]
            mew2 += p_b * data[i][0]

        tau1 = temp_tau1 / len(data)
        tau2 = temp_tau2 / len(data)

        theta1 = mew1 / (10 * len(data) * tau1)
        theta2 = mew2 / (10 * len(data) * tau2)

        iters += 1

    return tau1, tau2, theta1, theta2


data = [[8, 2],
        [5, 5],
        [9, 1],
        [4, 6],
        [7, 3]]

theta1 = 0.6
theta2 = 0.5

tau1 = 0.7
tau2 = 0.3

tau1, tau2, theta1, theta2 = EM(tau1, tau2, theta1, theta2, data)
print "\n"
print "Value of theta1 is: ", theta1
print "Value of theta2 is: ", theta2
print "\n"
print "Value of tau1 is: ", tau1
print"Value of tau2 is: ", tau2
