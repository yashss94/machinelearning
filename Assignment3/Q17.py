import numpy as np
import math

u1 = np.array([[0.1641],
               [0.6278],
               [-0.2604],
               [-0.5389],
               [0.4637],
               [0.0752]])

u2 = np.array([[0.2443],
               [0.1070],
               [-0.8017],
               [0.4277],
               [-0.1373],
               [-0.2904]])

u3 = np.array([[-0.0710],
               [0.2934],
               [0.3952],
               [0.3439],
               [0.3644],
               [-0.7083]])

lambda1 = 4.0414
lambda2 = 2.2239
lambda3 = 1.7237

root_lambda1 = math.sqrt(lambda1)
root_lambda2 = math.sqrt(lambda2)
root_lambda3 = math.sqrt(lambda3)

max_positive_corr = max(u1)
max_negative_corr = min(u1)
print ("The feature with greatest effect with positive correlation on u1 is: ")
print max_positive_corr
print ("The feature with greatest effect with negative correlation on u1 is: ")
print max_negative_corr
print ("The feature with greatest overall effect on u1 is: ")
overall_max_effect = max((abs(max_positive_corr)), (abs(max_negative_corr)))
print overall_max_effect

CLV1 = root_lambda1 * u1
print("The component loading vector of component u1 is:")
print CLV1
print("The component loading vector of component u2 is:")
CLV2 = root_lambda2 * u2
print CLV2
CLV3 = root_lambda3 * u3
print("The component loading vector of component u3 is:")
print CLV3

significant_CLV = CLV1 + CLV2


def sorting(num_array):
    return sorted(num_array, key=abs)


def find_rank(to_check):
    temp_list = []
    temp_dict = {}
    order = 1
    rank = 1
    weight = 0
    for i in to_check:
        temp_list.append(float(i[0]))
        temp_dict[float([i][0])] = [order, rank, weight]
        order = order + 1
    temp_list = sorting(temp_list)
    rank = 6
    for i in temp_list:
        temp_dict[i][1] = rank
        rank = rank - 1
    for i in temp_dict:
        temp_dict[i][2] = 500 / temp_dict[i][1]
        print("Rank of", i, "is ", temp_dict[i][1], "and the relative weight is", temp_dict[i][2])


print("The relative important features using the two most significant components are:")
find_rank(significant_CLV)

