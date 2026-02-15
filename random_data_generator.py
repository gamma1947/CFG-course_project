from pathlib import Path
import numpy as np
import numpy.random as rand

output_file = Path("/home/photon/CFG-course_project/random_data.txt")
output_file.parent.mkdir(parents=True, exist_ok=True)

# transition probabilities

# indexes = {"A":0, "T":1, "G":2, "C":3}
indexes_1 = {0:'A', 1:'T', 2:'G', 3:'C'}


transition_matrix = np.zeros((4,4))
for j in range(len(transition_matrix)):
    whole = 100
    for i in range(len(transition_matrix)-1):
        # print(i)
        cut = rand.uniform(0, whole/2)
        # print(cut)
        left_over = whole - cut
        # print(whole)
        transition_matrix[j,i] = (cut/100)
        whole = left_over
        # print(whole)
    transition_matrix[j,-1] = 1-sum(transition_matrix[j])

print(transition_matrix)


c_idx = 0
length = 100000
string = []

for _ in range(length):
    c_char = indexes_1[c_idx]
    string.append(c_char)
    next_idx = rand.choice([0, 1, 2, 3], p = transition_matrix[c_idx,])
    c_idx = next_idx
    # next_char = indexes_1[c_idx]

final_str = "".join(string)
# print(final_str)
with open(output_file, 'w') as o_f:
    o_f.write(final_str)

