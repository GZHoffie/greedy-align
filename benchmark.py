import Levenshtein

with open("/mnt/c/Users/zhenh/trace_recon/Output_greedy.txt") as f:
    reconstructed_strands = f.readlines()

reconstructed_strands = list(map(lambda x:x.strip(), reconstructed_strands))


#reconstructed_strands = # your reconstructed strands go here
answer_file = open("/mnt/c/Users/zhenh/trace_recon/our_nanopore_refs.txt","r")


answer = answer_file.readlines()
answer = list(map(lambda x:x.strip(), answer))

print(len(reconstructed_strands), len(answer))


offset = 0

total_num_bases = 0
total_edit_dist = 0
total_num_strands = 0
total_errors = 0
for i in range(len(answer)):
        total_num_strands += 1
        total_edit_dist += Levenshtein.distance(reconstructed_strands[i-offset], answer[i])
        total_num_bases += len(answer[i])
        for j in range(len(answer[i])):
            if j>= len(reconstructed_strands[i-offset]):
                total_errors+= len(answer[i]) - len(reconstructed_strands[i-offset])
                break
            elif reconstructed_strands[i-offset][j] != answer[i][j]:
                total_errors+= 1

print("(Avg Accuracy, Avg Edit Dist):")
print('\t',1- (total_errors/total_num_bases), total_edit_dist/total_num_strands)
print('\n--------------------------------------------------\n')




# # Below For Plotting Error Rate

# answer_file = open("EncodedStrands.txt",'r')
# answer = answer_file.readlines()
# answer = list(map(lambda x:x.strip(), answer))

print('Actual:\t',answer[0])
print('Recon: \t',reconstructed_strands[0])

correct = 0
for i in range(len(reconstructed_strands)):
    if answer[i] == reconstructed_strands[i]:
        correct+= 1
    #else:
    #    print(i)
    #    print('Actual:\t',answer[i])
    #    print('Recon: \t', reconstructed_strands[i])
        
print("number exact match:", str(correct)+'/'+str(len(answer)))

index_wrong_counts = [0]*len(answer[0])
for i in range(len(reconstructed_strands)):
    for j in range(len(answer[0])):
        if j >= len(reconstructed_strands[i]):
            index_wrong_counts[j] += 1
            continue
        # print(i,j)
        # print(len(reconstructed_strands[i]),len(answer[i]))
        # print(answer[i])
        if reconstructed_strands[i][j] != answer[i][j]:
            index_wrong_counts[j] += 1

import matplotlib.pyplot as plt
import numpy as np

xpoints = list(range(len(answer[0])))
ypoints = list(map(lambda x: x/len(reconstructed_strands), index_wrong_counts))
plt.plot(xpoints,ypoints)
plt.savefig('error_rates.png')