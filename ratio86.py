import numpy as np
from matplotlib import pyplot as plt

#Read data
data = np.loadtxt("toiee.dat",usecols=(0,1,2,3,4,5,6))

A = data[:,1]
Z = data[:,2]
J = data[:,3]
n = data[:,5]
p = data[:,4]
E = data[:,6]

ratio = []
N = []
a =[]

for i in range(len(A)):
    if J[i] == 8 and n[i]==1 and p[i]==1:                                                 #Look for yrast 8+ state
        for j in range(len(A)): 
            if j<i and J[j] == 6 and n[i]==1 and p[i]==1 and A[i]==A[j] and Z[i]==Z[j]:   #Look for yrast 6+ state in same nucleus
                ratio.append(E[i]/E[j])                                                   #Calculate ratio
                a.append(A[i])
                N.append(A[i]-Z[i])
                break
                

plt.xlabel('A')
plt.ylabel('E(8+)/E(6+)')            
plt.plot(a,ratio,'ro')
plt.savefig('ratio.png')
plt.show()
