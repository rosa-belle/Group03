import numpy as np
from matplotlib import pyplot as plt
import sys

#Read data and put columns into lists
data = np.loadtxt("aud16.dat")

Z = data[:,0]
A = data[:,1]
N = data[:,5]
BE = data[:,2]

sys.stdout=open('pdecay_out','w')
print('Z      N      S_p      S_2p')

for i in range(len(Z)):                                    
    #print(Z[i])
    for j in range(len(N)):
        S_p = 0
        S_2p = 0
        if Z[i]-1==Z[j] and N[j]==N[i]:                       #Find isotope with Z-1
            S_p = BE[i]-BE[j]                                 #Calculate proton separation energy
            #print(S_p)
            if S_p>0.:                                        #Only continue if it's bound
                #print(S_p)
                for k in range(len(A)):
                    #print(A[k])
                    if Z[i]-2==Z[k] and N[i]==N[k]:           #Find isotope with Z-2
                        S_2p = BE[i]-BE[k]                    #Calculate 2-proton separation energy
                        #print(S_2p)
                        if S_2p < 0:                          #Only continue if it's unbound
                            #print(S_2p)
                            #print(Z[i],N[i],A[i],S_p,S_2p)   
                            sys.stdout=open('pdecay_out','a')                   #Print result to file
                            print("%s   %s   %s   %s"%(Z[i],N[i],S_p,S_2p))

sys.stdout.close()

                    
