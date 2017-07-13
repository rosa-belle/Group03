import numpy as np
import sys
from matplotlib import pyplot as plt

#Read data file and create lists from columns
data = np.loadtxt("aud16.dat")

Z = data[:,0]
N = data[:,5]
BE = data[:,2]

#Set initial nucleus
n = 135
z = 88

for i in range(len(Z)):
    if N[i]==(n-8) and Z[i]==(z-6):                    #find daughter nucleus
        for j in range(len(Z)):
            if N[j]==n and Z[j]==z:                    #find mother nucleus
                Q = (BE[i] + BE[53] - BE[j])/1000      #calculate Q-value for 14C decay
                #print(BE[53])
                print(Q)                               #print Q in MeV
                T_log = Z[i]/np.sqrt(Q)                #calculate log(T1/2) of 14C decay
                print(T_log)
                #print(np.power(10,T_log))
                Tlit_a = 984960                        #half-life of alpha decay of 223Ra
                B = Tlit_a/np.power(10,T_log)          #branching ratio of 14C deacy
                print(B)
                sys.stdout=open('Q_results','w')
                print("Q_c = %s MeV" %Q)
                print("logT_c = %s" %T_log)
                print("B = %s" %B)
                sys.stdout.close()
                

                
