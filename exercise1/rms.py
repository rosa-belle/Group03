import numpy as np
from matplotlib import pyplot as plt

#Read data and put columns in lists
data = np.loadtxt("rms13.dat",usecols=(0,1,2,3))

Z = data[:,0]
N = data[:,1]
rms = data[:,3]

diff = []
neutr = []
#print(Z)

for i in range(len(N)):                         
    if N[i]>N[i-1]:                           #Only use data from the same element
        delta = rms[i]-rms[i-1]               #calculate delta(rms)
        neutr.append(N[i])
        diff.append(delta)



plt.plot(neutr,diff,'ro')
plt.axis([0,160,-0.2,0.2])
plt.xlabel('N')
plt.ylabel('rms(N)-rms(N-1)')
plt.title('Brix-Kopfermann plot')
plt.savefig('rms.png')
plt.show()

