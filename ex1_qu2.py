#!/usr/bin/env python
# -*- coding: utf-8 ; ispell-local-dictionary: "english" -*-
"""
================================================================================
GET THE NEUTRON DRIP LINE FOR A RANGE OF Z VALUES

This program first defines the binding energy using the Bethe-Weizsäcker
formula without the pairing term, with mass and proton numbers A and Z as
variables. The following fixed values are assumed for the parameters of this
binding energy:
a1 = 15.49 MeV
a2 = 17.23 MeV
a3 = 0.697 MeV
a4 = 22.6 MeV
Then Z is varied over a specified range and for each value of Z, A is
progressively increased (starting from the value of Z + 1) until the neutron
separation energy (based on the binding energy) becomes strictly negative. In
fact, given the set of constants of the Bethe-Weizsächer formula, the neutron
separation energy should be always positive for first values of A greater or
equal to Z (verify it...), because a single neutron is always bound to one or
more protons.

USAGE (command line)
*****
$ ThisProgram zmin zmax

zmin: [integer] lower boundary (included) of the Z range. HAS TO BE GREATER THAN 0.
zmax: [integer] upper boundary (included) of the Z range.


04/07/2017
================================================================================
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker     # Needed to activate minor ticks.
import sys     # For command line file input.

nrange = []
zrange = []
drip_sn = []


# ==============================================================================
# Avoiding stupid cases
# ==============================================================================
zmin = int(sys.argv[1])
zmax = int(sys.argv[2])

if zmin <= 0:
    print "\n ERROR: zmin="+str(sys.argv[1])+""\
                         "\n        The lower boundary zmin cannot be lower or equal to zero."\
                         "\n        Exiting."
    sys.exit()
    
if zmin > zmax:
    print "\n ERROR: zmin="+str(sys.argv[1])+"; "+"zmax="+str(sys.argv[2])+""\
                         "\n        The lower boundary zmin cannot be greater than the upper boundary."\
                         "\n        Exiting."
    sys.exit()

    
# ==============================================================================
# Parameters of the liquid drop model binding energy (no pairing term)
# ==============================================================================
a1 = 15.49    # MeV
a2 = 17.23    # MeV
a3 = 0.697    # MeV
a4 = 22.6     # MeV


# ==============================================================================
# Define the liquid drop model binding energy
# ==============================================================================
def be(A, Z):
    vol = a1 * A     # Volume term.
    sur = - a2 * A**(2/3)     # Surface term.
    cou = - a3 * (Z*(Z-1)) * A**(-1/3)     # Coulomb term.
    asy = - a4 * (A-2*Z)**2 * A**(-1)     # Asymmetry term.
    return vol + sur + cou + asy


# ==============================================================================
# Variation over the Z range: find the dripping nuclei and store this
# information in a output file
# ==============================================================================
# Determine the Z range from input values:
zrange = range(zmin, zmax+1)

output_basename = "NDL_"+str(sys.argv[1])+"-"+str(sys.argv[2])

output = open(output_basename+".out","w")
output.write("zmin = "+str(sys.argv[1])+"\n"+"zmax = "+str(sys.argv[2])+"\n\n"+"N     Z     Neutron separation energy")

for z in zrange:
    a = z+1     # Starting point: one neutron.
    sn = be(a, z) - be(a-1, z)     # Separation energy of the (A,Z) nucleus.
    while sn >= 0:
        a += 1
        if a > 300:
            print "\n No positive neutron separation energy found before reaching A=300."\
                "\n Exiting."
            sys.exit()
        sn = be(a, z) - be(a-1, z)     # Separation energy of the (A,Z) nucleus.
    drip_sn.append( be(a-1, z) - be(a-2, z) )     # At this point, the last calculated separation erergy is the first (strictly) negative, therefore the last positive or null is the previous one.
    nrange.append((a-1)-z)     # Store the value of n (si comment on the previous line to know why 'a-1').
    output.write("\n"+str((a-1)-z)+"     "+str(z)+"     "+str(drip_sn[z-zmin]))     # Write in the output file. 'z-zmin' is the index of the last appended drip separation energy.

output.close()
        
# ==============================================================================
# Plot the neutron drip line
# ==============================================================================
#plt.style.use("systematics")

ax = plt.gca()

plt.scatter(nrange, zrange, marker="s", color=[0.15,0.35,1])

plt.title("Last neutron-bound nuclei for Z in ["+str(sys.argv[1])+","+str(sys.argv[2])+"]")
plt.xlabel("N")
plt.ylabel("Z")
# Minor ticks:
ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

plt.savefig(output_basename+".pdf")

plt.show()
plt.close()
