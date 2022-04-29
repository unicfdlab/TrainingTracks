"""
Analytical solution for pulsating sphere
----------------------------------------

Sphere vibrates along Ox axis:

u(r) = p0 / r * exp(1j*k*r),

p0   = rho*c*1j*k*a*a*U0*exp(-1j*k*a) / (1j*k*a - 1)

p0 	--- amplitude of pressure pulsation
U0 	--- amplitude of velocity pulsation
k  	--- wave number
a 	--- radius of sphere
rho	--- air density

"""

import numpy as np
import sys

def analyticalSPU(U0,k,a,r):
    p = rho*c*1j*k*a*a*U0*np.exp(-1j*k*a) / (1j*k*a - 1)
    return abs(p * np.exp(-1j*k*r)) / r
    
def analyticalSP(p0,a,r):
    return p0*a / r
    
#------------------------------------------------------
    
def analitR(U0,k,a,R1,R2,h,outputDirName):

    print ("Analytical solution per radius...")

    points = np.arange(R1,R2+h,h,dtype=float)
    res = analyticalSPU(U0,k,a,points)

    # write

    fName =  outputDirName + "analitR.txt"

    fOutput = open(fName,"w")
    
    i = 0
    for num in res:
        fOutput.write('%s' % points[i] + ' ' + '%s' % num)
        fOutput.write('\n')
        i = i + 1

    fOutput.close()

def analitPolar(U0,k,a,R,h,outputDirName):

    print ("Polar analytical solution...")

    theta = np.arange(0,2*np.pi,h,dtype = np.float)
    
    res = [analyticalSPU(U0,k,a,R)]*len(theta)

    # write

    fName =  outputDirName + "analitPolar.txt"

    fOutput = open(fName,"w")
    
    i = 0
    for num in res:
        fOutput.write('%s' % theta[i] + ' ' + '%s' % abs(num))
        fOutput.write('\n')
        i = i + 1

    fOutput.close()
                