
# import python modules
import bempp.api
from bempp.api.file_interfaces import general_interface
import numpy as np
import sys

# input data

meshFile = "sphere2.msh"
freq    = 100                   # frequency
U0      = 0.01                  # amplitude
c       = 100                   # speed of sound
rho     = 14.1855               # density of gas
epsilon = 1E-5                  # solution accuracy
k       = freq * 2 * np.pi / c  # wave number
muD     = 1.0/k                 # numerical constant for combined BIE formulation
a       = 0.1                   # radius of sphere

# pressure amplitude                        
p0      = rho*c*1j*k*a*U0*np.exp(-1j*k*a) / (1j*k*a - 1)


print ("k = " + "%s" % k)

# output directory

import os
outputDirName = "output/"
if not os.path.exists(outputDirName):
    os.makedirs(outputDirName)

# mesh: unit sphere with standard shapes

reader = general_interface.import_grid(file_name = meshFile)
grid = reader

# time control

import time
time1 = time.time()

# spaces

piecewise_lin_space = bempp.api.function_space(grid, "P", 1)

domain_space = piecewise_lin_space 
range_space = piecewise_lin_space 
dual_to_range_space = piecewise_lin_space 

# operators

identity = bempp.api.operators.boundary.sparse.identity(
    domain_space, range_space, dual_to_range_space)
dlb = bempp.api.operators.boundary.helmholtz.double_layer(
    domain_space, range_space, dual_to_range_space,k)
slb = bempp.api.operators.boundary.helmholtz.single_layer(
    domain_space, range_space, dual_to_range_space,k)

# Dirichlet BC as function

def dir_data(x, n, domain_index, result):
    result[0] = p0
    
dirichlet_fun = bempp.api.GridFunction(piecewise_lin_space, fun = dir_data)

# CBIE

lhs = (.5*identity + dlb) - 1j * muD * slb

time2 = time.time()

print ("Prepared in {0} s".format(time2 - time1))

# solve

from bempp.api.linalg import gmres

it_count = 0 
def iteration_counter(x):
    global it_count
    it_count += 1
    
w_fun,info,resid, it_count = gmres(lhs, dirichlet_fun, tol=epsilon, return_residuals = True, return_iteration_count = True)

time3 = time.time()

print("Solved in {0} iterations, {1} s".format(it_count, time3 - time2))

# post-processing

print ("--------------")

exec(open("subscript/postProcess.py").read())
exec(open("subscript/analit.py").read())

# microphones
x = np.ndarray((6,), buffer=np.array([1.0,2.0,3.0,4.0,5.0,10.0]),dtype=np.float)
M = np.vstack((x.ravel(), 0*np.ones(x.size), 0*np.ones(x.size) ))

resM = abs(result(M))
analitM = analyticalSPU(U0,k,a,x)

for mic,num,an in zip(x,resM,analitM):
    print("Microphone: ({0},0,0), pressure amplitude: {1}, analytical: {2}".format(mic,num,an))
                    
# for gnuplot output
print ("--------------")    

R1 = 1
R2 = 5
h = 0.01

# optional postprocessing -- see graph/ and ouput/ folders
#plotPerRadius (R1,R2,h,outputDirName)
#analitR(U0,k,a,R1,R2,h,outputDirName)

#plotPerTime((2,0,0),0.0,0.01,1e-3,outputDirName)

#polarPlot(R2,outputDirName)
#analitPolar(U0,k,a,R2,h,outputDirName)

#matrixPlot(-10,-10,20,20,200,200,outputDirName)

print ("Total time = {0} s".format(time3 - time1))
print ("End")


