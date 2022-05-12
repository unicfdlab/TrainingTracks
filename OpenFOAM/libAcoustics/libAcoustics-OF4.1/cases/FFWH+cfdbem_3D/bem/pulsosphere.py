# import python modules
import sys
import bempp.api
import numpy as np
from cfdbem.file_interfaces import FileReader

# input data
c = int(sys.argv[1])
freq = int(sys.argv[2])
dataName = sys.argv[3]

epsilon = 1E-5

k = freq * 2 * np.pi / c;
print ("k = {0}".format(k))


U0  = 0.01 # amplitude
c   = 100 # speed of sound
rho = 14.1855 # density of gas
a   = 0.1 # radius of initial sphere


# numerical constant for combined formulation
muD = 1.0/k

# mesh
reader = FileReader(file_name = dataName)

grid = reader.grid

import time
time1 = time.time()

# spaces
piecewise_lin_space = bempp.api.function_space(grid, "P", 1)
piecewise_const_space = bempp.api.function_space(grid, "DP", 0)

domain_space = piecewise_lin_space 
range_space = piecewise_lin_space 
dual_to_range_space = piecewise_lin_space 


# operators
identity = bempp.api.operators.boundary.sparse.identity(
    domain_space, range_space, dual_to_range_space)
dlp = bempp.api.operators.boundary.helmholtz.double_layer(
    domain_space, range_space, dual_to_range_space,k)
slp = bempp.api.operators.boundary.helmholtz.single_layer(
    domain_space, range_space, dual_to_range_space,k)

        
#print ("ok")

# dirichlet BC
dirichlet_fun = bempp.api.GridFunction(piecewise_lin_space, coefficients = reader.node_data)

# equation in combined formulation!!!!!
lhs = (.5*identity + dlp) - 1j * muD * slp

time2 = time.time()

print ("prepared to ", time2 - time1, " s")


# solve
from bempp.api.linalg import gmres

it_count = 0
def iteration_counter(x):
    global it_count
    it_count += 1
    
w_fun,info,resid, it_count = gmres(lhs, dirichlet_fun, tol=epsilon, return_residuals = True, return_iteration_count = True)
print("The linear system was solved in {0} iterations".format(it_count))

time3 = time.time()

print ("solved to ", time3 - time2, " s")



# compute and output result

execfile("bem/subscript/analit.py")
execfile("bem/subscript/postProcess.py")
    
# microphones
print ("------------")

x = np.ndarray((6,), buffer=np.array([1.0,2.0,3.0,4.0,5.0,10.0]),dtype=np.float)
M = np.vstack((x.ravel(), 0*np.ones(x.size), 0*np.ones(x.size) ))

resM = abs(result(M))
analitM = analyticalSPU(U0,k,a,x)

for mic,num,an in zip(x,resM,analitM):
    print("Mic: (0,0,{0}), pressure amplitude: {1}, analytical solution: {2}".format(mic,num,an))
    
    
    
#plotPerTime((0,0,1),0,0.1,0.0001,"pressureTime1.dat")
#plotPerTime((0,0,2),0,0.1,0.0001,"pressureTime2.dat")
#plotPerTime((0,0,3),0,0.1,0.0001,"pressureTime3.dat")
#plotPerTime((0,0,4),0,0.1,0.0001,"pressureTime4.dat")
#plotPerTime((0,0,5),0,0.1,0.0001,"pressureTime5.dat")




