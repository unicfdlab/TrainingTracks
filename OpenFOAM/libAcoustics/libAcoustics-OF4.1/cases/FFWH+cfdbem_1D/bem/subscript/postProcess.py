def result (points):
    from bempp.api.operators.potential import helmholtz as helmholtz_potential

    slp_pot=helmholtz_potential.single_layer(piecewise_lin_space, points, k)
    dlp_pot=helmholtz_potential.double_layer(piecewise_lin_space, points, k)

    res =  1j * muD * slp_pot.evaluate(w_fun) - dlp_pot.evaluate(w_fun)
    
    return res[0]
    
def SPL(p):
    pRef = 2e-5
    return 20*np.log10(abs(p)/p0)
    
# write file
def writetxt(fileName,array):
    with open(fileName, 'wb') as f:
        np.savetxt(f, array, delimiter=' ', newline='\n', header='', footer='', comments='# ')


# plot along X-axis in given time point
def plotPerRadius (R1, R2, h, outputDirName):
    
    print ("Plot per radius...")
    
    coordx = np.arange(R1,R2+h,h,dtype=np.float)

    points = np.vstack((coordx.ravel(), 0*np.ones(coordx.size), 0*np.ones(coordx.size)))

    res = zip(coordx,abs(result(points)))

    fileName = outputDirName + "pressureR.txt"
    
    writetxt(fileName,res)
    
    
# sound pressure per time [t0, T + t0] in given point
def plotPerTime (point, t0, T, tau, outputDirName):

    print ("Plot per time...")
    
    p = np.vstack(point)
    res = result(p)
    
    fName =  outputDirName + "pressureTime.txt"

    fOutput = open(fName,"w")

    for t in np.arange(t0, T + t0 + tau, tau):
	
	cur1 = np.real(res[0] * np.exp (- 1j * 2 * np.pi * freq * (t)))
	fOutput.write('%s' % t + ' ' + '%s' % cur1)
	fOutput.write('\n')
    
    fOutput.close()

# sound pressure: polar plot in Oxz
def polarPlot(R, outputDirName):

    print ("Polar plot...")
    
    theta = np.arange(0,2*np.pi,0.01,dtype = np.float)

    x1 = R*np.cos(theta)
    z1 = R*np.sin(theta)
    points = np.vstack(( x1.ravel(), 0*np.ones(theta.size), z1.ravel() ))
    res = zip(theta,abs(result(points)))
    
    fileName = outputDirName + "polarR" + '%s'%R + ".txt"
    
    writetxt(fileName,res)

def matrixPlot (x1, y1, L1, L2, n1, n2, outputDirName):
    print ("Matrix plot...")
    plot_grid = np.mgrid[x1:x1+L1:n1*1j,y1:y1+L2:n2*1j]
    points = np.vstack((plot_grid[0].ravel(),np.zeros(plot_grid[0].size),plot_grid[1].ravel()))
    res = np.real(result(points))

    matrix = np.split(res,n1)
    
    fileName = outputDirName + "matrixPlot.txt"
    
    writetxt(fileName,matrix)         