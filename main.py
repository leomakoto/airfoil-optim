import matplotlib.pyplot as plt
import numpy as np
import numpy.random as rng
import scipy.special as scp
import os
import writer
from timeit import default_timer as timer

#CST airfoil parametrization.	
def cst(Au,Al,par):
	x=np.linspace(0.0,1.0,par['n'])
	y_upper=np.zeros(par['n'])
	y_lower=np.zeros(par['n'])
	for i in range(par['n']):
		for k in range(par['deg']+1):
			y_upper[i]=y_upper[i]+scp.binom(par['deg']+1,k)*(x[i]**k)*((1.0-x[i])**(par['deg']+1-k))*Au[k]
			y_lower[i]=y_lower[i]+scp.binom(par['deg']+1,k)*(x[i]**k)*((1.0-x[i])**(par['deg']+1-k))*Al[k]
		y_upper[i]=y_upper[i]*(x[i]**par['N1'])*((1.0-x[i])**par['N2'])+x[i]*par['trail_gap']
		y_lower[i]=y_lower[i]*(x[i]**par['N1'])*((1.0-x[i])**par['N2'])+x[i]*par['trail_gap']
	x_vals=np.concatenate((np.flip(x),x[1:par['n']]))
	y_vals=np.concatenate((np.flip(y_upper), -y_lower[1:par['n']]))
	#plt.plot(x_vals,y_vals)
	#plt.show()
	return [x_vals,y_vals]
	
def airfoil_data(Au,Al,par,param,boundary,control,forces_control,fvsolution):
	[x,y]=cst(Au,Al,par)
	writer.write_blockMeshDict(x,y,param)
	
	dirlist=[x.path[2:] for x in os.scandir() if x.is_dir()]
	if '0' not in dirlist: os.system("mkdir 0")
	
	writer.write_boundaryCond(boundary)
	writer.write_controlDict(control)
	writer.write_forceCoeffs(forces_control)
	writer.write_fvSolution(fvsolution)
	os.system("blockMesh > out_blockMesh.txt")
	os.system("checkMesh > out_checkMesh.txt")
	#os.system("checkMesh")
	#os.system("paraFoam")
	os.system("simpleFoam > out_simpleFoam.txt")
	
	#List all directories
	dirlist=[x.path[2:] for x in os.scandir() if x.is_dir()]
			
	#List completed iterations
	it=[x for x in dirlist if x.isnumeric()]
	it=sorted([int(x) for x in it])
			
	#Save data and clean garbage
	#n_it=int(np.round((control['tf']-control['t0'])/control['writeint']))
	#for i in range(n_it): 
	#	if control['writeint']*(i+1)<=it[-1]: os.system("rm -r {:d}".format(control['writeint']*(i+1)))
	for i in range(control['purge']): os.system("rm -r {:d}".format(it[i+1]))
	
	#Read data from the output of simpleFoam
	data=np.genfromtxt("postProcessing/forceCoeffs/0/coefficient_0.dat", delimiter='\t')
	labels=['Iteration', 'Cd', 'Cs', 'Cl', 'CmRoll', 'CmPitch', 'CmYaw', 'Cd(f)', 'Cd(r)', 'Cs(f)', 'Cs(r)', 'Cl(f)', 'Cl(r)']
	return dict(zip(labels,data[-1]))

#Evaluates Cd and Cl. Reads from feval.in and writes in feval.out. 
def feval(n,param,boundary,control,forces_control,fvsolution):
	#Read airfoil from txt. Order: N1 N2 Au Al
	A=np.genfromtxt("eval/feval.in", delimiter=' ')
	n_cont=int((len(A)-2)/2)
	
	#Very important observation about python's crop: A[1:1]=[] and A[1:2]=A[1]

	#Airfoil contour
	Au=A[2:(n_cont+2)]; Al=A[(n_cont+2):]

	#Parameters for CST airfoil parametrization
	par={'n':n, 'deg':(n_cont-1), 'N1':A[0], 'N2':A[1], 'trail_gap':0.0}

	#Simulate airfoil
	D=airfoil_data(Au,Al,par,param,boundary,control,forces_control,fvsolution)

	#Print aerodynamical coefficients for last iteration
	Mdict=open('eval/feval.out','w')
	Mdict.write("{:.8f} {:.8f}".format(D['Cd'], D['Cl']))
	
	#Print log - comment it if you want
	#print("\n.:Airfoil CST data:.\nN1 = {:.8f}\nN2 = {:.8f}".format(A[0],A[1]))
	#Au_string="Upper coeffs: "
	#for x in Au: Au_string += str(x) + " "
	#Al_string="Lower coeffs: "
	#for x in Al: Al_string += str(x) + " "
	#print(Au_string)
	#print(Al_string)
	print("\n.:Current point:.\n")
	print(A[0],A[1],Au,Al)
	print("\n.:Forces:.\nCd = {:.8f}\nCl = {:.8f}\n".format(D['Cd'], D['Cl']))
