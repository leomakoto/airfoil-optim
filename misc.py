import matplotlib.pyplot as plt
import numpy as np
import numpy.random as rng
import scipy.special as scp
import os
import writer
import main
	
#Same as above, but random
def random_airfoil_eval(n,param,boundary,control,forces_control):
	#Random stuff
	par={'n':n, 'deg':int(rng.rand(1)*10), 'N1':rng.rand(1), 'N2':rng.rand(1), 'trail_gap':0.0}
	
	#Fix some variables
	par['N1']=0.5
	par['N2']=1.0
	par['deg']=2
	
	#Generate
	scale=rng.rand(1)
	Au=np.multiply(rng.rand(par['deg']+1), scale)
	Al=np.multiply(rng.rand(par['deg']+1), scale*0.2)
	
	#Write in feval.in (optional, just to keep track of it)
	input_string="{:.8f} {:.8f}".format(par['N1'], par['N2'])
	for val in Au: input_string=input_string+" {:.8f}".format(val)
	for val in Al: input_string=input_string+" {:.8f}".format(val)
	Mdict=open('eval/feval.in','w')
	Mdict.write(input_string)
	
	#Simulate airfoil
	D=main.airfoil_data(Au,Al,par,param,boundary,control,forces_control)

	#Print aerodynamical coefficients for last iteration
	Mdict=open('eval/feval.out','w')
	Mdict.write("{:.8f} {:.8f}".format(D['Cd'], D['Cl']))

#Same as airfoil_data + feval, but for a circle
def circle_airfoil_eval(n,param,boundary,control,forces_control):
	x_aux=np.linspace(0.0,1.0,n)
	y_aux=np.sqrt(0.25-(x_aux-0.5)**(2))
	x=np.concatenate((np.flip(x_aux),x_aux[1:n]))
	y=np.concatenate((np.flip(y_aux), -y_aux[1:n]))
	
	writer.write_blockMeshDict(x,y,param)
	
	dirlist=[x.path[2:] for x in os.scandir() if x.is_dir()]
	if '0' not in dirlist: os.system("mkdir 0")
	
	writer.write_boundaryCond(boundary)
	writer.write_controlDict(control)
	writer.write_forceCoeffs(forces_control)
	os.system("blockMesh > out_blockMesh.txt")
	os.system("checkMesh > out_checkMesh.txt")
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
	D = dict(zip(labels,data[-1]))
	
	#Print aerodynamical coefficients for last iteration
	Mdict=open('eval/feval.out','w')
	Mdict.write("{:.8f} {:.8f}".format(D['Cd'], D['Cl']))
	
#A test function
def test_dir(d,scale,n,param,boundary,control,forces_control,fvsolution):
	A=np.genfromtxt("eval/feval.in", delimiter=' ')
	n_cont=int((len(A)-2)/2)

	#Ignore first 2 components
	dA=d[2:]

	#Airfoil contour
	Au=A[2:(n_cont+2)]; Al=A[(n_cont+2):]

	#Parameters for CST airfoil parametrization
	par={'n':n, 'deg':(n_cont-1), 'N1':A[0], 'N2':A[1], 'trail_gap':0.0}
	
	Cd=[]; Cl=[]; CdCl=[]; delta=[]
	
	dA=np.multiply(dA,scale)
	
	Imax=30
	h=1.0/Imax
	#step=-0.5
	step=0.0
	
	gr, ind = plt.subplots(4,figsize=(6,8))
	gr.suptitle('Variation of Cd, Cl, and Cd/Cl along gradient desc. dir.')
	
	for i in range(Imax):
		Au = Au + np.multiply(dA[0:n_cont], step)
		Al = Al + np.multiply(dA[n_cont:], step)
		D = main.airfoil_data(Au,Al,par,param,boundary,control,forces_control,fvsolution)
		Cd.append(D['Cd'])
		Cl.append(D['Cl'])
		CdCl.append(D['Cd']/D['Cl'])
		
		[x,y]=main.cst(Au,Al,par)
		ind[0].plot(x,y,"0.{:d}".format(Imax**2-i**2))
		
		delta.append(step*scale)
		print("Step = {:d} --- Delta = {:.8f} --- Cd = {:.8f} --- Cl = {:.8f}".format(i, step*scale, D['Cd'], D['Cl']))
		step=step+h
		
	ind[1].plot(delta, Cd, color="blue")
	ind[2].plot(delta, Cl, color="orange")
	ind[3].plot(delta, CdCl, color="green")
	plt.show()

#Gradient of f, central differences, variables: Au, Al	
def grad_feval(n,param,boundary,control,forces_control,fvsolution,eps):
	A=np.genfromtxt("eval/feval.in", delimiter=' ')
	n_cont=int((len(A)-2)/2)

	#Airfoil contour for fixed N1 and N2
	x=A[2:]
	
	#Parameters for CST airfoil parametrization
	par={'n':n, 'deg':(n_cont-1), 'N1':A[0], 'N2':A[1], 'trail_gap':0.0}
	
	#Gradient of Cd and Cl in Au, Al
	gradCd=np.zeros(2*n_cont)
	gradCl=np.zeros(2*n_cont)
	gradf=np.zeros(2*n_cont)
	
	step=eps/2.0
	d=np.zeros(2*n_cont);
	
	for i in range(2*n_cont):
		d[i]=1.0
		x_plus = x + np.multiply(d, step)
		D = main.airfoil_data(x_plus[0:n_cont],x_plus[n_cont:],par,param,boundary,control,forces_control,fvsolution)
		Cd_plus=D['Cd']; Cl_plus=D['Cl']
		x_minus = x + np.multiply(d, -step)
		D = main.airfoil_data(x_minus[0:n_cont],x_minus[n_cont:],par,param,boundary,control,forces_control,fvsolution)
		Cd_minus=D['Cd']; Cl_minus=D['Cl']
		gradCd[i]=(Cd_plus-Cd_minus)/eps
		gradCl[i]=(Cl_plus-Cl_minus)/eps
		gradf[i]=(Cd_plus/Cl_plus - Cd_minus/Cl_minus)/eps
		d[i]=0.0
		print("{:d}-th partial derivative computed".format(i+1))
	
	gradCd=np.multiply(gradCd,1.0/np.linalg.norm(gradCd))
	gradCl=np.multiply(gradCl,1.0/np.linalg.norm(gradCl))
	
	#Write in file
	#input_string="{:.8f} {:.8f}".format(A[0], A[1])
	input_string="0.0 0.0"
	for val in gradf: input_string=input_string+" {:.8f}".format(val)
	Mdict=open('eval/gradfeval.out','w')
	Mdict.write(input_string)
	
	return [gradCd,gradCl,gradf]
		
def plot_airfoil(n,d,Au,Al,colscale):
	par={'n':n, 'deg':d, 'N1':A[0], 'N2':A[1], 'trail_gap':0.0}
	[x,y]=cst(Au,Al,par)
	plt.plot(x,y,"0.{:d}".format(9-i**2))
	plt.show()
	

