import matplotlib.pyplot as plt
import numpy as np
import os
		
#controlDict writer
def write_controlDict(control):
	Mdict=open('system/controlDict','w')
 
	#Header format
	header='FoamFile\n{{\n   version  2.0;\n   format  ascii;\n   class  dictionary;\n   object  controlDict;\n}}\n\napplication   simpleFoam;\nstartFrom   startTime;\nstartTime   {:d};\nstopAt   endTime;\nendTime   {:d};\ndeltaT   {:d};\nwriteControl   timeStep;\nwriteInterval   {:d};\npurgeWrite   2;\nwriteFormat   ascii;\nwritePrecision   6;\nwriteCompression   off;\ntimePrecision   6;\nrunTimeModifiable   true;\n\n#include "forceCoeffs.H"\n\nfunctions\n{{\n   '.format(control['t0'],control['tf'],control['dt'],control['writeint'])
	Mdict.write(header)
	
	#momErr and contErr format
	momErr='   momErr\n   {\n      type   mometumError;\n      libs   (fieldFunctionObjects);\n      executeControl   writeTime;\n      writeControl   writeTime;\n   }\n\n'
	contErr='   contErr\n   {\n      type   div;\n      libs   (fieldFunctionObjects);\n      field   phi;\n      executeControl   writeTime;\n      writeControl   writeTime;\n   }\n\n'
	
	#turbulence format
	turbulence='   turbulenceFields1\n   {\n      type   turbulenceFields;\n      libs   (fieldFunctionObjects);\n      fields\n      (\n         k\n         epsilon\n         nut\n         nuEff\n         R\n         devReff\n         L\n         I\n      );\n      executeControl   writeTime;\n      writeControl   writeTime;\n   }\n\n'
	
	#mag format
	mag='   mag1\n   {\n      type   mag;\n      libs   (fieldFunctionObjects);\n      field   turbulenceProperties:R;\n      result   magR;\n      executeControl   writeTime;\n      writeControl   writeTime;\n   }\n\n'
	
	Mdict.write(momErr)
	Mdict.write(contErr)
	Mdict.write(turbulence)
	Mdict.write(mag)
	Mdict.write("}\n")
	
def write_forceCoeffs(forces_control):
	Mdict=open('system/forceCoeffs.H','w')
	
	Mdict.write("functions\n{\n")
	
	CofR=[0.0, 0.0, 0.0]
	
	#forces format
	forces='   forces\n   {{\n      type   forces;\n      libs   ("libforces.so");\n      patches   (airfoil);\n      pName   p;\n      UName   U;\n      rho   rhoInf;\n      rhoInf   {:.3f};\n      CofR   ({:.2f} {:.2f} {:.2f});\n      writeControl   timeStep;\n      timeInterval   {:d};\n      log   {:s};\n   }}\n\n'.format(forces_control['rhoInf'],CofR[0],CofR[1],CofR[2],forces_control['timeint'],forces_control['log'])
	
	#Convert AoA to radians
	radAoA=(forces_control['AoA']*np.pi)/180.0
	liftDir=[-np.sin(radAoA), 0.0, np.cos(radAoA)];
	dragDir=[np.cos(radAoA), 0.0, np.sin(radAoA)];
	pitchAxis=[0.0, 1.0, 0.0];
	
	#forceCoeffs format
	forceCoeffs='   forceCoeffs\n   {{\n      type   forceCoeffs;\n      libs   ("libforces.so");\n      patches   (airfoil);\n      pName   p;\n      UName   U;\n      rho   rhoInf;\n      rhoInf   {:.3f};\n      CofR   ({:.2f} {:.2f} {:.2f});\n      liftDir   ({:.2f} {:.2f} {:.2f});\n      dragDir   ({:.2f} {:.2f} {:.2f});\n      pitchAxis   ({:.2f} {:.2f} {:.2f});\n      magUInf   {:.2f};\n      lRef   {:.2f};\n      Aref   {:.2f};\n      writeControl   timeStep;\n      timeInterval   {:d};\n      log   {:s};\n   }}\n\n'.format(forces_control['rhoInf'],CofR[0],CofR[1],CofR[2],liftDir[0],liftDir[1],liftDir[2],dragDir[0],dragDir[1],dragDir[2],pitchAxis[0],pitchAxis[1],pitchAxis[2],forces_control['magUInf'],forces_control['lRef'],forces_control['Aref'],forces_control['timeint'],forces_control['log'])
	
	Mdict.write(forces)
	Mdict.write(forceCoeffs)
	
	Mdict.write("}\n")
	
#boundary writer.
def write_boundaryCond(boundary):
	Mdict=open('0/U','w')
	
	#Dimensions: [mass (kg), space (m), time (s), temperature (K), concentration (mol), electric current (A), luminous intensity (cd)]
	#Example: [0 2 -1 0 0 0 0] means m^2 * s^(-1) = m^2/s
	
	#U format - m/s
	U='FoamFile\n{{\n   version  2.0;\n   format  ascii;\n   class  volVectorField;\n   object  U;\n}}\n\ndimensions   [0 1 -1 0 0 0 0];\n\ninternalField   uniform ({:.2f} {:.2f} {:.2f});\n\nboundaryField\n{{\n   freestream_front\n   {{\n      type   freestreamVelocity;\n      freestreamValue   $internalField;\n   }}\n\n   freestream_back\n   {{\n      type   freestreamVelocity;\n      freestreamValue   $internalField;\n   }}\n\n   airfoil\n   {{\n      type   noSlip;\n   }}\n\n   empty1\n   {{\n      type   empty;\n   }}\n\n   empty2\n   {{\n      type   empty;\n   }}\n}}\n\n'.format(boundary['u1'],boundary['u2'],boundary['u3'])
	Mdict.write(U)
	
	Mdict=open('0/p','w')
	
	#p format - m^2/s^2 - kinematic pressure (makes sense in incompressible case)
	p='FoamFile\n{{\n   version  2.0;\n   format  ascii;\n   class  volScalarField;\n   object  p;\n}}\n\ndimensions   [0 2 -2 0 0 0 0];\n\ninternalField   uniform {:.2f};\n\nboundaryField\n{{\n   freestream_front\n   {{\n      type   freestreamPressure;\n      freestreamValue   $internalField;\n   }}\n\n   freestream_back\n   {{\n      type   freestreamPressure;\n      freestreamValue   $internalField;\n   }}\n\n   airfoil\n   {{\n      type   zeroGradient;\n   }}\n\n   empty1\n   {{\n      type   empty;\n   }}\n\n   empty2\n   {{\n      type   empty;\n   }}\n}}\n\n'.format(boundary['p'])
	Mdict.write(p)
	
	Mdict=open('0/nut','w')
	
	#nut format - m^2/s
	nut='FoamFile\n{{\n   version  2.0;\n   format  ascii;\n   class  volScalarField;\n   object  nut;\n}}\n\ndimensions   [0 2 -1 0 0 0 0];\n\ninternalField   uniform {:.7f};\n\nboundaryField\n{{\n   freestream_front\n   {{\n      type   freestream;\n      freestreamValue   uniform {:.7f};\n   }}\n\n   freestream_back\n   {{\n      type   freestream;\n      freestreamValue   uniform {:.7f};\n   }}\n\n   airfoil\n   {{\n      type   nutUSpaldingWallFunction;\n      value   uniform 0.0;\n   }}\n\n   empty1\n   {{\n      type   empty;\n   }}\n\n   empty2\n   {{\n      type   empty;\n   }}\n}}\n\n'.format(boundary['nut'],boundary['nut'],boundary['nut'])
	Mdict.write(nut)
	
	Mdict=open('0/nuTilda','w')
	
	#nuTilda format - m^2/s
	nuTilda='FoamFile\n{{\n   version  2.0;\n   format  ascii;\n   class  volScalarField;\n   object  nuTilda;\n}}\n\ndimensions   [0 2 -1 0 0 0 0];\n\ninternalField   uniform {:.7f};\n\nboundaryField\n{{\n   freestream_front\n   {{\n      type   freestream;\n      freestreamValue   uniform {:.7f};\n   }}\n\n   freestream_back\n   {{\n      type   freestream;\n      freestreamValue   uniform {:.7f};\n   }}\n\n   airfoil\n   {{\n      type   fixedValue;\n      value   uniform 0.0;\n   }}\n\n   empty1\n   {{\n      type   empty;\n   }}\n\n   empty2\n   {{\n      type   empty;\n   }}\n}}\n\n'.format(boundary['nuTilda'],boundary['nuTilda'],boundary['nuTilda'])
	Mdict.write(nuTilda)
	
def write_fvSolution(fvsolution):
	Mdict=open('system/fvSolution','w')
 
	#Header format
	header='FoamFile\n{\n   version  2.0;\n   format  ascii;\n   class  dictionary;\n   object  fvSolution;\n}\n\nsolvers\n{\n'
	Mdict.write(header)
	
	#p format
	p='   p\n   {{\n      solver   GAMG;\n      smoother   GaussSeidel;\n      tolerance   1e-{:d};\n      relTol   1e-{:d};\n   }}\n\n'.format(fvsolution['ptol'],fvsolution['preltol'])
	
	#U format
	U='   U\n   {{\n      solver   GAMG;\n      smoother   GaussSeidel;\n      tolerance   1e-{:d};\n      relTol   1e-{:d};\n      nSweeps   {:d};\n   }}\n\n'.format(fvsolution['Utol'],fvsolution['Ureltol'],fvsolution['Usweeps'])
	
	#nuTilda format
	nuTilda='   nuTilda\n   {{\n      solver   GAMG;\n      smoother   GaussSeidel;\n      tolerance   1e-{:d};\n      relTol   1e-{:d};\n      nSweeps   {:d};\n   }}\n\n'.format(fvsolution['nuTildatol'],fvsolution['nuTildareltol'],fvsolution['nuTildasweeps'])
	
	Mdict.write(p)
	Mdict.write(U)
	Mdict.write(nuTilda)
	Mdict.write("}\n\n")
	
	SIMPLE='SIMPLE\n{{\n   nNonOrthogonalCorrectors   {:d};\n\n   residualControl\n   {{\n      p   1e-{:d};\n      U   1e-{:d};\n      nuTilda   1e-{:d};\n   }}\n}}\n\n'.format(fvsolution['nonortho'],fvsolution['pres'],fvsolution['Ures'],fvsolution['nuTildares'])
	Mdict.write(SIMPLE)
	
	relaxFactors='relaxationFactors\n{{\n   fields\n   {{\n      p   {:.3f};\n   }}\n\n   equations\n   {{\n      U   {:.3f};\n      nuTilda   {:.3f};\n   }}\n}}\n\n'.format(fvsolution['pfact'],fvsolution['Ufact'],fvsolution['nuTildafact'])
	Mdict.write(relaxFactors)

def test_run(cases,control,forces_control,boundary,fvsolution,AoAini,AoAend,n):	
	if n>1: AoA=np.linspace(AoAini,AoAend,n)
	else: AoA=AoAini
	
	Cd=np.zeros(n)
	Cl=np.zeros(n)
	
	for case in cases:
		#List all directories
		dirlist=[x.path[2:] for x in os.scandir() if x.is_dir()]
		
		#Replace current polyMesh
		if 'polyMesh' in dirlist: os.system("rm -r constant/polyMesh")
		os.system("cp -r airfoilMeshes/{:s}/polyMesh constant".format(case))
		
		write_controlDict(control)
		write_fvSolution(fvsolution)
		
		#Run different attack angles
		Res=open("test-results/results-{:s}.txt".format(case),'w')
		
		for j in range(n):
			if n>1: 
				print("---Running case: {:s} with AoA = {:.2f}---\n".format(case,AoA[j]))
				forces_control['AoA']=AoA[j]
				radAoA=(AoA[j]*np.pi)/180.0
				boundary['u1']=forces_control['magUInf']*np.cos(radAoA)
				boundary['u2']=0.0
				boundary['u3']=forces_control['magUInf']*np.sin(radAoA) 
				
				#Run simpleFoam
				write_forceCoeffs(forces_control)
				write_boundaryCond(boundary)
				#os.system("mpirun -np 4 simpleFoam -parallel > out/out_simpleFoam_test_{:s}_AoA{:.2f}.txt".format(case,AoA[j]))
				os.system("simpleFoam > out/out_simpleFoam_test_{:s}_AoA{:.2f}.txt".format(case,AoA[j]))
			else: 
				print("---Running case: {:s} with AoA = {:.2f}---\n".format(case,AoA))
				forces_control['AoA']=AoA
				radAoA=(AoA*np.pi)/180.0
				boundary['u1']=forces_control['magUInf']*np.cos(radAoA)
				boundary['u2']=0.0
				boundary['u3']=forces_control['magUInf']*np.sin(radAoA) 
				
				#Run simpleFoam
				write_forceCoeffs(forces_control)
				write_boundaryCond(boundary)
				#os.system("mpirun -np 4 simpleFoam -parallel > out/out_simpleFoam_test_{:s}_AoA{:.2f}.txt".format(case,AoA))
				os.system("simpleFoam > out/out_simpleFoam_test_{:s}_AoA{:.2f}.txt".format(case,AoA))
			
			#List all directories again (maybe unnecessary)
			dirlist=[x.path[2:] for x in os.scandir() if x.is_dir()]
			
			#List completed iterations
			it=[x for x in dirlist if x.isnumeric()]
			it=sorted([int(x) for x in it])
			
			#Save data and clean garbage
			#n_it=int(np.round((control['tf']-control['t0'])/control['writeint']))
			#for i in range(n_it): 
			#	if control['writeint']*(i+1)<=it[-1]: os.system("rm -r {:d}".format(control['writeint']*(i+1)))
			os.system("rm -r {:d}".format(it[1]))
			os.system("rm -r {:d}".format(it[2]))
			
			#Read data from the output of simpleFoam
			data=np.genfromtxt("postProcessing/forceCoeffs/0/coefficient_0.dat", delimiter='\t')
			labels=['Iteration', 'Cd', 'Cs', 'Cl', 'CmRoll', 'CmPitch', 'CmYaw', 'Cd(f)', 'Cd(r)', 'Cs(f)', 'Cs(r)', 'Cl(f)', 'Cl(r)']
			print(dict(zip(labels,data[-1])))
			print("\n\n")
			
			Cd[j]=data[-1][1]
			Cl[j]=data[-1][3]
			
			Res.write("{:.4f} {:.9f} {:.9f}\n".format(AoA[j],Cd[j],Cl[j]))
		Res.close()
	
def plot_results(cases):
	gr, ind = plt.subplots(2,figsize=(6,8))
	gr.suptitle('AoA vs CD and CL')
	Data=np.genfromtxt("experimental.txt", delimiter="  ")
	ind[1].scatter(Data[:,0], Data[:,1], color="black")
	ind[0].scatter(Data[:,0], Data[:,2], color="black")

	for case in cases:
		G=np.genfromtxt("test-results/results-{:s}.txt".format(case), delimiter=" ")
		ind[0].plot(G[:,0], G[:,1])
		ind[1].plot(G[:,0], G[:,2])

	plt.show()

#solver parameters
fvsolution={
	#tolerance (power 1e-X)
	'ptol':6, 
	'preltol':1,
	
	'Utol':6, 
	'Ureltol':1,
	'Usweeps':2,
	
	'nuTildatol':6, 
	'nuTildareltol':1,
	'nuTildasweeps':2,
	
	#Nonorthogonal correctors (how many times it corrects nonorthogonality)
	'nonortho':2, 
	
	#residual control (stopping criteria, power 1e-X)
	'pres':6,
	'Ures':6,
	'nuTildares':6,
	
	#relaxation factors
	'pfact':0.3,
	'Ufact':0.7,
	'nuTildafact':0.7,
}

#Maybe change gradSchemes to default cellMDLimited Gauss linear 0.5

#control parameters
control={
	#Initial time
	't0':0, 
	
	#Final time
	'tf':5000, 
	
	#Time step
	'dt':1,
	
	#Log
	'writeint':50
}

#forces control parameters (Mach 0.85 = magUInf 291.55)
forces_control={
	#Angle of attack
	'AoA':0.0,
	
	#Air density
	'rhoInf':1.225,
	
	#Free stream velocity in m/s
	'magUInf':51.4815,#291.55,#30.8, 
	
	#I don't know what this is
	'lRef':1.00, 
	'Aref':1.00, 
	'timeint':100, 
	
	#Log
	'log':'true'
}

#boundary parameters
radAoA=(forces_control['AoA']*np.pi)/180.0
kin_visc=8.58*10**(-6) #1.8*10**(-5)

boundary={
	#Velocity (3D)
	'u1':forces_control['magUInf']*np.cos(radAoA), 
	'u2':0.0,
	'u3':forces_control['magUInf']*np.sin(radAoA),
	
	#Pressure
	'p':0.0, 
	
	#Turbulence variables - ideal value = (3~5) * <kinematic viscosity>
	#<kinematic viscosity> of air in this case = 1.8e-5
	'nut':0.22*kin_visc,
	'nuTilda':4.0*kin_visc
}

#Test case list	
cases=[
	'113-33', 
	#'225-65', 
	#'449-129', 
	#'897-257', 
	#'1793-513'
]

#Main run - plot results
plot_results(cases)

