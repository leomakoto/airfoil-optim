import numpy as np
import os

def pt3(a,b,c):
	return "   ({:.6f} {:.6f} {:.6f})\n".format(a,b,c)

#blockMeshDict writer
def write_blockMeshDict(x,y,param):
	Mdict=open('system/blockMeshDict','w')
 
	#Header format
	header="FoamFile\n{{\n   version  2.0;\n   format  ascii;\n   class  dictionary;\n   object  blockMeshDict;\n}}\n\nscale  {:.2f};\n\ngeometry\n{{\n}}\n\n".format(param['scale'])
	Mdict.write(header)
	
	#Vertices format
	Mdict.write("vertices\n(\n")
	Mdict.write(pt3(0,0,0)) #1
	Mdict.write(pt3(1,0,0)) #2
	Mdict.write(pt3(1,param['radius'],0)) #3
	Mdict.write(pt3(1-param['radius'],0,0)) #4
	Mdict.write(pt3(0,0,param['depth'])) #5
	Mdict.write(pt3(1,0,param['depth'])) #6
	Mdict.write(pt3(1,param['radius'],param['depth'])) #7 
	Mdict.write(pt3(1-param['radius'],0,param['depth'])) #8
	Mdict.write(pt3(1+param['dist_x_out'],np.sin((np.pi/180)*param['angle'])*(param['dist_x_out']+1),0)) #9
	Mdict.write(pt3(1+param['dist_x_out'],param['radius'],0)) #10
	Mdict.write(pt3(1+param['dist_x_out'],np.sin((np.pi/180)*param['angle'])*(param['dist_x_out']+1),param['depth'])) #11
	Mdict.write(pt3(1+param['dist_x_out'],param['radius'],param['depth'])) #12
	Mdict.write(pt3(1,-param['radius'],0)) #13
	Mdict.write(pt3(1,-param['radius'],param['depth'])) #14
	Mdict.write(pt3(1+param['dist_x_out'],-param['radius'],0)) #15
	Mdict.write(pt3(1+param['dist_x_out'],-param['radius'],param['depth'])) #16
	Mdict.write(pt3(1,0,0)) #17
	Mdict.write(pt3(1,0,param['depth'])) #18
	Mdict.write(");\n\n")

	#Blocks format	
	
	#Number of mesh on boundary layer
	N10=int(np.round(np.log10(1+(param['exp_ratio']-1)*param['bound_thick']/param['first_thick'])/np.log10(param['exp_ratio']), 0)); 
	
	#Expansion ratio on boundary layer
	O10=param['exp_ratio']**N10
	
	#Expansion ratio out of boundary layer
	O13=param['max_cell_size_in']/(param['first_thick']*param['exp_ratio']**N10)
	
	#Expansion ratio at tail
	O16=param['max_cell_size_out']/param['cell_size_trail']
	
	#Expansion ratio 1
	H11=(param['first_thick']*param['exp_ratio']**N10)/param['radius']*(O13-1)+1
	
	#Number of mesh out of boundary layer
	N13=int(np.round(np.log10(O13)/np.log10(H11), 0)); 
	
	#Expansion ratio 2
	H13=param['cell_size_trail']/param['dist_x_out']*(O16-1)+1
	
	#Number of mesh at tail
	N16=int(np.round(np.log10(O16)/np.log10(H13), 0))
	
	#Expansion ratio at outlet
	O18=param['max_cell_size_in_out']*O13/param['max_cell_size_out']*(N13+N10)/N13
	
	#Expansion ratio at leading
	O22=param['cell_size_middle']/param['cell_size_lead']
	
	#Expansion ratio at trailing
	O24=param['cell_size_middle']/param['cell_size_trail']
	
	#Expansion ratio 3
	H15=param['cell_size_lead']/param['sep_point']*(O22-1)+1
	
	#Expansion ratio 4
	H17=param['cell_size_trail']/(1-param['cell_size_middle'])*(O24-1)+1
	
	#Number of mesh at leading
	O21=int(np.round(np.log10(O22)/np.log10(H15), 0))
	
	#Number of mesh at trailing
	O23=int(np.round(np.log10(O24)/np.log10(H17), 0)); 
	
	#Expansion ratio inlet 1
	O28=param['cell_size_trail']/param['max_cell_size_in']; 
	
	Mdict.write("blocks\n(\n   hex  (0 1 2 3 4 5 6 7)  ({:d} {:d} 1)\n   edgeGrading\n   (\n   (\n".format(O21+O23, N10+N13))
	Mdict.write(pt3(param['sep_point'],O21/(O21+O23),O22))
	Mdict.write(pt3(param['sep_point'],1-O21/(O21+O23),1/O24))
	Mdict.write(")\n")
	Mdict.write("{:.6f} {:.6f}\n(\n".format(O28,O28))
	Mdict.write(pt3(param['sep_point'],O21/(O21+O23),O22))
	Mdict.write(pt3(1-param['sep_point'],1-O21/(O21+O23),1/O24))
	Mdict.write(")\n(\n")
	Mdict.write(pt3(param['bound_thick']/param['radius'],N10/(N13+N10),O10))
	Mdict.write(pt3(1-param['bound_thick']/param['radius'],1-N10/(N13+N10),O13))
	Mdict.write(")\n(\n")
	Mdict.write(pt3(param['bound_thick']/param['radius'],N10/(N13+N10),O10)) 
	Mdict.write(pt3(1-param['bound_thick']/param['radius'],1-N10/(N13+N10),O13))
	Mdict.write(")\n(\n")
	Mdict.write(pt3(param['bound_thick']/param['radius'],N10/(N13+N10),O10)) 
	Mdict.write(pt3(1-param['bound_thick']/param['radius'],1-N10/(N13+N10),O13))
	Mdict.write(")\n(\n")
	Mdict.write(pt3(param['bound_thick']/param['radius'],N10/(N13+N10),O10)) 
	Mdict.write(pt3(1-param['bound_thick']/param['radius'],1-N10/(N13+N10),O13))
	Mdict.write(")\n\n\n")
	Mdict.write("   {:.6f} {:.6f} {:.6f} {:.6f} \n".format(1,1,1,1))  
	Mdict.write(") \n \n hex  (1 8 9 2 5 10 11 6) ( {:d}  {:d}  1)\n edgeGrading \n( \n \n ".format(N16,N13+N10))
	Mdict.write("{:.6f} {:.6f} {:.6f} {:.6f} ".format(O16,O16,O16,O16))
	Mdict.write(" \n( \n")
	Mdict.write(pt3(param['bound_thick']/param['radius'],N10/(N13+N10),O10))
	Mdict.write(pt3(1-param['bound_thick']/param['radius'],1-N10/(N13+N10),O13))
	Mdict.write("\n) \n")
	Mdict.write("{:.6f} {:.6f}".format(O18,O18))
	Mdict.write("\n( \n")
	Mdict.write(pt3(param['bound_thick']/param['radius'],N10/(N13+N10),O10))
	Mdict.write(pt3(1-param['bound_thick']/param['radius'],1-N10/(N13+N10),O13))
	Mdict.write("\n) \n \n \n ")
	Mdict.write(" 1  1  1  1 ")
	Mdict.write("\n) \n \n hex  (3 12 16 0 7 13 17 4)  ( {:d} {:d} {:d} )  \n edgeGrading \n(\n \n {:.6f} \n( \n".format(O21+O23,N10+N13,1,O28))
	Mdict.write(pt3(param['sep_point'],O21/(O21+O23),O22))   
	Mdict.write(pt3(1-param['sep_point'],1-O21/(O21+O23),1/O24))
	Mdict.write(")\n(\n")
	Mdict.write(pt3(param['sep_point'],O21/(O21+O23),O22))
	Mdict.write(pt3(1-param['sep_point'],1-O21/(O21+O23),1/O24))
	Mdict.write(") \n {:.6f} \n  \n(\n".format(O28))
	Mdict.write(pt3(1-param['bound_thick']/param['radius'],1-N10/(N13+N10),1/O13))
	Mdict.write(pt3(param['bound_thick']/param['radius'],N10/(N13+N10),1/O10))
	Mdict.write(") \n(\n")
	Mdict.write(pt3(1-param['bound_thick']/param['radius'],1-N10/(N13+N10),1/O13))
	Mdict.write(pt3(param['bound_thick']/param['radius'],N10/(N13+N10),1/O10))
	Mdict.write(") \n( \n")    
	Mdict.write(pt3(1-param['bound_thick']/param['radius'],1-N10/(N13+N10),1/O13))
	Mdict.write(pt3(param['bound_thick']/param['radius'],N10/(N13+N10),1/O10))
	Mdict.write(") \n( \n")    
	Mdict.write(pt3(1-param['bound_thick']/param['radius'],1-N10/(N13+N10),1/O13))
	Mdict.write(pt3(param['bound_thick']/param['radius'],N10/(N13+N10),1/O10))
	Mdict.write(") \n \n \n 1 1 1 1 \n ) \n \n \n \n \n \n ") 
	Mdict.write("hex  (12  14  8  16  13  15  10  17)  ( {:d}  {:d}  1) \n edgeGrading \n( \n {:.6f}  {:.6f}  {:.6f}  {:.6f} \n  \n(\n".format(N16,N10+N13,O16,O16,O16,O16))
	Mdict.write(pt3(1-param['bound_thick']/param['radius'],1-N10/(N13+N10),1/O13)) 
	Mdict.write(pt3(param['bound_thick']/param['radius'],N10/(N13+N10),1/O10)) 
	Mdict.write(") \n {:.6f}  {:.6f}".format(1/O18,1/O18))
	Mdict.write("\n(\n")
	Mdict.write(pt3(1-param['bound_thick']/param['radius'],1-N10/(N13+N10),1/O13))
	Mdict.write(pt3(param['bound_thick']/param['radius'],N10/(N13+N10),1/O10))
	Mdict.write(")\n\n 1  1  1  1 \n) \n); \n \n \n")
	Mdict.write("edges \n(\n    arc  3 2 ( {:.6f} {:.6f} {:.6f} )\n".format(1-param['radius']*np.sin(np.pi/4),param['radius']*np.sin(np.pi/4),0))
	Mdict.write("    arc  7 6 ( {:.6f} {:.6f} {:.6f} )\n".format(1-param['radius']*np.sin(np.pi/4),param['radius']*np.sin(np.pi/4),param['depth']))
	Mdict.write("\nspline  1  0 \n(\n")
	for i in range(param['n']-2):
		Mdict.write(pt3(x[i+1],y[i+1],0))
	Mdict.write(")\n")
	Mdict.write("\nspline  5  4 \n(\n")
	for i in range(param['n']-2):
		Mdict.write(pt3(x[i+1],y[i+1],param['depth']))
	Mdict.write(")\n\n")
	Mdict.write("    arc  3 12 ( {:.6f} {:.6f} {:.6f} )\n".format(1-param['radius']*np.sin(np.pi/4),-param['radius']*np.sin(np.pi/4),0))
	Mdict.write("    arc  7 13 ( {:.6f} {:.6f} {:.6f} )\n".format(1-param['radius']*np.sin(np.pi/4),-param['radius']*np.sin(np.pi/4),param['depth']))
	Mdict.write("\nspline  0  16 \n(\n")
	for i in range(param['n']-2):
		Mdict.write(pt3(x[i+param['n']],y[i+param['n']],0))
	Mdict.write(")\n")
	Mdict.write("\nspline  4  17 \n(\n")
	for i in range(param['n']-2):
		Mdict.write(pt3(x[i+param['n']],y[i+param['n']],param['depth']))
	Mdict.write(")\n);\n\n")
	Mdict.write("faces\n(\n\n);\n\nfaces\n(\n\n);\n\n\ndefaultPatch\n{\n   name frontAndBack;\n   type empty;\n}\n\n")
	Mdict.write("boundary\n(\ninlet\n   {\n    type patch;\n    faces\n    (\n     (9 2 6 11)\n     (2 3 7 6)\n     (3 12 13 7)\n     (12 15 14 13)\n    );\n   }\n\n")
	Mdict.write("outlet\n   {\n    type patch;\n    faces\n    (\n     (8 9 10 11)\n     (15 8 10 14)\n    );\n   }\n\n")
	Mdict.write("walls\n   {\n    type wall;\n    faces\n    (\n     (0 1 5 4)\n     (0 4 17 16)\n    );\n   }\n\n")
	Mdict.write("interface1\n   {\n    type patch;\n    faces\n    (\n     (1 8 10 5)\n    );\n   }\n\n")
	Mdict.write("interface2\n   {\n    type patch;\n    faces\n    (\n     (16 17 10 8)\n    );\n   }\n);\n\n")
	Mdict.write("mergePatchPairs\n(\n   (interface1 interface2)\n);")

#controlDict writer
def write_controlDict(control):
	Mdict=open('system/controlDict','w')
 
	#Header format
	header='FoamFile\n{{\n   version  2.0;\n   format  ascii;\n   class  dictionary;\n   object  controlDict;\n}}\n\napplication   simpleFoam;\nstartFrom   startTime;\nstartTime   {:d};\nstopAt   endTime;\nendTime   {:d};\ndeltaT   {:d};\nwriteControl   timeStep;\nwriteInterval   {:d};\npurgeWrite   {:d};\nwriteFormat   ascii;\nwritePrecision   6;\nwriteCompression   off;\ntimePrecision   6;\nrunTimeModifiable   true;\n\n#include "forceCoeffs.H"\n\nfunctions\n{{\n   '.format(control['t0'],control['tf'],control['dt'],control['writeint'],control['purge'])
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
	
#forceCoeffs writer
def write_forceCoeffs(forces_control):
	Mdict=open('system/forceCoeffs.H','w')
	
	Mdict.write("functions\n{\n")
	
	CofR=[0.0, 0.0, 0.0]
	
	#forces format
	forces='   forces\n   {{\n      type   forces;\n      libs   ("libforces.so");\n      patches   (walls);\n      pName   p;\n      UName   U;\n      rho   rhoInf;\n      rhoInf   {:.3f};\n      CofR   ({:.2f} {:.2f} {:.2f});\n      writeControl   timeStep;\n      timeInterval   {:d};\n      log   {:s};\n   }}\n\n'.format(forces_control['rhoInf'],CofR[0],CofR[1],CofR[2],forces_control['timeint'],forces_control['log'])
	
	#Convert AoA to radians
	radAoA=(forces_control['AoA']*np.pi)/180.0
	liftDir=[-np.sin(radAoA), np.cos(radAoA), 0.0];
	dragDir=[np.cos(radAoA), np.sin(radAoA), 0.0];
	pitchAxis=[0.0, 0.0, 1.0];
	
	#forceCoeffs format
	forceCoeffs='   forceCoeffs\n   {{\n      type   forceCoeffs;\n      libs   ("libforces.so");\n      patches   (walls);\n      pName   p;\n      UName   U;\n      rho   rhoInf;\n      rhoInf   {:.3f};\n      CofR   ({:.2f} {:.2f} {:.2f});\n      liftDir   ({:.2f} {:.2f} {:.2f});\n      dragDir   ({:.2f} {:.2f} {:.2f});\n      pitchAxis   ({:.2f} {:.2f} {:.2f});\n      magUInf   {:.2f};\n      lRef   {:.2f};\n      Aref   {:.2f};\n      writeControl   timeStep;\n      timeInterval   {:d};\n      log   {:s};\n   }}\n\n'.format(forces_control['rhoInf'],CofR[0],CofR[1],CofR[2],liftDir[0],liftDir[1],liftDir[2],dragDir[0],dragDir[1],dragDir[2],pitchAxis[0],pitchAxis[1],pitchAxis[2],forces_control['magUInf'],forces_control['lRef'],forces_control['Aref'],forces_control['timeint'],forces_control['log'])
	
	Mdict.write(forces)
	Mdict.write(forceCoeffs)
	
	Mdict.write("}\n")

#boundary writer.
def write_boundaryCond(boundary):
	Mdict=open('0/U','w')
	
	#Dimensions: [mass (kg), space (m), time (s), temperature (K), concentration (mol), electric current (A), luminous intensity (cd)]
	#Example: [0 2 -1 0 0 0 0] means m^2 * s^(-1) = m^2/s
	
	#U format - m/s
	U='FoamFile\n{{\n   version  2.0;\n   format  ascii;\n   class  volVectorField;\n   object  U;\n}}\n\ndimensions   [0 1 -1 0 0 0 0];\n\ninternalField   uniform ({:.2f} {:.2f} {:.2f});\n\nboundaryField\n{{\n   inlet\n   {{\n      type   freestreamVelocity;\n      freestreamValue   $internalField;\n   }}\n\n   outlet\n   {{\n      type   freestreamVelocity;\n      freestreamValue   $internalField;\n   }}\n\n   walls\n   {{\n      type   noSlip;\n   }}\n\n   frontAndBack\n   {{\n      type   empty;\n   }}\n}}\n\n'.format(boundary['u1'],boundary['u2'],boundary['u3'])
	Mdict.write(U)
	
	Mdict=open('0/p','w')
	
	#p format - m^2/s^2 - kinematic pressure (makes sense in incompressible case)
	p='FoamFile\n{{\n   version  2.0;\n   format  ascii;\n   class  volScalarField;\n   object  p;\n}}\n\ndimensions   [0 2 -2 0 0 0 0];\n\ninternalField   uniform {:.2f};\n\nboundaryField\n{{\n   inlet\n   {{\n      type   freestreamPressure;\n      freestreamValue   $internalField;\n   }}\n\n   outlet\n   {{\n      type   freestreamPressure;\n      freestreamValue   $internalField;\n   }}\n\n   walls\n   {{\n      type   zeroGradient;\n   }}\n\n   frontAndBack\n   {{\n      type   empty;\n   }}\n}}\n\n'.format(boundary['p'])
	Mdict.write(p)
	
	Mdict=open('0/nut','w')
	
	#nut format - m^2/s
	nut='FoamFile\n{{\n   version  2.0;\n   format  ascii;\n   class  volScalarField;\n   object  nut;\n}}\n\ndimensions   [0 2 -1 0 0 0 0];\n\ninternalField   uniform {:.7f};\n\nboundaryField\n{{\n   inlet\n   {{\n      type   freestream;\n      freestreamValue   uniform {:.7f};\n   }}\n\n   outlet\n   {{\n      type   freestream;\n      freestreamValue   uniform {:.7f};\n   }}\n\n   walls\n   {{\n      type   nutUSpaldingWallFunction;\n      value   uniform 0.0;\n   }}\n\n   frontAndBack\n   {{\n      type   empty;\n   }}\n}}\n\n'.format(boundary['nut'],boundary['nut'],boundary['nut'])
	Mdict.write(nut)
	
	Mdict=open('0/nuTilda','w')
	
	#nuTilda format - m^2/s
	nuTilda='FoamFile\n{{\n   version  2.0;\n   format  ascii;\n   class  volScalarField;\n   object  nuTilda;\n}}\n\ndimensions   [0 2 -1 0 0 0 0];\n\ninternalField   uniform {:.7f};\n\nboundaryField\n{{\n   inlet\n   {{\n      type   freestream;\n      freestreamValue   uniform {:.7f};\n   }}\n\n   outlet\n   {{\n      type   freestream;\n      freestreamValue   uniform {:.7f};\n   }}\n\n   walls\n   {{\n      type   fixedValue;\n      value   uniform 0.0;\n   }}\n\n   frontAndBack\n   {{\n      type   empty;\n   }}\n}}\n\n'.format(boundary['nuTilda'],boundary['nuTilda'],boundary['nuTilda'])
	Mdict.write(nuTilda)
	
#fvSolution writer
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
	
