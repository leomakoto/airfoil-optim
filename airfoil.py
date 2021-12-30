import matplotlib.pyplot as plt
import numpy as np
import numpy.random as rng
import scipy.special as scp
import os
import writer
import misc
import main
from timeit import default_timer as timer

#Create eval directory
dirlist=[x.path[2:] for x in os.scandir() if x.is_dir()]
if 'eval' not in dirlist: os.system("mkdir eval")

# :::::::::::::SOLVER:::::::::::::
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
	'pres':5,
	'Ures':5,
	'nuTildares':5,
	
	#relaxation factors
	'pfact':0.3,
	'Ufact':0.7,
	'nuTildafact':0.7,
}

# :::::::::::::CONTROL:::::::::::::
 
#control parameters
control={
	#Initial time
	't0':0, 
	
	#Final time (in stationary case this is iteration number)
	'tf':20000, 
	
	#Time step
	'dt':1,
	
	#Log
	'writeint':50,
	
	#Purge (how many iterations to save)
	'purge':2
}

# :::::::::::::FORCES:::::::::::::

#forces control parameters (Mach 0.85 = magUInf 291.55)
forces_control={
	#Angle of attack
	'AoA':5,
	
	#Air density
	'rhoInf':1.225,
	
	#Free stream velocity in m/s
	'magUInf':291.55,#30.8, 
	
	#I don't know what this is
	'lRef':1.00, 
	'Aref':1.00, 
	'timeint':100, 
	
	#Log
	'log':'true'
}

# :::::::::::::BOUNDARY:::::::::::::

#boundary parameters
radAoA=(forces_control['AoA']*np.pi)/180.0
kin_visc=1.8*10**(-5)

boundary={
	#Velocity (3D)
	'u1':forces_control['magUInf']*np.cos(radAoA), 
	'u2':forces_control['magUInf']*np.sin(radAoA), 
	'u3':0.0,
	
	#Pressure
	'p':0.0, 
	
	#Turbulence variables - ideal value = (3~5) * <kinematic viscosity>
	#<kinematic viscosity> of air in this case = 1e-5
	'nut':0.22*kin_visc,
	'nuTilda':3.0*kin_visc
}

# :::::::::::::MESH:::::::::::::

#Airfoil points
n=101
preset_scale=0.5

#blockMeshDict parameters:
param={
	#n = number of points in the airfoil contour
	'n':n,
	
	#scale = scale factor of the mesh
	'scale':1, 
	
	#radius = horizontal distance between the front of the airfoil and source of air, in particular, radius of the front semicircle
	'radius':15, 
	
	#dist_x_out = horizontal distance between the front of the airfoil and sink of air
	'dist_x_out':20, 
	
	#depth = z axis lenght
	'depth':0.3, 
	
	#angle = angle of attack
	'angle':forces_control['AoA'], 
	
	#exp_ratio = expansion ratio near airfoil (>1) close to 1 means that grid quality decays slower near airfoil
	'exp_ratio':1.05, 
	
	#bound_thick = boundary thickness (>0) close to 0 means denser mesh near airfoil (this is the spacement between several contours copies around the airfoil, like soundwaves)
	'bound_thick':0.5, 
	
	#first_thick = first layer thickness (total mystery)
	'first_thick':0.003, 
	
	#max_cell_size (in, out, inout) = maximum cell size (inlet, outlet, inlet x outlet) refinement level of grid at the circle and rectangle boundaries 
	'max_cell_size_in':1, 
	'max_cell_size_out':1, 
	'max_cell_size_in_out':1, 
	
	#cell_size (middle, trail, lead) = cell size at (middle, trailing edge, leading edge), defines refinement (>0) close to zero means more refined mesh
	'cell_size_middle':0.05, #Cell size around the airfoil
	'cell_size_trail':0.15, #Cell size near end of trail
	'cell_size_lead':0.003, #Cell size near nose of airfoil
	
	#sep_point = separating point from leading in (0.1, 0.9), appears to be somehow related with curvature of grid, 0.5 seems optimal.
	'sep_point':0.5
}

# :::::::::::::CALL:::::::::::::

#Run feval
main.feval(n,param,boundary,control,forces_control,fvsolution)

