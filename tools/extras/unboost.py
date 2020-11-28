'''Transform boosted frame data to the lab frame.
Creates a new data file and adds to the existing metadata file.
Please note this routine still has hard coded parameters.'''
import sys
import numpy as np
import scipy.interpolate
import json
import twutils.plot as twplot

if len(sys.argv)<6:
	print('Usage: python unboost.py <gamma> <x> <y> <file_base> <component> coup=<coupled component> vb=<vb> vu=<vu>')
	print('<gamma> = gamma of the boosted frame (if negative reverse the boost direction)')
	print('<x>,<y> = transverse indices of the lineout to be transformed')
	print('<file_base> = expect filenames like <file_base>_<component>.npy')
	print('<component> = the component to transform')
	print('<coupled component> = Some unboosts are done in pairs.  This is the other component in the pair.')
	print('  for pairs, <component> is the time-like one, examples are (phi,Az) and (Ex,By).')
	print('<vb>,<vu> = galilean velocity in the boosted and lab frames, respectively.')
	exit(0)

Nlab = (512,2048)
g = np.double(sys.argv[1])
x = np.int(sys.argv[2])
y = np.int(sys.argv[3])
filename = sys.argv[4]
tcomp = sys.argv[5]
zcomp = None
vb = 0.0
vu = 0.0
valid_keys = ['coup','vb','vu']
for pair in sys.argv[6:]:
	key = pair.split('=')[0]
	val = pair.split('=')[1]
	if key not in valid_keys:
		raise KeyError('<'+key+'> is not a valid argument key.')
	if key=='coup':
		zcomp = val
	if key=='vb':
		vb = float(val)
	if key=='vu':
		vu = float(val)

gb = np.sign(g)*(g**2-1)**0.5
g = np.abs(g)

def Unboost(F,boosted_nodes,vb,unboosted_nodes,vu):
	'''Unboost in two spacetime dimensions.

	:param numpy.ndarray F: 2-tuple of fields in boosted frame, each field having shape (Nt,Nz), e.g., (Ax,None) or (Ex,By).
	:param tuple boosted_nodes: tuple of arrays (t',z'), where t' has shape (Nt,) and z' has shape (Nz,)
	:param float vb: Interpret z' as z'- vb*t'.
	:param tuple unboosted_nodes: tuple of arrays (t,z), where the new grid has shape (Ut,Uz), and t.shape = (Ut,) and z.shape = (Uz,)
	:param bool vu: Interpret z as z - vu*t.

	:return: unboosted field.  If F[1]==None, only coordinates are transformed.  Otherwise the value is transformed as g*F[0] + gb*F[1].
	:rtype: numpy.ndarray'''
	# Strategy is to interpolate on the uniform boosted grid.
	# We merely have to put the unboosted nodes in boosted frame coordinates, then interpolate.
	# In order to avoid expensive triangulation we demand the boosted grid be rectangular.
	tp = boosted_nodes[0]
	zp = boosted_nodes[1]
	U = (unboosted_nodes[0].shape[0] , unboosted_nodes[1].shape[0])
	x0 = np.zeros(U+(2,))
	x0[...,0] = unboosted_nodes[0][:,np.newaxis]
	x0[...,1] = unboosted_nodes[1][np.newaxis,:]
	x = np.copy(x0)
	# Boost the unboosted nodes
	x0[...,1] = x0[...,1] + vu*x0[...,0]
	x[...,0] = g*x0[...,0] - gb*x0[...,1]
	x[...,1] = g*x0[...,1] - gb*x0[...,0]
	x[...,1] = x[...,1] - vb*x[...,0]
	# Report the new domain as it appears in the boosted domain
	print('Earliest t sampled in boosted domain',np.min(x[...,0]))
	print('Latest t sampled in boosted domain',np.max(x[...,0]))
	print('Lowest z sampled in boosted domain',np.min(x[...,1]))
	print('Highest z sampled in boosted domain',np.max(x[...,1]))
	# unboost the field
	if type(F[1])==type(None):
		A = F[0]
	else:
		A = g*F[0] + gb*F[1]
	# Now we have the unboosted values and nodes on the boosted grid.
	# The interpolation moves this over to the unboosted grid.
	func = scipy.interpolate.RegularGridInterpolator((tp,zp),A,bounds_error=False,fill_value=0.0)
	pts = x.reshape(-1,2)
	return func(pts).reshape(U)

print('Get boosted frame data...')
p_ex = twplot.plotter(filename+'_'+tcomp+'.npy')
p_ex.display_info()
Ex,plot_dict = p_ex.falsecolor2d('ztxy',(x,y),dyn_range=0.0)
if zcomp!=None:
	p_by = twplot.plotter(filename+'_'+zcomp+'.npy')
	By,plot_dict = p_by.falsecolor2d('ztxy',(x,y),dyn_range=0.0)
else:
	By = None

# N.b. the plotter gives the data back with axes swapped for compatibility with pyplot.imshow routine
N = (plot_dict['ypts'].shape[0],plot_dict['xpts'].shape[0]) # t,z
tp = plot_dict['ypts']
zp = plot_dict['xpts']

print('Unboost...')
t = np.linspace(-2000,2000,Nlab[0])
z = np.linspace(-20,0,Nlab[1])
unboosted_data = Unboost((Ex,By),(tp,zp),vb,(t,z),vu)

new_file = 'unboost_'+filename+'_'+tcomp+'.npy'
ref_file = filename+'_'+tcomp+'.npy'
with open('tw_metadata.json','r') as f:
	meta = json.load(f)
	if new_file not in meta:
		meta[new_file] = {}
	meta[new_file]['axes'] = meta[ref_file]['axes']
	meta[new_file]['native units'] = 'plasma'
	meta[new_file]['grid'] = 'unboost_grid_warp.txt'
with open('tw_metadata.json','w') as f:
    json.dump(meta,f,indent=4)

with open('unboost_grid_warp.txt','w') as f:
	for i,t0 in enumerate(t):
		f.write('t = '+str(t0)+'\n')
		if i==0:
			f.write('axis1 = 0.0\n')
			f.write('axis2 = 0.0\n')
			f.write('axis3 = ')
			for x0 in z:
				if x0==z[-1]:
					f.write(str(x0)+'\n')
				else:
					f.write(str(x0)+' ')

np.save(new_file,unboosted_data[:,np.newaxis,np.newaxis,:])
