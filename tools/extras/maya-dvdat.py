import sys
import glob
import subprocess
import numpy as np
from mayavi import mlab
import twutils.plot as twplot

# This program must operate in a turboWAVE tools environment with mayavi installed as well.
# To install mayavi into a conda environment use pip:
# execute: pip install mayavi
# As of this writing do NOT use conda for installing mayavi (but you can operate in a conda environment).

# To create animation:
# Install ImageMagick suite if not already installed (e.g., sudo apt install ImageMagick)
# Give one slice index as x, this will make a sequence of image files varying that slice.

# To create subset of the animation, do the above and then from the command line (e.g.):
# convert -delay 20 test.gif'[4-6]' sub.gif
# To extend the last frame by one frame (can be repeated, effectively pauses on the last frame)
# convert -delay 20 test.gif test.gif[-1] test.gif

if len(sys.argv)<3:
	print('Usage: maya-dvdat.py [slicing=slices] [file] [type=1] [dr=0.0] [ask=yes]')
	print('slicing,slices, and file are positional arguments, the rest are key/value pairs.')
	print('slicing = 4-character string, such as xyzt.')
	print('  First axes appearing are plotted, remaining are sliced.')
	print('slices = comma delimited list of slice indices, NO SPACES.')
	print('  putting x for a slice index produces an animation sequence over that slice')
	print('  number of slices given determines dimension of plot (surface or 3D)')
	print('type = 1 or 2, 1 gives 2d-falsecolor or 3d-isosurface, 2 gives 2d-surface or 3d-volume')
	print('dr = dynamic range (0 signals linear scale)')
	print('ask = yes or no , whether to ask before overwriting temporary image files')
	exit()

# note that changing to log plot will likely require changing volume_contrast
small_pos = 1e-25
my_contours = [0.3,0.6,0.9]
volume_contrast = [0.8,0.9]

def uniform_grid_metrics(extent):
	low = np.array(extent[0::2])
	high = np.array(extent[1::2])
	origin = 0.5*(low+high)
	sizes = high-low
	return origin.astype(np.double),sizes.astype(np.double)

def uniform_grid(ax,dims,extent):
	low = np.array(extent[0::2])
	high = np.array(extent[1::2])
	ds = (high[ax]-low[ax]) / dims[ax]
	return np.linspace(low[ax]+0.5*ds,high[ax]-0.5*ds,dims[ax])

# Process command line arguments

slicing_spec = sys.argv[1].split('=')[0]
if len(slicing_spec)!=4:
	raise ValueError('slicing_spec must be 4 characters')
primitive_slices = (sys.argv[1].split('=')[1]).split(',')
num_slice_axes = len(primitive_slices)
if num_slice_axes<1 or num_slice_axes>2:
	raise ValueError('number of slice axes given should be 1 or 2')
file_to_plot = sys.argv[2]
plot_type = 1
dyn_range = 0.0
ask = 'yes'
for arg in sys.argv[3:]:
	akey = arg.split('=')[0]
	if akey not in ['type','dr','ask']:
		raise KeyError('Invalide key '+akey+' in arguments')
	aval = arg.split('=')[1]
	if akey=='type':
		plot_type = int(aval)
	if akey=='dr':
		dyn_range = np.double(aval)
	if akey=='ask':
		ask = aval

if plot_type!=1 and plot_type!=2:
	raise ValueError('Plot type must be 1 or 2')
if ask not in ['yes','no']:
	raise ValueError('ask key must have value of yes or no')

plotter = twplot.plotter(file_to_plot,buffered=True)
plotter.display_info()

# Set up animation slices

slices = []
axes = twplot.get_axis_info(slicing_spec)
for i,s in enumerate(primitive_slices):
	ax = axes[4-num_slice_axes+i]
	num = plotter.dims4()[ax]
	if s=='x':
		slices.append(np.arange(0,num))
	else:
		slices.append(np.array([np.int(s)]))

# Check existing image files and clean

img_files = glob.glob('tempfile*.png')
if len(img_files)>0 and ask=='yes':
	ans = ''
	while ans!='y' and ans!='n':
		ans = input('Found some tempfile*.png files, OK to clean (y/n) ?')
	if ans=='n':
		print('STOPPED. Please run script in a directory where there are no important files of the form tempfile*.png.')
		exit(1)

for img_file in img_files:
	subprocess.run(['rm',img_file])

# Make Plots

fig = mlab.figure(size=(600,600),bgcolor=(0,0,0))

# 2D PLOTS

if num_slice_axes==2:
	for s0 in slices[0]:
		for s1 in slices[1]:
			print('Rendering slice point',s0,',',s1,'...')
			data_slice,plot_dict = plotter.falsecolor2d(slicing_spec,(s0,s1),dyn_range,lib='mayavi')
			data_slice = data_slice.astype(np.double)
			if len(slices[0])>1 or len(slices[1])>1:
				v0,v1 = plotter.get_global_minmax(dyn_range)
			else:
				v0 = plot_dict['vmin']
				v1 = plot_dict['vmax']
			if dyn_range!=0.0:
				data_slice[np.where(data_slice<v0)] = v0
			ext = np.concatenate((plot_dict['extent'],[v0,v1]))
			origin,sizes = uniform_grid_metrics(ext)
			warp_factor = np.max(sizes[:2])/sizes[2]
			ext[4] *= warp_factor
			ext[5] *= warp_factor
			data_slice *= warp_factor
			origin,sizes = uniform_grid_metrics(ext)
			dims = plot_dict['dims']
			x1 = uniform_grid(0,dims,ext)
			x2 = uniform_grid(1,dims,ext)
			mlab.clf()
			src = mlab.pipeline.array2d_source(x1,x2,data_slice)
			# Plot a slice as a falsecolor image
			if plot_type==1:
				obj = mlab.pipeline.image_actor(src,colormap='gist_ncar')
				mlab.view(azimuth=0,elevation=0,distance=2*np.max(sizes),focalpoint=origin)
			# Plot a slice as a 3D surface
			if plot_type==2:
				warp = mlab.pipeline.warp_scalar(src,warp_scale=1.0)
				normals = mlab.pipeline.poly_data_normals(warp)
				obj = mlab.pipeline.surface(normals,colormap='gist_ncar')
				mlab.outline(extent=ext)
				#WARNING: the z-axis will include the warp factor and is highly misleading
				#mlab.axes(src,z_axis_visibility=False,extent=ext,xlabel=plot_dict['xlabel'],ylabel=plot_dict['ylabel'],zlabel=plotter.name)
				mlab.view(azimuth=90,elevation=60,distance=5*np.max(sizes),focalpoint=origin)
			s = s0*len(slices[1])+s1
			mlab.savefig('tempfile{:03d}.png'.format(s))

# 3D PLOTS

if num_slice_axes==1:
	for s in slices[0]:
		print('Rendering frame',s,'...')
		data_slice,plot_dict = plotter.volume3d(slicing_spec,s,dyn_range)
		data_slice = data_slice.astype(np.double)
		ext = plot_dict['extent']
		dims = plot_dict['dims']
		origin,sizes = uniform_grid_metrics(ext)
		x1 = np.einsum('i,j,k',uniform_grid(0,dims,ext),np.ones(dims[1]),np.ones(dims[2])).astype(np.double)
		x2 = np.einsum('i,j,k',np.ones(dims[0]),uniform_grid(1,dims,ext),np.ones(dims[2])).astype(np.double)
		x3 = np.einsum('i,j,k',np.ones(dims[0]),np.ones(dims[1]),uniform_grid(2,dims,ext)).astype(np.double)
		v0 = plot_dict['vmin']
		v1 = plot_dict['vmax']
		dv = v1-v0
		contour_list = []
		for x in my_contours:
			contour_list.append(v0+x*dv)
		mlab.clf()
		# CAUTION: the extent key has to appear in just the right places or we get confusing results
		src = mlab.pipeline.scalar_field(x1,x2,x3,data_slice,extent=ext)
		if plot_type==1:
			obj = mlab.pipeline.iso_surface(src,contours=contour_list,opacity=0.3)
		if plot_type==2:
			obj = mlab.pipeline.volume(src,vmin=v0+volume_contrast[0]*dv,vmax=v0+volume_contrast[1]*dv)
		mlab.outline(extent=ext)
		mlab.view(azimuth=-80,elevation=30,distance=3*np.max(sizes),focalpoint=origin)
		mlab.savefig('tempfile{:03d}.png'.format(s))

# Produce animation

try:
	img_files = glob.glob('tempfile*.png')
	if len(img_files)>1:
		subprocess.run(['convert','-delay','20','tempfile*.png','tempfile.gif'])
		for img_file in img_files:
			subprocess.run(['rm',img_file])
except OSError:
	print('WARNING: Failed to create animated GIF, is ImageMagick "convert" installed?')
	print('INFO: We will leave image sequence in directory.')

print('Processing complete...to quit close Maya window.')
mlab.show()
