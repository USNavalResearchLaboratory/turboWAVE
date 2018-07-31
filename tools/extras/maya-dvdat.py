import sys
import numpy as np
from mayavi import mlab

# To create a conda environment for mayavi:
#--------------------------------------------
# conda create -n maya python=2.7 vtk=6.3.0 mayavi
# source activate maya
#-----------------------------
# Note mayavi has to use python2
# Therefore twutils may not install correctly
# As a workaround we include dvdat functions directly in this file

# To create animation:
# Give the frame to plot as -1, this will make a sequence of image files.
# Install ImageMagick suite if not already installed (sudo apt install ImageMagick)
# convert -delay 20 test*.png test.gif
# To create subset of the animation
# convert -delay 20 test.gif'[4-6]' sub.gif
# To extend the last frame by one frame (can be repeated, effectively pauses on the last frame)
# convert -delay 20 test.gif src.gif[-1] test.gif

if len(sys.argv)!=4:
	print('Usage: maya-dvdat.py filename frame type')
	print('frame = -1 induces file sequence for all frames')
	print('type = 1,2,3,4')
	print('1 = falsecolor image of a slice')
	print('2 = surface plot of a slice')
	print('3 = 3D isocontours')
	print('4 = 3D volumetric rendering')
	exit()

filename=sys.argv[1]
frame_to_plot = np.int(sys.argv[2])
plot_type = np.int(sys.argv[3])
# axis numbers go as t,x,y,z
# only matters for 2D plots
haxis = 1
vaxis = 2
faxis = 0
saxis = 3
# note that changing to log plot will likely require changing volume_contrast
log_plot = False
dyn_range = 5.0
my_contours = [0.02,0.05,0.08,0.12]
volume_contrast = [0.02,0.4]
hlbl = 'Axis'+str(haxis)
vlbl = 'Axis'+str(vaxis)
datalbl = filename.split('.')[0]

# define the data types occuring in a dvdat file
dt1 = np.dtype(np.int32)
dt1 = dt1.newbyteorder('>')
dt2 = np.dtype(np.float32)
dt2 = dt2.newbyteorder('>')

class dvdatParams:
	"""Contains parameters describing a dvdat file"""
	header = 'DataViewer 2.0.0'
	file_size = 0
	frames = 0
	floats_per_frame = 0
	dims = np.array([1,1,1],dtype=dt1)
	phys_dims = np.array([0.,1.,0.,1.,0.,1.],dtype=dt2)

def GetParameters(file_name):
	"""Returns the parameters of the dvdat file
	inputs:
	file_name (string)"""
	p = dvdatParams()
	with open(file_name,"rb") as f:
		f.seek(16,0)
		p.dims = np.fromfile(f,dt1,3)
		p.phys_dims = np.fromfile(f,dt2,6)
		f.seek(0,0)
		# for some reason, in python3, the following cannot come before np.fromfile(...)
		p.header = f.read(16) # read header string "DataViewer 2.0.0"
		f.seek(0,2) # go to end of file
		p.file_size = f.tell() # get length of file
		p.floats_per_frame = p.dims[0]*p.dims[1]*p.dims[2]
		p.frames = (p.file_size-52)/(p.floats_per_frame*4)
	return p

def PrintParameters(p):
	"""Prints the parameters of the dvdat file
	inputs:
	p (dvdatParams)"""
	print("------------------------------------------------------------------")
	print("File properties:")
	print("------------------------------------------------------------------")
	print("header string: ",p.header)
	print("length: ",p.file_size)
	print("frames: ",p.frames)
	print("index dimensions: ",p.dims)
	print('X Bounds: {0:.3f} , {1:.3f}'.format(float(p.phys_dims[0]),float(p.phys_dims[1])))
	print('Y Bounds: {0:.3f} , {1:.3f}'.format(float(p.phys_dims[2]),float(p.phys_dims[3])))
	print('Z Bounds: {0:.3f} , {1:.3f}'.format(float(p.phys_dims[4]),float(p.phys_dims[5])))
	print("------------------------------------------------------------------")

def GetFrameArray(file_name,frame_to_plot,c_order=False,dbl=False):
	"""Returns numpy array containing dvdat frame
	inputs:
	file_name (string)
	frame_to_plot (int)"""
	p = GetParameters(file_name)
	with open(file_name,"rb") as f:
		f.seek(52+p.floats_per_frame*4*frame_to_plot,0) # go to requested frame
		the_data = np.fromfile(f,dt2,p.floats_per_frame,"") # read the frame
		the_data = the_data.reshape(p.dims[2],p.dims[1],p.dims[0]) # dimension the resulting numpy array
	if c_order:
		the_data = the_data.swapaxes(0,2)
	if dbl:
		the_data = the_data.astype(np.double)
	return the_data

params = GetParameters(filename)
PrintParameters(params)

if frame_to_plot==-1:
	data = np.zeros((params.frames,params.dims[0],params.dims[1],params.dims[2]))
	for frame_now in range(params.frames):
		data[frame_now,...] = GetFrameArray(filename,frame_now,c_order=True,dbl=True)
else:
	data = np.zeros((1,params.dims[0],params.dims[1],params.dims[2]))
	data[0,...] = GetFrameArray(filename,frame_to_plot,c_order=True,dbl=True)

if log_plot==True:
	data = np.log10(np.abs(data))
	cutoff = np.max(data)-dyn_range
	data[np.where(data<cutoff)]=cutoff

# SET UP SOME VARIABLES
v1 = np.max(data)
v0 = np.min(data)
dv = v1-v0
slice_to_plot = np.int(data.shape[saxis]/2)
contour_list = []
for x in my_contours:
	contour_list.append(v0+x*dv)
data_slice = np.take(data,[slice_to_plot],axis=saxis)
if haxis>vaxis:
	data_slice = data_slice.swapaxes(haxis,vaxis)
pd = np.concatenate(([0.0,1.0],params.phys_dims))
sizes = np.array([pd[3]-pd[2],pd[5]-pd[4],pd[7]-pd[6]])
extent_2d = [pd[haxis*2],pd[haxis*2+1],pd[vaxis*2],pd[vaxis*2+1],np.min(data_slice),np.max(data_slice)]
extent_2d = np.array(extent_2d)
warp_factor = (extent_2d[1]-extent_2d[0])/np.max(np.abs(data_slice))
extent_2d[4] *= warp_factor
extent_2d[5] *= warp_factor
center_2d = (extent_2d[0]/2+extent_2d[1]/2,extent_2d[2]/2+extent_2d[3]/2,extent_2d[4]/2+extent_2d[5]/2)
center_3d = (pd[2]/2+pd[3]/2,pd[4]/2+pd[5]/2,pd[6]/2+pd[7]/2)

# CREATE THE FIGURES

fig = mlab.figure(size=(600,600),bgcolor=(0,0,0))

# Plot a slice as a falsecolor image
if plot_type==1:
	frame_slice = np.take(data_slice,[0],axis=faxis)
	frame_slice = np.squeeze(frame_slice,axis=(faxis,saxis))
	src = mlab.pipeline.array2d_source(frame_slice)
	obj = mlab.pipeline.image_actor(src,colormap='gist_ncar',extent=extent_2d)
	mlab.axes(src,z_axis_visibility=False,extent=extent_2d,xlabel=hlbl,ylabel=vlbl,zlabel=datalbl)

# Plot a slice as a 3D surface
if plot_type==2:
	frame_slice = np.take(data_slice,[0],axis=faxis)
	frame_slice = np.squeeze(frame_slice,axis=(faxis,saxis))
	src = mlab.pipeline.array2d_source(frame_slice * warp_factor)
	warp = mlab.pipeline.warp_scalar(src,warp_scale=1.0)
	normals = mlab.pipeline.poly_data_normals(warp)
	obj = mlab.pipeline.surface(normals,colormap='gist_ncar',extent=extent_2d)
	mlab.outline()
	mlab.axes(src,extent=extent_2d,xlabel=hlbl,ylabel=vlbl,zlabel=datalbl)

# Plot 3D data as isocontours
if plot_type==3:
	obj = mlab.contour3d(data[0,...],contours=contour_list,opacity = 0.3,extent=params.phys_dims)

# Full volumetric rendering
# May have to tweak vmin and vmax to see inside effectively
# May have to tweak color scale as well (how?)
if plot_type==4:
	src = mlab.pipeline.scalar_field(data[0,...])
	obj = mlab.pipeline.volume(src,vmin=v0+volume_contrast[0]*dv,vmax=v0+volume_contrast[1]*dv)

if plot_type==1 or plot_type==2:
	mlab.view(azimuth=-80,elevation=45,distance=3*np.max(sizes),focalpoint=center_2d)
else:
	mlab.view(azimuth=-80,elevation=30)#,distance=3*np.max(sizes),focalpoint=center_3d)

if frame_to_plot==-1:
	for frame_now in range(data.shape[faxis]):
		print('Rendering frame',frame_now,'...')
		if plot_type==1 or plot_type==2:
			frame_slice = np.take(data_slice,[frame_now],axis=faxis)
			frame_slice = np.squeeze(frame_slice,axis=(faxis,saxis))
			obj.mlab_source.scalars = frame_slice * warp_factor
		else:
			obj.mlab_source.scalars = data[frame_now,...]
		#mlab.text3d(120,50,0.0,'e-',scale=20.0)
		#mlab.text3d(30,-180,-70,'e+',scale=20.0)
		mlab.savefig('test{:03d}.png'.format(frame_now))

print('Processing complete...to quit close Maya window.')
mlab.show()
