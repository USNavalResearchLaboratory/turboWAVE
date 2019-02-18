"""Functions to process dvdat files
Contents:
class dvdatParams
def GetParameters
def PrintParameters
def GetFrameArray
def GetTimeArray"""

import numpy as np

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
		p.frames = int((p.file_size-52)/(p.floats_per_frame*4))
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

def GetFullArray(file_name,c_order=False,dbl=False):
	"""Returns numpy array containing dvdat full 4D data
	inputs:
	file_name (string)
	c_order (opt bool, false) - if true reorder array for numpy access as ary[t,x,y,z]
	dbl (opt bool, false) - if true cast to double precision after reading"""
	p = GetParameters(file_name)
	with open(file_name,"rb") as f:
		f.seek(52,0) # go to start
		the_data = np.fromfile(f,dt2,p.floats_per_frame*p.frames,"")
		the_data = the_data.reshape(p.frames,p.dims[2],p.dims[1],p.dims[0]) # dimension the resulting numpy array
	if c_order:
		the_data = the_data.swapaxes(1,3)
	if dbl:
		the_data = the_data.astype(np.double)
	return the_data

def GetFrameArray(file_name,frame_to_plot,c_order=False,dbl=False):
	"""Returns numpy array containing dvdat frame
	inputs:
	file_name (string)
	frame_to_plot (int) - can be negative per usual python rules
	c_order (opt bool, false) - if true reorder array for numpy access as ary[x,y,z]
	dbl (opt bool, false) - if true cast to double precision after reading"""
	p = GetParameters(file_name)
	with open(file_name,"rb") as f:
		if frame_to_plot<0:
			frame_to_plot = p.frames + frame_to_plot
		f.seek(52+p.floats_per_frame*4*frame_to_plot,0) # go to requested frame
		the_data = np.fromfile(f,dt2,p.floats_per_frame,"") # read the frame
		the_data = the_data.reshape(p.dims[2],p.dims[1],p.dims[0]) # dimension the resulting numpy array
	if c_order:
		the_data = the_data.swapaxes(0,2)
	if dbl:
		the_data = the_data.astype(np.double)
	return the_data

def GetTimeArray(file_name,idx):
	"""Returns numpy array containing time series from dvdat file
	inputs:
	file_name (string)
	idx (int) = index relative to each frame"""
	p = GetParameters(file_name)
	with open(file_name,"rb") as f:
		the_data = np.zeros([p.frames],dt2)
		for j in range(0,p.frames):
			f.seek(52 + j*p.floats_per_frame*4 + idx*4,0)
			the_data[j] = np.fromfile(f,dt2,1,"")[0]
	return the_data

def Create(file_name,data_array,phys_dims,c_order=False):
	"""Creates dvdat file using numpy array
	inputs:
	file_name (string)
	data_array (3D or 4D numpy array)
	phys_dims (6 element numpy array with x0,x1,y0,y1,z0,z1)
	c_order (opt bool, false) - if true assume numpy ordering : data_array[x,y,z]"""
	p = dvdatParams()
	with open(file_name,"w") as f:
		f.write(p.header)
		if c_order:
			p.dims[0] = data_array.shape[0]
			p.dims[1] = data_array.shape[1]
			p.dims[2] = data_array.shape[2]
		else:
			p.dims[0] = data_array.shape[2]
			p.dims[1] = data_array.shape[1]
			p.dims[2] = data_array.shape[0]
		p.phys_dims = np.array(phys_dims).astype(dt2)
		p.dims.tofile(f)
		p.phys_dims.tofile(f)
		if c_order:
			if len(data_array.shape)==3:
				temp = np.copy(data_array.swapaxes(0,2))
			else:
				temp = np.copy(data_array.swapaxes(1,3))
			temp.astype(dt2).tofile(f)
		else:
			data_array.astype(dt2).tofile(f)

def Append(file_name,data_array,c_order=False):
	"""Add frames to a dvdat file given a numpy array
	inputs:
	file_name (string)
	data_array (3D or 4D numpy array)
	c_order (opt bool, false) - if true assume numpy ordering : data_array[x,y,z]"""
	with open(file_name,"a") as f:
		if c_order:
			if len(data_array.shape)==3:
				temp = np.copy(data_array.swapaxes(0,2))
			else:
				temp = np.copy(data_array.swapaxes(1,3))
			temp.astype(dt2).tofile(f)
		else:
			data_array.astype(dt2).tofile(f)
