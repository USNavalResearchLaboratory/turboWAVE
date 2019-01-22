import sys
import glob
import subprocess
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import twutils.plot as twplot
import twutils.dvdat as dv
import twutils.pre as twpre
from scipy import constants as C

if len(sys.argv)<4:
	print('Usage: wigner.py [slicing=slices] [real_file] [imag_file] [opt:dynamic range (induces log plot):0] [opt:color map:viridis]')
	print('slicing is a 4-character string, such as xyzt, where the first axis is the one to transform.')
	print('slices is a comma delimited list of 3 slice indices, NO SPACES.')
	print('Optional arguments must be in the order given or else not given at all.')
	print('Dynamic range = 0 signals full range on linear scale.')
	print('Colors: viridis,magma,plasma,inferno,Spectral,bwr,seismic,prism,ocean,rainbow,jet,nipy_spectral')
	print('Color maps may be inverted by adding "_r" to the name')
	print('Note any Matplotlib color maps can be used.')
	print('Animations:')
	print('Put a python range as one of the slices to generate animated GIF.')
	print('For example, zxyt=0,0,2:5 would animate time slices 2,3,4.')
	print('Note: ImageMagick suite must be installed for animations.')
	exit()

def WignerTransform(A,ds):
	N = A.shape[0]
	M = np.int(N/2) + 1
	corr = np.zeros((N,M)).astype(np.complex)
	Ai = np.zeros(N*2-1).astype(np.complex)
	Ai[::2] = A
	Ai[1::2] = 0.5*(np.roll(Ai,1)+np.roll(Ai,-1))[1::2]
	for j in range(M):
		corr[:,j] = (np.conj(np.roll(Ai,j))*np.roll(Ai,-j))[::2]
	wig = np.fft.hfft(corr,axis=1)*ds/(2*np.pi)
	return np.fft.fftshift(wig,axes=1)

def cleanup(wildcarded_path):
	cleanstr = glob.glob(wildcarded_path)
	if len(cleanstr)>0:
		cleanstr.insert(0,'rm')
		subprocess.run(cleanstr)

def ParseSlices(dims,ax_list,slice_str_list):
	'''Function to generate a list of slice tuples for the movie.
	dims = dimensions of all 4 axes
	ax_list = slicing_spec as list of integer axis identifiers.
	slice_str_list = list of slice strings, can be indices or ranges.
	Returns slice_tuples,movie.'''
	slice_tuples = []
	range_tuples = []
	movie = False
	# Construct list of range tuples
	sax = ax_list[4-len(slice_str_list):]
	for saxi,slice_str in enumerate(slice_str_list):
		rng = slice_str.split(':')
		tup = ()
		for i,el in enumerate(rng):
			if el=='' and i==0:
				el = '0'
			if el=='' and i==1:
				el = str(dims[sax[saxi]])
			if el=='' and i==2:
				el = '1'
			tup += (int(el),)
		range_tuples.append(tup)
	# Determine the range of the movie frames
	frame_rng = range(1)
	for rng in range_tuples:
		movie = movie or len(rng)>1
		if len(rng)==2:
			frame_rng = range(rng[0],rng[1])
		if len(rng)==3:
			frame_rng = range(rng[0],rng[1],rng[2])
	# Construct list of slice tuples
	for r in frame_rng:
		tup = ()
		for rng in range_tuples:
			if len(rng)>1:
				tup += (r,)
			else:
				tup += rng
		slice_tuples.append(tup)
	return slice_tuples,movie

# normalization constants in mks

n1 = 2.8e19*1e6
su = twpre.SimUnits(n1*1e-6)
t1 = su.t1
x1 = su.x1
E1 = su.E1
U1 = C.m_e*C.c*C.c
N1 = n1*x1**3

# Matplotlib setup

mpl.rcParams.update({'text.usetex' : False , 'font.size' : 14})
my_color_map = 'viridis'
proportional = False
if proportional:
	my_aspect = 'equal'
else:
	my_aspect = 'auto'

# Process command line arguments and setup plotter object

slicing_spec = sys.argv[1].split('=')[0]
primitive_slices = (sys.argv[1].split('=')[1]).split(',')
if len(primitive_slices)!=3:
	raise ValueError('Need three slices.')
real_data_file = sys.argv[2]
imag_data_file = sys.argv[3]
if len(sys.argv)>4:
	dyn_range = np.double(sys.argv[4])
else:
	dyn_range = 0.0
if len(sys.argv)>5:
	my_color_map = sys.argv[5]
ask = 'yes'

plotter_r = twplot.plotter(real_data_file)
plotter_i = twplot.plotter(imag_data_file)
plotter_r.display_info()

# Set up animation slices

axes = twplot.get_axis_info(slicing_spec)
dims = plotter_r.dims4()
slice_tuples,movie = ParseSlices(dims,axes,primitive_slices)

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

# Determine the global color scale bounds
# If a movie we have to do all the transforms first
# We don't save the data, just do redundant transforms later

global_min = 1e50
global_max = -1e50
for file_idx,slice_now in enumerate(slice_tuples):
	abcissa,real_data,plot_dict = plotter_r.lineout(slicing_spec,slice_now,dyn_range)
	abcissa,imag_data,plot_dict = plotter_i.lineout(slicing_spec,slice_now,dyn_range)
	envelope = real_data + 1j*imag_data
	wigner = WignerTransform(envelope,abcissa[1]-abcissa[0])
	local_min = np.min(wigner)
	local_max = np.max(wigner)
	if local_min<global_min:
		global_min = local_min
	if local_max>global_max:
		global_max = local_max

# Make a movie or display single frame

for file_idx,slice_now in enumerate(slice_tuples):

	plt.figure(file_idx,figsize=(10,8),dpi=100)

	abcissa,real_data,plot_dict = plotter_r.lineout(slicing_spec,slice_now,dyn_range)
	abcissa,imag_data,plot_dict = plotter_i.lineout(slicing_spec,slice_now,dyn_range)
	envelope = real_data + 1j*imag_data
	wigner = WignerTransform(envelope,abcissa[1]-abcissa[0])
	time_extent = list(plot_dict['extent'][:2])
	dt = (time_extent[1]-time_extent[0]) / envelope.shape[0]
	freq_extent = [-np.pi/dt,np.pi/dt]
	plt.imshow(wigner.swapaxes(0,1),
		origin='lower',
		aspect=my_aspect,
		extent=time_extent+freq_extent,
		vmin=global_min,
		vmax=global_max,
		cmap=my_color_map)
	plt.colorbar()
	plt.xlabel(r'$t$',fontsize=18)
	plt.ylabel(r'$\omega$',fontsize=18)

	if movie:
		img_file = 'frame{:03d}.png'.format(file_idx)
		print('saving',img_file,'...')
		plt.savefig(img_file)
		plt.close()

if movie:
	try:
		print('Consolidating into movie file...')
		subprocess.run(["convert","-delay","30","frame*.png","mov.gif"])
		cleanup('frame*.png')
		print('Done.')
	except:
		cleanup('frame*.png')
		raise OSError("The convert program from ImageMagick may not be installed.")
else:
	plt.show()
