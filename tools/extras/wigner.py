import sys
import glob
import os
import PIL.Image
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import twutils.plot as twplot
import twutils.pre as twpre
from scipy import constants as C

if len(sys.argv)<3:
	print('Usage: python wigner.py slicing=slices real_file,imag_file [panels=a,b] [layout=1x2]')
	print('   [dr=0.0,0.0] [color=viridis,jet] [roi=h0,h1,v0,v1/h0,h1,v0,v1]')
	print('------------------Examples----------------')
	print('Envelope: python wigner.py xyzt=0,0,0 e_real.dvdat,e_imag.dvdat')
	print('Carrier resolved: python wigner.py xyzt=0,0,0 Ex.dvdat,1.0')
	print('-------------------General Notes--------------------')
	print('Extra spaces (e.g. around commas, semicolons, or equals) are not allowed.')
	print('Displays two panels, the field in real space, and the field in Wigner phase space.')
	print('Required arguments are positional, optional arguments are key=value pairs.')
	print('Double values refer to the two plots.')
	print('------------------Arguments-----------------------')
	print('slicing: 4-character string, such as xyzt, where the first axis is the one to transform.')
	print('slices: is a comma delimited list of 3 slice indices.')
	print('   The real space plot only uses the last 2 slices.')
	print('real_file: name of the file with the real part of the field')
	print('imag_file: EITHER name of the file with the imaginary part of the field')
	print('   or the carrier frequency if the data is real.')
	print('layout: can be 1x2 or 2x1')
	print('dr: 0.0 signals full range on linear scale, any other float is the number of decades spanned.')
	print('color: viridis,magma,plasma,inferno,Spectral,bwr,seismic,prism,ocean,rainbow,jet,nipy_spectral')
	print('   Color maps may be inverted by adding "_r" to the name')
	print('   Note any Matplotlib color maps can be used.')
	print('roi: select a subdomain to plot, otherwise the full domain is plotted.')
	print('----------------------Animations----------------------')
	print('Put a python range as one of the slices to generate animated GIF.')
	print('For example, zxyt=0,0,2:5 would animate time slices 2,3,4.')
	print('Note: ImageMagick suite must be installed for animations.')
	print('On Windows hard coding the path to the correct convert.exe may be required.')
	exit()

def WignerTransform(A,ds,eta0):
	# If A is in V/m, ds in secs, and eta0 in ohms, returns J/m^2.
	# Then, integration over consistent phase space units (dimensionless product) gives the actual fluence.
	# The frequency variable is assumed to be an angular frequency (e.g. rad/s)
	N = A.shape[0]
	M = np.int(N/2) + 1
	corr = np.zeros((N,M)).astype(np.complex)
	Ai = np.zeros(N*2-1).astype(np.complex)
	Ai[::2] = A
	Ai[1::2] = 0.5*(np.roll(Ai,1)+np.roll(Ai,-1))[1::2]
	for j in range(M):
		corr[:,j] = (np.conj(np.roll(Ai,j))*np.roll(Ai,-j))[::2]
	wig = np.fft.hfft(corr,axis=1)*ds/(2*np.pi)
	wig = np.fft.fftshift(wig,axes=1)
	# Fix the units
	dw = (2*np.pi/ds) / N
	fluence = 0.5*np.sum(np.abs(A)**2)*ds/eta0
	wigfluence = np.sum(wig)*ds*dw
	return wig*fluence/wigfluence

def cleanup(wildcarded_path):
	cleanstr = glob.glob(wildcarded_path)
	for f in cleanstr:
		os.remove(f)

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

n1 = 2.65e17*1e6
su = twpre.SimUnits(n1*1e-6)
t1 = su.t1
x1 = su.x1
E1 = su.E1
U1 = C.m_e*C.c*C.c
N1 = n1*x1**3
eta0 = np.sqrt(C.mu_0/C.epsilon_0)

# Matplotlib setup and default args

mpl.rcParams.update({'text.usetex' : False , 'font.size' : 10})
color = ['viridis','viridis']
proportional = False
if proportional:
	my_aspect = 'equal'
else:
	my_aspect = 'auto'
dyn_range = [0.0,0.0]
roi = [[],[]]
ask = 'yes'
layout = '1x2'
panels = ''

# Process command line arguments and setup plotter object

slicing_spec = sys.argv[1].split('=')[0]
primitive_slices = (sys.argv[1].split('=')[1]).split(',')
if len(primitive_slices)!=3:
	raise ValueError('Need three slices.')
real_data_file = sys.argv[2].split(',')[0]
imag_data_file = sys.argv[2].split(',')[1]
for keyval in sys.argv[3:]:
	key = keyval.split('=')[0]
	val = keyval.split('=')[1]
	if key=='panels':
		panels = val.split(',')
	if key=='layout':
		layout = val
	if key=='dr':
		dyn_range = []
		dyn_range.append(float(val.split(',')[0]))
		dyn_range.append(float(val.split(',')[1]))
	if key=='color':
		color = val.split(',')
	if key=='roi':
		for s in val.split('/')[0].split(','):
			roi[0].append(int(s))
		for s in val.split('/')[1].split(','):
			roi[1].append(int(s))

plotter_r = twplot.plotter(real_data_file,buffered=False)
try:
	carrier = float(imag_data_file)
	plotter_i = 0.0
except ValueError:
	carrier = 0.0
	plotter_i = twplot.plotter(imag_data_file,buffered=False)
plotter_r.display_info()

# Set up animation slices

axes = twplot.get_axis_info(slicing_spec)
dims = plotter_r.dims4()
slice_tuples,movie = ParseSlices(dims,axes,primitive_slices)

# Check existing image files and clean

img_files = glob.glob('frame*.png')
if len(img_files)>0 and ask=='yes':
	ans = ''
	while ans!='y' and ans!='n':
		ans = input('Found some frame*.png files, OK to clean (y/n) ?')
	if ans=='n':
		print('STOPPED. Please run script in a directory where there are no important files of the form frame*.png.')
		exit(1)

for img_file in img_files:
	os.remove(img_file)

def form_envelope(real_field,carrier,dz,ax):
	k = 2*np.pi*np.fft.fftfreq(real_field.shape[ax],dz)
	dk = k[1]-k[0]
	kc = -k[int(real_field.shape[ax]/2)]
	carrier_idx = int(real_field.shape[ax] * carrier / (2*kc))
	env = np.fft.fft(real_field,axis=ax)
	if ax==0:
		env[int(env.shape[0]/2):,...] = 0.0
	else:
		env[...,int(env.shape[1]/2):] = 0.0
	env = np.roll(env,-carrier_idx,axis=ax)
	return 2*E1*np.fft.ifft(env,axis=ax)

def extract_plot_data(plotter_r,plotter_i,slice_now):
	if carrier!=0.0:
		real2d,dict2d = plotter_r.falsecolor2d(slicing_spec,slice_now[1:],dyn_range[0])
		abcissa,real1d,dict1d = plotter_r.lineout(slicing_spec,slice_now,dyn_range[1])
		envelope2d = form_envelope(real2d,carrier,abcissa[1]-abcissa[0],1)
		envelope1d = form_envelope(real1d,carrier,abcissa[1]-abcissa[0],0)
	else:
		real2d,dict2d = plotter_r.falsecolor2d(slicing_spec,slice_now[1:],dyn_range[0])
		imag2d,dict2d = plotter_i.falsecolor2d(slicing_spec,slice_now[1:],dyn_range[0])
		abcissa,real1d,dict1d = plotter_r.lineout(slicing_spec,slice_now,dyn_range[1])
		abcissa,imag1d,dict1d = plotter_i.lineout(slicing_spec,slice_now,dyn_range[1])
		envelope2d = E1*(real2d + 1j*imag2d)
		envelope1d = E1*(real1d + 1j*imag1d)
	z_extent = list(dict1d['extent'][:2])
	dz = (z_extent[1]-z_extent[0]) / envelope1d.shape[0]
	k_extent = [-np.pi/dz,np.pi/dz]
	wig_ext = z_extent + k_extent
	return envelope2d,dict2d['extent'],WignerTransform(envelope1d,dz*x1/C.c,eta0),wig_ext


# Determine the global color scale bounds for both plots
# If a movie we have to do all the transforms first
# We don't save the data, just do redundant transforms later

global_min1 = 1e50
global_max1 = -1e50
global_min2 = 1e50
global_max2 = -1e50
for slice_now in slice_tuples:
	envelope2d,ext1,wigner,ext2=extract_plot_data(plotter_r,plotter_i,slice_now)
	local_min = np.min(np.abs(envelope2d))
	local_max = np.max(np.abs(envelope2d))
	if local_min<global_min1:
		global_min1 = local_min
	if local_max>global_max1:
		global_max1 = local_max
	local_min = np.min(wigner)
	local_max = np.max(wigner)
	if local_min<global_min2:
		global_min2 = local_min
	if local_max>global_max2:
		global_max2 = local_max

# Make a movie or display single frame

for file_idx,slice_now in enumerate(slice_tuples):

	if layout=='1x2':
		plt.figure(file_idx,figsize=(8,3.5),dpi=150)
	else:
		plt.figure(file_idx,figsize=(3.5,7),dpi=150)

	envelope2d,ext1,wigner,ext2=extract_plot_data(plotter_r,plotter_i,slice_now)
	if layout=='1x2':
		plt.subplot(121)
	else:
		plt.subplot(211)
	plt.imshow(np.abs(envelope2d)*1e-12*1e-2,
		origin='lower',
		aspect=my_aspect,
		extent=ext1,
		vmin=global_min1*1e-12*1e-2,
		vmax=global_max1*1e-12*1e-2,
		cmap=color[0])
	b = plt.colorbar()
	b.set_label(r'${\cal E}$ (TV/cm)',size=12)
	if len(roi[0])==4:
		plt.xlim(roi[0][0],roi[0][1])
		plt.ylim(roi[0][2],roi[0][3])
	else:
		roi[0] = ext1
	plt.xlabel(r'$\omega_p(z/c - t)$',size=12)
	plt.ylabel(r'$\omega_p\rho/c$',size=12)
	if not panels=='':
		plt.text(roi[0][0],roi[0][3]+0.03*(roi[0][3]-roi[0][2]),'('+panels[0]+')')

	if layout=='1x2':
		plt.subplot(122)
	else:
		plt.subplot(212)
	plt.imshow(wigner.swapaxes(0,1)*1e-6*1e-4,
		origin='lower',
		aspect=my_aspect,
		extent=ext2,
		vmin=global_min2*1e-6*1e-4,
		vmax=global_max2*1e-6*1e-4,
		cmap=color[1])
	b = plt.colorbar()
	b.set_label(r'${\cal N}$ (MJ/cm$^2$)',size=12)
	if len(roi[1])==4:
		plt.xlim(roi[1][0],roi[1][1])
		plt.ylim(roi[1][2],roi[1][3])
	else:
		roi[1] = ext2
	plt.xlabel(r'$\omega_p(z/c - t)$',size=12)
	plt.ylabel(r'$\delta\omega/\omega_p$',size=12)
	if not panels=='':
		plt.text(roi[1][0],roi[1][3]+0.03*(roi[1][3]-roi[1][2]),'('+panels[1]+')')

	plt.tight_layout()

	if movie:
		img_file = 'frame{:03d}.png'.format(file_idx)
		print('saving',img_file,'...')
		plt.savefig(img_file)
		plt.close()

if movie:
	print('Consolidating into movie file...')
	images = []
	frameRateHz = 5
	img_files = sorted(glob.glob('frame*.png'))
	for f in img_files:
		images.append(PIL.Image.open(f))
	images[0].save('mov.gif',save_all=True,append_images=images[1:],duration=int(1000/frameRateHz),loop=0)
	cleanup('frame*.png')
	print('Done.')
else:
	plt.show()
