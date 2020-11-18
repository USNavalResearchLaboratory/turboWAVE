import sys
import glob
import os
import PIL.Image
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import twutils.plot as twplot
from scipy import constants as C
eta0 = np.sqrt(C.mu_0/C.epsilon_0)

def print_usage():
	print('Usage: wigner slicing=slices/... real_file,... imag_file,... [panels=a,b,...] [layout=2x2]')
	print('   [dr=0.0,0.0,...] [color=viridis,jet,...] [roi=h0,h1,v0,v1/...]')
	print('   [range=low,high/...] [type=<type1>,...]')
	print('------------------Examples----------------')
	print('Envelope: python wigner.py xyzt=0,0,0 e_real.npy e_imag.npy')
	print('Carrier resolved: python wigner.py xyzt=0,0,0 Ex.npy 1.0')
	print('-------------------General Notes--------------------')
	print('Extra spaces (e.g. around commas, semicolons, or equals) are not allowed.')
	print('Displays two panels, the field in real space, and the field in Wigner phase space.')
	print('Required arguments are positional, optional arguments are key=value pairs.')
	print('Double values refer to the two plots.')
	print('------------------Arguments-----------------------')
	print('slicing: 4-character string, such as xyzt, where the first axis is the one to transform.')
	print('slices: is a comma delimited list of 2 (pulse) or 3 (spectrum,wigner) slice indices.')
	print('real_file: name of the file with the real part of the field')
	print('imag_file: EITHER name of the file with the imaginary part of the field')
	print('   or the carrier frequency if the data is real.')
	print('layout: grid to layout the plots on')
	print('dr: 0.0 signals full range on linear scale, any other float is the number of decades spanned.')
	print('color: viridis,magma,plasma,inferno,Spectral,bwr,seismic,prism,ocean,rainbow,jet,nipy_spectral')
	print('   Color maps may be inverted by adding "_r" to the name')
	print('   Note any Matplotlib color maps can be used.')
	print('roi: select a subdomain to plot, otherwise the full domain is plotted.')
	print('range: range of the color bar data.')
	print('type: can be pulse,spectrum,wigner')
	print('----------------------Animations----------------------')
	print('Put a python range as one of the slices to generate animated GIF.')
	print('For example, zxyt=0,0,2:5 would animate time slices 2,3,4.')

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

def get_plotters(real_data_file,imag_data_file):
	plotter_r = twplot.plotter(real_data_file,buffered=False,units='mks')
	try:
		carrier = float(imag_data_file)
		plotter_i = 0.0
	except ValueError:
		carrier = 0.0
		plotter_i = twplot.plotter(imag_data_file,buffered=False,units='mks')
	return plotter_r,plotter_i,carrier

def form_envelope(real_field,carrier,dz,ax):
	k = 2*np.pi*np.fft.fftfreq(real_field.shape[ax],dz)
	dk = k[1]-k[0]
	kc = -k[int(real_field.shape[ax]/2)]
	carrier_idx = int(real_field.shape[ax] * carrier / (2*C.c*kc))
	env = np.fft.fft(real_field,axis=ax)
	if ax==0:
		env[int(env.shape[0]/2):,...] = 0.0
	else:
		env[...,int(env.shape[1]/2):] = 0.0
	env = np.roll(env,-carrier_idx,axis=ax)
	return 2*np.fft.ifft(env,axis=ax)

def get_spectrum(plotter_r,plotter_i,carrier,slicing_spec,slice_now,dyn_range):
	if carrier!=0.0:
		abcissa,real1d,dict1d = plotter_r.lineout(slicing_spec,slice_now,dyn_range)
		envelope1d = form_envelope(real1d,carrier,abcissa[1]-abcissa[0],0)
	else:
		abcissa,real1d,dict1d = plotter_r.lineout(slicing_spec,slice_now,dyn_range)
		abcissa,imag1d,dict1d = plotter_i.lineout(slicing_spec,slice_now,dyn_range)
		envelope1d = (real1d + 1j*imag1d)
	spectrum = np.fft.fftshift(np.abs(np.fft.fft(envelope1d))**2)
	dwpts = np.fft.fftshift(2*np.pi*C.c*np.fft.fftfreq(envelope1d.shape[0],abcissa[1]-abcissa[0]))
	return spectrum,dwpts

def get_pulse(plotter_r,plotter_i,carrier,slicing_spec,slice_now,dyn_range):
	if carrier!=0.0:
		real2d,dict2d = plotter_r.falsecolor2d(slicing_spec,slice_now[1:],dyn_range)
		abcissa,real1d,dict1d = plotter_r.lineout(slicing_spec,slice_now,dyn_range)
		envelope2d = form_envelope(real2d,carrier,abcissa[1]-abcissa[0],1)
	else:
		real2d,dict2d = plotter_r.falsecolor2d(slicing_spec,slice_now[1:],dyn_range)
		imag2d,dict2d = plotter_i.falsecolor2d(slicing_spec,slice_now[1:],dyn_range)
		envelope2d = (real2d + 1j*imag2d)
	# Put result in TV/cm
	return envelope2d*1e-12*1e-2,dict2d['extent']

def get_wigner(plotter_r,plotter_i,carrier,slicing_spec,slice_now,dyn_range):
	if carrier!=0.0:
		abcissa,real1d,dict1d = plotter_r.lineout(slicing_spec,slice_now,dyn_range)
		envelope1d = form_envelope(real1d,carrier,abcissa[1]-abcissa[0],0)
	else:
		abcissa,real1d,dict1d = plotter_r.lineout(slicing_spec,slice_now,dyn_range)
		abcissa,imag1d,dict1d = plotter_i.lineout(slicing_spec,slice_now,dyn_range)
		envelope1d = (real1d + 1j*imag1d)
	z_extent = [abcissa[0],abcissa[-1]]
	dz = (z_extent[1]-z_extent[0]) / envelope1d.shape[0]
	k_extent = [-np.pi*C.c/dz,np.pi*C.c/dz]
	wig_ext = z_extent + k_extent
	# Put result in MJ/cm2
	return WignerTransform(envelope1d,dz/C.c,eta0)*1e-6*1e-4,wig_ext

def get_plot_data(plotter_r,plotter_i,carrier,plot_type,slicing_spec,slice_now,dyn_range):
	if plot_type=='pulse':
		return get_pulse(plotter_r,plotter_i,carrier,slicing_spec,slice_now,dyn_range)
	if plot_type=='spectrum':
		return get_spectrum(plotter_r,plotter_i,carrier,slicing_spec,slice_now,dyn_range)
	if plot_type=='wigner':
		return get_wigner(plotter_r,plotter_i,carrier,slicing_spec,slice_now,dyn_range)

def main():
	# Matplotlib setup and default args

	mpl.rcParams.update({'text.usetex' : False , 'font.size' : 10})
	proportional = False
	if proportional:
		my_aspect = 'equal'
	else:
		my_aspect = 'auto'

	if len(sys.argv)<4:
		print_usage()
		exit()

	# Required arguments

	slicing_spec = []
	primitive_slices = []
	real_data_file = []
	imag_data_file = []
	for s in sys.argv[1].split('/'):
		slicing_spec += [s.split('=')[0]]
		primitive_slices += [(s.split('=')[1]).split(',')]
	for f in sys.argv[2].split(','):
		real_data_file += [f]
	for f in sys.argv[3].split(','):
		imag_data_file += [f]
	N = len(real_data_file)
	if len(imag_data_file)!=N:
		raise ValueError('Number of files does not match.')

	# Optional arguments

	color = []
	dyn_range = []
	roi = []
	val_range = []
	ask = 'yes'
	panels = []
	type = []
	cols = 1
	rows = N
	for keyval in sys.argv[4:]:
		key = keyval.split('=')[0]
		val = keyval.split('=')[1]
		if key=='panels':
			panels = val.split(',')
		if key=='layout':
			layout = val.split('x')
			rows = int(layout[0])
			cols = int(layout[1])
		if key=='dr':
			for dr in val.split(','):
				dyn_range += [float(dr)]
		if key=='color':
			color = val.split(',')
		if key=='roi':
			roi = val.split('/')
		if key=='range':
			val_range = val.split('/')
		if key=='type':
			type = val.split(',')

	# Fill in defaults
	for i in range(N-len(slicing_spec)):
		slicing_spec.append(slicing_spec[-1])
		primitive_slices.append(primitive_slices[-1])
	for i in range(N-len(panels)):
		panels.append('')
	for i in range(N-len(dyn_range)):
		dyn_range.append(0.0)
	for i in range(N-len(color)):
		color.append('viridis')
	for i in range(N-len(roi)):
		roi.append('')
	for i in range(N-len(val_range)):
		val_range.append('')
	for i in range(N-len(type)):
		type.append('wigner')

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


	# Determine the global color scale bounds for all plots
	# If a movie we have to do all the transforms first
	# We don't save the data, just do redundant transforms later

	num_frames = []
	global_min = [1e50]*N
	global_max = [-1e50]*N
	print('First pass...')
	print('-------------')
	for i in range(N):
		plotter_r,plotter_i,carrier = get_plotters(real_data_file[i],imag_data_file[i])
		axes = twplot.get_axis_info(slicing_spec[i])
		dims = plotter_r.dims4()
		slice_tuples,movie = ParseSlices(dims,axes,primitive_slices[i])
		num_frames += [len(slice_tuples)]

		if val_range[i]=='':
			for slice_now in slice_tuples:
				dat,ext = get_plot_data(plotter_r,plotter_i,carrier,type[i],slicing_spec[i],slice_now,dyn_range[i])
				local_min = np.min(dat)
				local_max = np.max(dat)
				if local_min<global_min[i]:
					global_min[i] = local_min
				if local_max>global_max[i]:
					global_max[i] = local_max
		else:
			global_min[i] = float(val_range[i].split(',')[0])
			global_max[i] = float(val_range[i].split(',')[1])

	# Make a movie or display single frame

	print('Generate images...')
	print('------------------')

	for file_idx in range(num_frames[0]):

		sz = (cols*5,rows*4)
		plt.figure(file_idx,figsize=sz,dpi=100)

		for i in range(N):
			plt.subplot(rows,cols,i+1)
			plotter_r,plotter_i,carrier = get_plotters(real_data_file[i],imag_data_file[i])
			axes = twplot.get_axis_info(slicing_spec[i])
			dims = plotter_r.dims4()
			slice_tuples,movie = ParseSlices(dims,axes,primitive_slices[i])
			slice_now = slice_tuples[file_idx]
			dat,ext = get_plot_data(plotter_r,plotter_i,carrier,type[i],slicing_spec[i],slice_now,dyn_range[i])

			if type[i]=='pulse':
				ext_scaled = list(1e15*np.array(ext[:2])/C.c) + list(1e6*np.array(ext[2:]))
				plt.imshow(np.abs(dat),
					origin='lower',
					aspect=my_aspect,
					extent=ext_scaled,
					vmin=global_min[i],
					vmax=global_max[i],
					cmap=color[i])
				b = plt.colorbar()
				b.set_label(r'${\cal E}$ (TV/cm)',size=12)
				plt.xlabel(r'$z/c - t$ (fs)',size=12)
				plt.ylabel(r'$\rho$ ($\mu$m)',size=12)

			if type[i]=='wigner':
				ext_scaled = list(1e15*np.array(ext[:2])/C.c) + list(1e-15*np.array(ext[2:]))
				plt.imshow(dat.swapaxes(0,1),
					origin='lower',
					aspect=my_aspect,
					extent=ext_scaled,
					vmin=global_min[i],
					vmax=global_max[i],
					cmap=color[i])
				b = plt.colorbar()
				b.set_label(r'${\cal N}$ (MJ/cm$^2$)',size=12)
				plt.xlabel(r'$z/c - t$ (fs)',size=12)
				plt.ylabel(r'$\delta\omega$ (fs$^{-1}$)',size=12)

			if type[i]=='spectrum':
				ext *= 1e-15
				ext_scaled = [ext[0],ext[-1],global_min[i],global_max[i]]
				plt.plot(ext,dat)
				plt.xlabel(r'$\delta\omega$ (fs$^{-1}$)',size=12)
				plt.ylabel(r'$|E(\omega)|^2$',size=12)

			if roi[i]=='':
				roi_i = ext_scaled
			else:
				roi_i = []
				for s in roi[i].split(','):
					roi_i.append(float(s))
				plt.xlim(roi_i[0],roi_i[1])
				plt.ylim(roi_i[2],roi_i[3])
			if not panels[i]=='':
				plt.text(roi_i[0],roi_i[3]+0.03*(roi_i[3]-roi_i[2]),'('+panels[i]+')')

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
