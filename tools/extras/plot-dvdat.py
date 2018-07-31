#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import twutils.plot as twplot
import twutils.dvdat as dv
import twutils.pre as twpre
from scipy import constants as C

if len(sys.argv)<3:
	print('Usage: plot-dvdat.py [slicing=slices] [file] [opt:dynamic range (induces log plot):0] [opt:color map:viridis]')
	print('slicing is a 4-character string, such as xyzt, where the first axes appearing are the ones that are plotted.')
	print('slices is a comma delimited list of slice indices, NO SPACES.')
	print('If 2 slice indices are given a 2D plot is produced, if 3 a 1D plot is produced.')
	print('Optional arguments must be in the order given or else not given at all.')
	print('Dynamic range = 0 signals full range on linear scale.')
	print('Colors: viridis,magma,plasma,inferno,Spectral,bwr,seismic,prism,ocean,rainbow,jet,nipy_spectral')
	print('Color maps may be inverted by adding "_r" to the name')
	print('Note any Matplotlib color maps can be used.')
	exit()

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

# Process command line arguments

slicing_spec = sys.argv[1].split('=')[0]
slices = (sys.argv[1].split('=')[1]).split(',')
file_to_plot = sys.argv[2]
if len(sys.argv)>3:
	dyn_range = np.double(sys.argv[3])
else:
	dyn_range = 0.0
if len(sys.argv)>4:
	my_color_map = sys.argv[4]

plotter = twplot.plotter(file_to_plot)
plotter.display_info()

plt.figure(1,figsize=(10,8),dpi=100)

if len(slices)==2:
	data_slice,plot_dict = plotter.falsecolor2d(slicing_spec,slices,dyn_range)
	plt.imshow(data_slice,
		origin='lower',
		aspect=my_aspect,
		extent=plot_dict['extent'],
		vmin=plot_dict['vmin'],
		vmax=plot_dict['vmax'],
		cmap=my_color_map)
	plt.colorbar()
	plt.xlabel(plot_dict['xlabel'],fontsize=18)
	plt.ylabel(plot_dict['ylabel'],fontsize=18)

if len(slices)==3:
	abcissa,ordinate,plot_dict = plotter.lineout(slicing_spec,slices,dyn_range)
	plt.plot(abcissa,ordinate)
	plt.xlabel(plot_dict['xlabel'],fontsize=18)
	plt.ylabel(plot_dict['ylabel'],fontsize=18)

plt.show()
