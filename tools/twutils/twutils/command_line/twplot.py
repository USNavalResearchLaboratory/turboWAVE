import sys
import glob
import os
import PIL.Image
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import twutils.plot as twplot
from scipy import constants as C

# Plotter for publication quality figures or presentation quality movies

def print_usage():
	print('Usage: twplot slicing=slices/... file1,... [panels=a,b,...] [mult=1.0,1.0,1.0/...]')
	print('    [range=low,high/...] [layout=2x2] [dr=0.0,...] [color=viridis,...] [roi=h0,h1,v0,v1/...]')
	print('    [units=mks,...] [[tex]labels=hlab,vlab,clab/...]')
	print('------------------Examples----------------')
	print('2D zx-plot: python twplot.py zxyt=0,0 rho.npy')
	print('1D x-plot: python twplot.py xyzt=0,0,0 Ex.npy')
	print('Hybrid: python twplot.py zxyt=0,0/zxyt=0,0,0 rho.npy,Ex.npy panels=a,b layout=1x2')
	print('2D zx-movie over all time slices: python twplot.py zxyt=0,: Ex.npy')
	print('-------------------General Notes--------------------')
	print('Extra spaces (e.g. around commas or equals) are not allowed.')
	print('Displays any number of panels.')
	print('Required arguments are positional, optional arguments are key=value pairs.')
	print('------------------Required Arguments-----------------------')
	print('slicing: 4-character string, such as xyzt, where the first 1 or 2 axes are plotted.')
	print('slices: is a comma delimited list of slice indices.')
	print('   Number of slices determines whether the plot is 1D or 2D.')
	print('file1: name of the first file to plot, etc. File names cannot contain commas.')
	print('------------------Optional Arguments-----------------------')
	print('panels: labels for the subplot panels, e.g., a,b,etc.')
	print('mult: rescaling multipliers, multiplying haxis,vaxis,data.')
	print('range: explicit setting of the low and high bounds of plotting range.')
	print('layout: layout of the grid on which to put multiple plots.')
	print('dr: 0.0 signals full range on linear scale, any other float is the number of decades spanned.')
	print('color: any Matplotlib color map, search internet for list (e.g., viridis, jet, inferno, etc.)')
	print('roi: select a subdomain to plot, otherwise the full domain is plotted.')
	print('units: select from plasma,atomic,natural,mks,cgs')
	print('[tex]labels: fully customizable labels for axes and colorbar.')
	print('   Including prefix tex causes each string to be automatically enclosed in $ tokens')
	print('   Some characters require special treatment if you want to use them in a label.')
	print('   Enter SPACE,COMMA,SEMICOLON,or SLASH to get the corresponding character.')
	print('   Depending on shell, various characters may need to be escaped (KNOW YOUR SHELL).')
	print('----------------------Animations----------------------')
	print('Put a python range as one of the slices to generate animated GIF.')
	print('For example, zxyt=0,0,2:5 would animate time slices 2,3,4.')

# Label scheme

def format_label(l):
	return l.replace('SPACE',' ').replace('COMMA',',').replace('SEMICOLON',';').replace('SLASH','/')

# Functions helpful for animations

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

def main():

	if len(sys.argv)<3:
		print_usage()
		exit(0)

	# Matplotlib setup

	mpl.rcParams.update({'text.usetex' : False , 'font.size' : 10})
	proportional = False
	if proportional:
		my_aspect = 'equal'
	else:
		my_aspect = 'auto'

	# Process command line arguments and setup plotter object

	slicing_spec = []
	primitive_slices = []
	for spec in sys.argv[1].split('/'):
		slicing_spec.append(spec.split('=')[0])
		primitive_slices.append(spec.split('=')[1].split(','))
	file_to_plot = sys.argv[2].split(',')
	N = len(file_to_plot)
	panels = []
	rows = 1
	cols = N
	dyn_range = []
	color = []
	roi = []
	val_range = []
	mult = []
	units = []
	labels=[]
	texlabels=[]
	ask = 'yes'
	keylist = ['panels','layout','dr','color','roi','range','mult','units','labels','texlabels']
	for keyval in sys.argv[3:]:
		key = keyval.split('=')[0]
		arg = keyval.split('=')[1]
		if key not in keylist:
			raise ValueError('Invalid key in optional arguments.')
		if key=='panels':
			panels = arg.split(',')
		if key=='layout':
			layout = arg.split('x')
			rows = int(layout[0])
			cols = int(layout[1])
		if key=='dr':
			dyn_range = arg.split(',')
		if key=='color':
			color = arg.split(',')
		if key=='roi':
			roi = arg.split('/')
		if key=='range':
			val_range = arg.split('/')
		if key=='mult':
			mult = arg.split('/')
		if key=='units':
			units = arg.split(',')
		if key=='labels':
			labels = arg.split('/')
		if key=='texlabels':
			texlabels = arg.split('/')

	for i in range(N-len(slicing_spec)):
		slicing_spec.append(slicing_spec[-1])
		primitive_slices.append(primitive_slices[-1])
	for i in range(N-len(panels)):
		panels.append('')
	for i in range(N-len(dyn_range)):
		dyn_range.append('0')
	for i in range(N-len(val_range)):
		val_range.append('')
	for i in range(N-len(color)):
		color.append('viridis')
	for i in range(N-len(roi)):
		roi.append('')
	for i in range(N-len(mult)):
		mult.append('1.0,1.0,1.0')
	for i in range(N-len(labels)):
		labels.append(',,')
	for i in range(N-len(texlabels)):
		texlabels.append(',,')
	for i in range(N-len(units)):
		units.append('plasma')

	plotter = []
	for i,f in enumerate(file_to_plot):
		needs_buffer = slicing_spec[i][0]=='t' or (slicing_spec[i][1]=='t' and len(primitive_slices[i])==2)
		plotter.append(twplot.plotter(f,buffered=needs_buffer,units=units[i]))
		plotter[-1].display_info()

	# Set up animation slices

	slice_tuples = []
	movie = []
	for i in range(N):
		axes = twplot.get_axis_info(slicing_spec[i])
		dims = plotter[i].dims4()
		st,m = ParseSlices(dims,axes,primitive_slices[i])
		slice_tuples.append(st)
		movie.append(m)

	# Check existing files and clean if there is a movie

	if movie[0]:
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

		mov_files = glob.glob('mov.gif')
		if len(mov_files)>0 and ask=='yes':
			ans = ''
			while ans!='y' and ans!='n':
				ans = input('Found mov.gif, OK to overwrite (y/n) ?')
			if ans=='n':
				print('STOPPED. Please do something with existing mov.gif and try again.')
				exit(1)

	# Determine the global color scale bounds

	global_min = 1e50*np.ones(N)
	global_max = -1e50*np.ones(N)
	for i,p in enumerate(plotter):
		if val_range[i]=='':
			for slice_now in slice_tuples[i]:
				dmult = float(mult[i].split(',')[2])
				if len(primitive_slices[i])==2:
					test_array,plot_dict = p.falsecolor2d(slicing_spec[i],slice_now,float(dyn_range[i]))
					local_min = plot_dict['vmin']*dmult
					local_max = plot_dict['vmax']*dmult
				if len(primitive_slices[i])==3:
					abcissa,test_array,plot_dict = p.lineout(slicing_spec[i],slice_now,float(dyn_range[i]))
					local_min = np.min(test_array)*dmult
					local_max = np.max(test_array)*dmult
				if local_min>local_max:
					local_min,local_max = local_max,local_min
				if local_min<global_min[i]:
					global_min[i] = local_min
				if local_max>global_max[i]:
					global_max[i] = local_max
		else:
			global_min[i] = float(val_range[i].split(',')[0])
			global_max[i] = float(val_range[i].split(',')[1])

	# Generate the movie or show a single frame.
	# The number of frames is dictated by the first panel.

	for file_idx in range(len(slice_tuples[0])):

		sz = (cols*5,rows*4)
		plt.figure(file_idx,figsize=sz,dpi=100)

		for i,p in enumerate(plotter):
			slice_now = slice_tuples[i][file_idx]
			exm = np.zeros(4)
			exm[:2] = float(mult[i].split(',')[0])
			exm[2:] = float(mult[i].split(',')[1])
			dmult = float(mult[i].split(',')[2])
			plt.subplot(rows,cols,i+1)
			if len(primitive_slices[i])==2:
				data_slice,plot_dict = p.falsecolor2d(slicing_spec[i],slice_now,float(dyn_range[i]))
				plt.imshow(data_slice*dmult,
					origin='lower',
					aspect=my_aspect,
					extent=np.array(plot_dict['extent'])*exm,
					vmin=global_min[i],
					vmax=global_max[i],
					cmap=color[i])
				b = plt.colorbar()
				l = labels[i].split(',')
				tl = texlabels[i].split(',')

				if l[0]=='' and tl[0]=='':
					plt.xlabel(plot_dict['xlabel'],size=12)
				if not l[0]=='':
					plt.xlabel(format_label(l[0]),size=12)
				if not tl[0]=='':
					plt.xlabel(format_label(r'$'+tl[0]+r'$'),size=12)

				if l[1]=='' and tl[1]=='':
					plt.ylabel(plot_dict['ylabel'],size=12)
				if not l[1]=='':
					plt.ylabel(format_label(l[1]),size=12)
				if not tl[1]=='':
					plt.ylabel(format_label(r'$'+tl[1]+r'$'),size=12)

				if l[2]=='' and tl[2]=='':
					b.set_label(plot_dict['blabel'],size=12)
				if not l[2]=='':
					b.set_label(format_label(l[2]),size=12)
				if not tl[2]=='':
					b.set_label(format_label(r'$'+tl[2]+r'$'),size=12)

				if roi[i]=='':
					roi_i = np.array(plot_dict['extent'])*exm
				else:
					roi_i = []
					for s in roi[i].split(','):
						roi_i.append(float(s))
					plt.xlim(roi_i[0],roi_i[1])
					plt.ylim(roi_i[2],roi_i[3])
				if not panels[i]=='':
					plt.text(roi_i[0],roi_i[3]+0.03*(roi_i[3]-roi_i[2]),'('+panels[i]+')')
			if len(primitive_slices[i])==3:
				abcissa,ordinate,plot_dict = p.lineout(slicing_spec[i],slice_now,float(dyn_range[i]))
				plt.plot(abcissa*exm[0],ordinate*exm[2])
				l = labels[i].split(',')
				tl = texlabels[i].split(',')

				if l[0]=='' and tl[0]=='':
					plt.xlabel(plot_dict['xlabel'],size=12)
				if not l[0]=='':
					plt.xlabel(format_label(l[0]),size=12)
				if not tl[0]=='':
					plt.xlabel(format_label(r'$'+tl[0]+r'$'),size=12)

				if l[1]=='' and tl[1]=='':
					plt.ylabel(plot_dict['ylabel'],size=12)
				if not l[1]=='':
					plt.ylabel(format_label(l[1]),size=12)
				if not tl[1]=='':
					plt.ylabel(format_label(r'$'+tl[1]+r'$'),size=12)

				if roi[i]=='':
					roi_i = [abcissa[0]*exm[0],abcissa[-1]*exm[0],global_min[i]*exm[2],global_max[i]*exm[2]]
				else:
					s = roi[i].split(',')
					roi_i = [float(s[0]),float(s[1]),float(s[2]),float(s[3])]
					plt.xlim(roi_i[0],roi_i[1])
					plt.ylim(roi_i[2],roi_i[3])
				if not panels[i]=='':
					plt.text(roi_i[0],roi_i[3]+0.03*(roi_i[3]-roi_i[2]),'('+panels[i]+')')

			plt.tight_layout()

		if movie[0]:
			img_file = 'frame{:03d}.png'.format(file_idx)
			print('saving',img_file,'...')
			plt.savefig(img_file)
			plt.close()

	if movie[0]:
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
