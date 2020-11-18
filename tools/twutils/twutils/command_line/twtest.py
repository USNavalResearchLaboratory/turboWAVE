import os
import glob
import ntpath
import sys
import traceback
import subprocess
import datetime
import numpy as np
import re
import twutils.plot as twplot
import matplotlib as mpl
import matplotlib.pyplot as plt
import PIL.Image

def print_usage():
	print('Usage: twtest.py tw_root cmd')
	print('tw_root = path to turboWAVE root directory.')
	print('cmd = command line to execute for each test.')
	print('To select explicit categories, append then to tw_root with double colons.')
	print('desktop example: python twtest.py ~/turboWAVE::hydro::pic tw3d -n 4 -c 5')
	print('cluster example: python twtest.py ~/turboWAVE mpirun -np 4 tw3d')
	print('N.b. as a corollary no double colons may appear in tw_root.')

def cleanup(wildcarded_path):
	cleanstr = glob.glob(wildcarded_path)
	for f in cleanstr:
		os.remove(f)

def gen_plot(img_file,plotter,slicing_spec,slices,dyn_range=0.0,val_range=(0.0,0.0),my_color_map='viridis',global_contrast=False):
	fig = plt.figure(figsize=(7,5),dpi=75)
	if len(slices)==2:
		data_slice,plot_dict = plotter.falsecolor2d(slicing_spec,slices,dyn_range)
		if global_contrast:
			minval,maxval = plotter.get_global_minmax(dyn_range)
		else:
			minval = plot_dict['vmin']
			maxval = plot_dict['vmax']
		if val_range!=(0.0,0.0):
			minval = val_range[0]
			maxval = val_range[1]
		plt.imshow(data_slice,
			origin='lower',
			aspect='auto',
			extent=plot_dict['extent'],
			vmin=minval,
			vmax=maxval,
			cmap=my_color_map)
		bar = plt.colorbar()
		plt.xlabel(plot_dict['xlabel'],fontsize=18)
		plt.ylabel(plot_dict['ylabel'],fontsize=18)
		bar.set_label(plot_dict['blabel'],fontsize=18)
	if len(slices)==3:
		abcissa,ordinate,plot_dict = plotter.lineout(slicing_spec,slices,dyn_range)
		plt.plot(abcissa,ordinate)
		plt.xlabel(plot_dict['xlabel'],fontsize=18)
		plt.ylabel(plot_dict['ylabel'],fontsize=18)
	plt.tight_layout()
	plt.savefig(img_file)
	plt.close()

def gen_movie(mov_file,plotter,slicing_spec,slices,dyn_range=0.0,val_range=(0.0,0.0),my_color_map='viridis'):
	frames = plotter.dims4()[0]
	for i in range(frames):
		img_file = 'frame{:03d}.png'.format(i)
		c = slicing_spec.find('t')
		slices[c-2] = str(i)
		gen_plot(img_file,plotter,slicing_spec,tuple(map(int,slices)),dyn_range,val_range,my_color_map,True)
	images = []
	frameRateHz = 5
	img_files = sorted(glob.glob('frame*.png'))
	for f in img_files:
		images.append(PIL.Image.open(f))
	images[0].save(mov_file,save_all=True,append_images=images[1:],duration=int(1000/frameRateHz),loop=0)
	cleanup('frame*.png')

def gen_viz(primitive_file,fig_dict):
	'''Figure out whether a still image or movie was requested and call appropriate routine'''
	file_to_plot = fig_dict['data']
	slicing_spec = fig_dict['slicing_spec']
	slices = fig_dict['slices']
	dyn_range = fig_dict['dyn_range']
	val_range = fig_dict['val_range']
	my_color_map = fig_dict['color']
	plotter = twplot.plotter(file_to_plot,units=fig_dict['units'])
	c = slicing_spec.find('t')
	if (c==2 or c==3) and slices[c-2]==':':
		gen_movie(primitive_file+'.gif',plotter,slicing_spec,slices,dyn_range,val_range,my_color_map)
		return primitive_file+'.gif'
	else:
		gen_plot(primitive_file+'.png',plotter,slicing_spec,tuple(map(int,slices)),dyn_range,val_range,my_color_map)
		return primitive_file+'.png'

def comment_remover(text):
	'''Function to strip comments from TW input file
	From stackoverflow.com, author Markus Jarderot'''
	def replacer(match):
		s = match.group(0)
		if s.startswith('/'):
			return " " # note: a space and not an empty string
		else:
			return s
	pattern = re.compile(
		r'//.*?$|/\*.*?\*/|\'(?:\\.|[^\\\'])*\'|"(?:\\.|[^\\"])*"',
		re.DOTALL | re.MULTILINE
	)
	return re.sub(pattern, replacer, text)

def parse_input_file(ex_path):
	'''Put important info from input file into a dictionary'''
	fig_dict = { }
	input_dict = { }
	with open(ex_path) as f:
		for line in f:
			c = line.find('TWTEST')
			if c>=0:
				args = line[c:].split()
				if args[1]!='matplotlib':
					raise ValueError('Example file with invalid plotter'+args[1])
				fig_dict['slicing_spec'] = args[2].split('=')[0]
				slices = []
				for s in args[2].split('=')[1].split(','):
					slices.append(s)
				fig_dict['slices'] = slices
				fig_dict['data'] = args[3]
				fig_dict['dyn_range'] = 0.0
				fig_dict['val_range'] = (0.0,0.0)
				fig_dict['color'] = 'viridis'
				fig_dict['units'] = 'plasma'
				if len(args)>4:
					for keyval in args[4:]:
						key = keyval.split('=')[0]
						val = keyval.split('=')[1]
						if key=='dr':
							fig_dict['dyn_range'] = np.float(val)
						if key=='range':
							r = val.split(',')
							fig_dict['val_range'] = ( np.float(r[0]) , np.float(r[1]) )
						if key=='color':
							fig_dict['color'] = val
						if key=='units':
							fig_dict['units'] = val
	if len(fig_dict)==0:
		return fig_dict,input_dict
	with open(ex_path) as f:
		data = f.read()
		data = comment_remover(data)
		words = data.replace("="," ").replace(","," ").replace("("," ").replace(")"," ").split()
		start_looking = False
		user_vars = {}
		for i,word in enumerate(words):
			if word=="#define":
				user_vars[words[i+1]] = words[i+2]
			if word=='new' and words[i+1]=='grid':
				start_looking = True
			if start_looking and word=='dimensions':
				subs = []
				for d in range(1,4):
					if words[i+d] in user_vars:
						subs += [user_vars[words[i+d]]]
					else:
						subs += [words[i+d]]
				input_dict['dims'] = (int(subs[0]),int(subs[1]),int(subs[2]))
				start_looking = False
	return fig_dict,input_dict

def parse_cmd(arg_list):
	cmd = ''
	mpi_procs = 1
	omp_threads = 1
	for i,arg in enumerate(arg_list):
		if arg=='-n' or arg=='-np':
			mpi_procs = int(arg_list[i+1])
			arg_list[i+1] = 'TWPROCS'
		if arg=='-c':
			omp_threads = int(arg_list[i+1])
			arg_list[i+1] = 'TWTHREADS'
		cmd = cmd + arg + ' '
	return cmd,mpi_procs,omp_threads

def optimize_parallel(num_procs,num_threads,dims):
	num_dims = 0
	dim1 = 0
	if dims[0]>1:
		dim1 = dims[0]
		num_dims += 1
	if dims[1]>1:
		dim1 = dims[1]
		num_dims += 1
	if dims[2]>1:
		dim1 = dims[2]
		num_dims += 1
	if num_dims==1:
		req_nodes = num_procs*num_threads
		mpi_nodes = req_nodes
		while dim1%(2*mpi_nodes)!=0:
			mpi_nodes -= 1
		return mpi_nodes,1
	else:
		return num_procs,num_threads

def main():

	if len(sys.argv)<2:
		print_usage()
		exit()

	subargs = sys.argv[1].split('::')
	tw_root = subargs[0]
	if not os.path.isdir(tw_root):
		print('You specified a root directory ('+tw_root+') which does not exist.')
		exit(1)
	if len(subargs)>1:
		req_categories = subargs[1:]
	else:
		req_categories = []

	mpl.rcParams.update({'text.usetex' : False , 'font.size' : 16})

	html_doc = ''
	html_doc += '<!DOCTYPE html>'
	html_doc += '<html><head><title>TurboWAVE Report</title></head>'
	html_doc += '<body>'

	html_doc += '<h1 style="background-color:black;color:white;">TurboWAVE Test Suite Report</h1>\n\n'

	html_doc += '<h2 style="background-color:rgb(0,20,100);color:white;">Version Information</h2>\n\n'
	html_doc += '<p>Date of report : ' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	# Get the command but insert placeholders for procs and threads ('TWPROCS' and 'TWTHREADS')
	# The requested procs and threads are stored separately
	tw_cmd,mpi_procs,omp_threads = parse_cmd(sys.argv[2:])
	html_doc += '<p>Requested MPI processes = ' + str(mpi_procs) + '</p>'
	html_doc += '<p>Requested OMP threads = ' + str(omp_threads) + '</p>'

	cleanup('*.gif')
	cleanup('*.png')
	cleanup('*.npy')

	save_cwd = os.getcwd()
	os.chdir(tw_root)
	try:
		compl = subprocess.run(["git","status"],stdout=subprocess.PIPE,universal_newlines=True)
		html_doc += '<p>Git status:</p>'
		html_doc += '<blockquote><pre>'+compl.stdout+'</pre></blockquote>'

		compl = subprocess.run(["git","log"],stdout=subprocess.PIPE,universal_newlines=True)
		html_doc += '<p>Git last commit:</p>'
		gitloglines = compl.stdout.splitlines(keepends=True)
		html_doc += '<blockquote><pre>'
		commit_count = 0
		for l in gitloglines:
			if l[:6]=='commit':
				commit_count += 1
			if commit_count==2:
					break
			html_doc += l
		html_doc += '</pre></blockquote>'
	except:
		html_doc += '<p>ERROR: something went wrong in trying to call git for versioning information.</p>'
	os.chdir(save_cwd)

	try:
		category_path_list = glob.glob(tw_root+'/core/examples/*')
		category_path_list = sorted(category_path_list,key=lambda s: ntpath.basename(s)[0])
		category_list = []
		for s in category_path_list:
			category_list.append(ntpath.basename(s))
		if len(req_categories)>0:
			for req in req_categories:
				if req not in category_list:
					raise ValueError('Requested category '+req+' not found.')
			category_list = req_categories

		for cat in category_list:
			print('=====================================')
			print('Category',cat)
			html_doc += '\n\n<h2 style="background-color:rgb(0,20,100);color:white;">"' + cat + '" Subdirectory</h2>\n\n'
			ex_path_list = glob.glob(tw_root+'/core/examples/'+cat+'/*')
			ex_path_list = sorted(ex_path_list,key=lambda s: ntpath.basename(s)[0])
			ex_list = []
			for s in ex_path_list:
				ex_list.append(ntpath.basename(s))
			for i,ex_path in enumerate(ex_path_list):
				print('--------------------------------------')
				print('Example',ex_list[i])
				html_doc += '\n<h3 style="background-color:rgb(100,250,200);">"' + ex_list[i] + '" Example</h3>\n'

				fig_dict,input_dict = parse_input_file(ex_path)

				if len(fig_dict)>0:
					# Run this case
					nprocs,nthreads = optimize_parallel(mpi_procs,omp_threads,input_dict['dims'])
					tw_cmd_now = tw_cmd.replace('TWPROCS',str(nprocs)).replace('TWTHREADS',str(nthreads))
					tw_cmd_now += '--no-interactive --input-file ' + ex_path
					print('Executing',tw_cmd_now,'...')
					html_doc += '<p>TW command line = <samp>'+tw_cmd_now+'</samp></p>'
					compl = subprocess.run(tw_cmd_now.split(),stdout=subprocess.PIPE,universal_newlines=True)

					# Record the results
					if compl.returncode==0:
						warnings = []
						for l in compl.stdout.splitlines():
							if l.find('WARNING')>=0:
								warnings.append(l)
						if len(warnings)>0:
							print('Completed successfully but with warnings:')
							html_doc += '<p>Completed successfully but with warnings:</p>'
							for w in warnings:
								print(w)
								html_doc += '<p><samp>' + w.replace('<','&lt;').replace('>','&gt;') + '</samp></p>'
						else:
							print('Completed successfully.')
							html_doc += '<p>Completed successfully.</p>'
						html_doc += '<p>Representative figure:</p>'
						print('Generating visualization...',end='',flush=True)
						primitive_file = cat + '-' + ex_list[i] + '-' + fig_dict['data']
						try:
							viz_file = gen_viz(primitive_file,fig_dict)
							html_doc += '<p><img src="' + viz_file + '" alt="Visualization is missing."</p>'
							print('OK')
						except OSError as err:
							print("Trouble")
							print(err)
							html_doc += '<p>There was a problem creating the visualization.</p>'
					else:
						print('Did not complete successfully.  Error report from TW follows:')
						html_doc += '<p>Did not complete successfully.  Error report from TW follows:</p>'
						for l in compl.stdout.splitlines():
							if l.find('ERROR')>=0:
								print(l)
								html_doc += '<p><samp>' + l.replace('<','&lt;').replace('>','&gt;') + '</samp></p>'
				else:
					print('  Test not requested, or not a TW input file.')
					html_doc += '<p>Test not requested, or not a TW input file.</p>'

				cleanup('*.npy')
				if len(glob.glob('twstat'))==1:
					os.remove('twstat')

	except:
		traceback.print_exc()
		print('Unrecoverable error, attempting to close report...')
		html_doc += '\n\n<h2>REPORT ENDS PREMATURELY</h2>\n\n'

	html_doc += '</body></html>'

	with open('twreport.html','w') as f:
		print(html_doc,file=f)
