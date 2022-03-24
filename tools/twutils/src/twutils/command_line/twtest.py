import os
import glob
import ntpath
import sys
import pkg_resources
import traceback
import subprocess
import datetime
import numpy as np
import re
import twutils.plot as twplot
import matplotlib as mpl
import matplotlib.pyplot as plt
import PIL.Image

class term:
	ok = '\u2713'
	err = '\u2717'
	red = '\u001b[31m'
	green = '\u001b[32m'
	yellow = '\u001b[33m'
	blue = '\u001b[94m'
	cyan = '\u001b[96m'
	reset_color = '\u001b[39;49m'
	reset_all = '\u001b[0m'
	bold = '\u001b[1m'

def print_usage():
	print('USAGE:')
	print('  Recommend running from a scratch directory.')
	print('  twtest [FLAGS] [OPTIONS] --root <tw_root> --command <cmd>')
	print()
	print('FLAGS:')
	print('  --unit         unit tests')
	print('  --integration  integration tests')
	print('  --sea-trials   full scale examples')
	print('  -h, --help     display help')
	print('  -v, --version  display version')
	print()
	print('OPTIONS:')
	print('  --categories <cat1,cat2,...>  limit to explicit categories')
	print()
	print('ARGS:')
	print('  <tw_root>   path to turboWAVE root directory.')
	print('  <cmd>       command line to execute for each test.')
	print()
	print('EXAMPLES:')
	print('  twtest --categories hydro,pic --sea-trials --root ~/turboWAVE --command tw3d -n 4 -c 5')
	print('  twtest --integration --unit --root ~/turboWAVE --command mpirun -np 4 tw3d')

def print_version():
    print('twtest is provided by twutils, version '+pkg_resources.get_distribution('twutils').version)

def cleanup(wildcarded_path):
	cleanstr = glob.glob(wildcarded_path)
	for f in cleanstr:
		os.remove(f)

def gen_plot(img_file,plotter,slicing_spec,slices,dyn_range=0.0,val_range=(0.0,0.0),my_color_map='viridis',global_contrast=False):
	plt.figure(figsize=(7,5),dpi=75)
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

def check_data(ci_dict):
	'''Read the data file and compare to expectations'''
	file_to_check = ci_dict['data']
	val_range = ci_dict['range']
	tolerance = ci_dict['tol']
	dat = np.load(file_to_check)
	if np.abs(np.max(dat[-1,...]) - val_range[1]) > tolerance:
		return 1
	if np.abs(np.min(dat[-1,...]) - val_range[0]) > tolerance:
		return 1
	return 0

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
	input_dict = {'ci' : [], 'fig': []}
	with open(ex_path) as f:
		for line in f:
			c = line.find('CITEST')
			if c>=0:
				args = line[c:].split()
				if len(args)!=4:
					raise ValueError('Example file with badly formed CITEST specification')
				input_dict['ci'] += [{}]
				input_dict['ci'][-1]['data'] = args[1]
				for s in args[2:]:
					if s.split('=')[0]=='range':
						input_dict['ci'][-1]['range'] = tuple(map(float,s.split('=')[1].split(',')))
					if s.split('=')[0]=='tolerance':
						input_dict['ci'][-1]['tol'] = float(s.split('=')[1])
				if len(input_dict['ci'][-1])!=3:
					raise ValueError('Example file with badly formed CITEST specification')
			c = line.find('TWTEST')
			if c>=0:
				args = line[c:].split()
				if args[1]!='matplotlib':
					raise ValueError('Example file with invalid plotter'+args[1])
				input_dict['fig'] += [{}]
				input_dict['fig'][-1]['slicing_spec'] = args[2].split('=')[0]
				slices = []
				for s in args[2].split('=')[1].split(','):
					slices.append(s)
				input_dict['fig'][-1]['slices'] = slices
				input_dict['fig'][-1]['data'] = args[3]
				input_dict['fig'][-1]['dyn_range'] = 0.0
				input_dict['fig'][-1]['val_range'] = (0.0,0.0)
				input_dict['fig'][-1]['color'] = 'viridis'
				input_dict['fig'][-1]['units'] = 'plasma'
				if len(args)>4:
					for keyval in args[4:]:
						key = keyval.split('=')[0]
						val = keyval.split('=')[1]
						if key=='dr':
							input_dict['fig'][-1]['dyn_range'] = np.float(val)
						if key=='range':
							r = val.split(',')
							input_dict['fig'][-1]['val_range'] = ( np.float(r[0]) , np.float(r[1]) )
						if key=='color':
							input_dict['fig'][-1]['color'] = val
						if key=='units':
							input_dict['fig'][-1]['units'] = val
	if len(input_dict['fig'])==0 and len(input_dict['ci'])==0:
		return input_dict
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
	return input_dict

def get_requested_concurrency(arg_list):
	mpi_procs = 1
	omp_threads = 1
	for i,arg in enumerate(arg_list):
		if arg=='-n' or arg=='-np':
			mpi_procs = int(arg_list[i+1])
		if arg=='-c':
			omp_threads = int(arg_list[i+1])
	return mpi_procs,omp_threads

def optimize_concurrency(num_procs,num_threads,dims):
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

def form_command(arg_list,num_procs,num_threads):
	cmd_list = []
	for arg in arg_list:
		cmd_list += [arg]
	for i,arg in enumerate(cmd_list):
		if arg=='-n' or arg=='-np':
			cmd_list[i+1] = str(num_procs)
		if arg=='-c':
			cmd_list[i+1] = str(num_threads)
	cmd = ''
	for s in cmd_list:
		cmd += s + ' '
	return cmd[:-1]

def SeaTrials(args):
	'''Runs examples with TWTEST line and generates human-readable report.
	Arguments:
	args = dictionary of preprocessed command line arguments'''

	mpl.rcParams.update({'text.usetex' : False , 'font.size' : 16})

	html_doc = ''
	html_doc += '<!DOCTYPE html>'
	html_doc += '<html><head><title>TurboWAVE Report</title></head>'
	html_doc += '<body>'

	html_doc += '<h1 style="background-color:black;color:white;">TurboWAVE Test Suite Report</h1>\n\n'

	html_doc += '<h2 style="background-color:rgb(0,20,100);color:white;">Version Information</h2>\n\n'
	html_doc += '<p>Date of report : ' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	mpi_procs,omp_threads = get_requested_concurrency(args['--command'])
	html_doc += '<p>Requested MPI processes = ' + str(mpi_procs) + '</p>'
	html_doc += '<p>Requested OMP threads = ' + str(omp_threads) + '</p>'

	cleanup('*.gif')
	cleanup('*.png')
	cleanup('*.npy')

	save_cwd = os.getcwd()
	os.chdir(args['--root'])
	try:
		print('Running Sea Trials')
		print('This may take a long time...')
		print()
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
	except FileNotFoundError:
		html_doc += '<p>ERROR: Could not find git (used to report versioning information).</p>'
	os.chdir(save_cwd)

	try:
		category_path_list = glob.glob(args['--root']+'/core/examples/*')
		category_path_list = sorted(category_path_list,key=lambda s: ntpath.basename(s)[0])
		category_list = []
		for s in category_path_list:
			category_list.append(ntpath.basename(s))
		if len(args['--categories'])>0:
			for req in args['--categories']:
				if req not in category_list:
					raise ValueError('Requested category '+req+' not found.')
			category_list = args['--categories']

		for cat in category_list:
			print('=====================================')
			print('Category',cat)
			html_doc += '\n\n<h2 style="background-color:rgb(0,20,100);color:white;">"' + cat + '" Subdirectory</h2>\n\n'
			ex_path_list = glob.glob(args['--root']+'/core/examples/'+cat+'/*')
			ex_path_list = sorted(ex_path_list,key=lambda s: ntpath.basename(s)[0])
			ex_list = []
			for s in ex_path_list:
				ex_list.append(ntpath.basename(s))
			for i,ex_path in enumerate(ex_path_list):
				print('--------------------------------------')
				print('Example',ex_list[i])
				html_doc += '\n<h3 style="background-color:rgb(100,250,200);">"' + ex_list[i] + '" Example</h3>\n'

				input_dict = parse_input_file(ex_path)

				if len(input_dict['fig'])>0:
					# Run this case
					nprocs,nthreads = optimize_concurrency(mpi_procs,omp_threads,input_dict['dims'])
					cmd = form_command(args['--command'],nprocs,nthreads)
					cmd += ' --no-interactive --input-file ' + ex_path
					print('Executing',cmd,'...')
					html_doc += '<p>TW command line = <samp>'+cmd+'</samp></p>'
					compl = subprocess.run(cmd.split(),stdout=subprocess.PIPE,universal_newlines=True)

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
						primitive_file = cat + '-' + ex_list[i] + '-' + input_dict['fig'][-1]['data']
						try:
							viz_file = gen_viz(primitive_file,input_dict['fig'][-1])
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
					cleanup('*.npy')
					cleanup('*tw_metadata*')
					cleanup('*grid_warp*')
					if len(glob.glob('twstat'))==1:
						os.remove('twstat')
				else:
					print('  Test not requested, or not a TW input file.')
					html_doc += '<p>Test not requested, or not a TW input file.</p>'


	except:
		traceback.print_exc()
		print('Unrecoverable error, attempting to close report...')
		html_doc += '\n\n<h2>REPORT ENDS PREMATURELY</h2>\n\n'

	html_doc += '</body></html>'

	with open('twreport.html','w') as f:
		print(html_doc,file=f)

def UnitTest(args):
	'''Runs turboWAVE unit tests.
	Arguments:
	args = dictionary of preprocessed command line arguments'''

	err_report = {}
	mpi_procs,omp_threads = get_requested_concurrency(args['--command'])

	try:
		print('Running all Unit Tests')
		print('----------------------')
		cmd = ['tw3d','-n','2','--unit-test','--all']
		compl = subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True,encoding='utf_8')

		# Check the results
		if compl.returncode==0:
			passing = compl.stdout.count(term.ok)
			failing = compl.stdout.count(term.err)
			print('    ' + term.ok + ' ' + term.green + 'passing' + term.reset_all)
			print('    {} passing, {} failing'.format(passing,failing))
		else:
			print('    ' + term.err + ' ' + term.red + 'failure' + term.reset_all)
			err_report['unit tests'] = {}
			err_report['unit tests']['stdout'] = compl.stdout
			err_report['unit tests']['stderr'] = compl.stderr

		# Error report
		for i,err in enumerate(err_report):
			print()
			print(str(i+1)+'. '+term.red + err + term.reset_all + ' is failing.')
			print('    stdout follows:')
			for line in err_report[err]['stdout'].splitlines():
				print('    ' + term.cyan + line + term.reset_all)
			print('    stderr follows:')
			for line in err_report[err]['stderr'].splitlines():
				print('    ' + term.yellow + line + term.reset_all)

	except:
		traceback.print_exc()
		print('unrecoverable error')
		return 1
	
	return len(err_report)

def CITest(args):
	'''Runs tests with CITEST line and spot checks results.
	Arguments:
	args = dictionary of preprocessed command line arguments'''

	err_report = {}
	mpi_procs,omp_threads = get_requested_concurrency(args['--command'])

	cleanup('*.gif')
	cleanup('*.png')
	cleanup('*.npy')

	try:
		print('Running Integration Tests')
		print('-------------------------')
		category_path_list = glob.glob(args['--root']+'/core/test/*')
		category_path_list = sorted(category_path_list,key=lambda s: ntpath.basename(s)[0])
		category_list = []
		for s in category_path_list:
			category_list.append(ntpath.basename(s))
		if len(args['--categories'])>0:
			for req in args['--categories']:
				if req not in category_list:
					raise ValueError('Requested category '+req+' not found.')
			category_list = args['--categories']

		for cat in category_list:
			print('Category',cat)
			ex_path_list = glob.glob(args['--root']+'/core/test/'+cat+'/*')
			ex_path_list = sorted(ex_path_list,key=lambda s: ntpath.basename(s)[0])
			ex_list = []
			for s in ex_path_list:
				ex_list.append(ntpath.basename(s))
			for i,ex_path in enumerate(ex_path_list):

				input_dict = parse_input_file(ex_path)

				if len(input_dict['ci'])>0:
					# Run this case
					nprocs,nthreads = optimize_concurrency(mpi_procs,omp_threads,input_dict['dims'])
					cmd = form_command(args['--command'],nprocs,nthreads)
					cmd += ' --no-interactive --input-file ' + ex_path
					compl = subprocess.run(cmd.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)

					# Check the results
					if compl.returncode==0:
						try:
							if check_data(input_dict['ci'][-1])==1:
								print('    ' + term.err + ' ' + term.red + ex_list[i] + term.reset_all)
								err_report[ex_list[i]] = {}
								err_report[ex_list[i]]['stdout'] = ''
								err_report[ex_list[i]]['stderr'] = input_dict['ci'][-1]['data'] + ' extremum delta exceeds tolerance'
							else:
								print('    ' + term.ok+ ' ' + term.green + ex_list[i] + term.reset_all)
						except OSError as err:
							err_report[ex_list[i]] = {}
							err_report[ex_list[i]]['stdout'] = ''
							err_report[ex_list[i]]['stderr'] = err
					else:
						print('    ' + term.err + ' ' + term.red + ex_list[i] + term.reset_all + ' (run failed)')
						err_report[ex_list[i]] = {}
						err_report[ex_list[i]]['stdout'] = compl.stdout
						err_report[ex_list[i]]['stderr'] = compl.stderr

					cleanup('*.npy')
					cleanup('*tw_metadata*')
					cleanup('*grid_warp*')
					if len(glob.glob('twstat'))==1:
						os.remove('twstat')

		# Error report
		for i,err in enumerate(err_report):
			print()
			print(str(i+1)+'. '+term.red + err + term.reset_all + ' is failing.')
			print('    stdout follows:')
			for line in err_report[err]['stdout'].splitlines():
				print('    ' + term.cyan + line + term.reset_all)
			print('    stderr follows:')
			for line in err_report[err]['stderr'].splitlines():
				print('    ' + term.yellow + line + term.reset_all)

	except:
		traceback.print_exc()
		print('unrecoverable error')
		return 1
	
	return len(err_report)

def panic(mess):
	print(term.red+'ERROR: '+term.reset_all+mess)
	print('suggestion '+term.bold+term.green+'twtest --help'+term.reset_all)
	exit(1)

def main():

	args = { '-v': None,
		'--version' : None,
		'-h' : None,
		'--help' : None,
		'--unit': False,
		'--integration': False,
		'--sea-trials': False,
		'--categories': [],
		'--root': '',
		'--command': [],
		'expecting_value': None
	}

	for arg in sys.argv[1:]:
		if arg not in args and args['expecting_value']==None:
			panic('argument '+arg+' not recognized')
		if args['expecting_value']!=None:
			if args['expecting_value']=='--categories':
				args['--categories'] = arg.split(',')
				args['expecting_value'] = None
			if args['expecting_value']=='--root':
				args['--root'] = arg
				args['expecting_value'] = None
			if args['expecting_value']=='--command':
				args['--command'] += [arg]
				# Consumes remaining arguments
		else:
			if arg=='--help' or arg=='-h':
				print_usage()
				exit(0)
			if arg=='--version' or arg=='-v':
				print_version()
				exit(0)
			if arg in ['--integration','--unit','--sea-trials']:
				args[arg] = True

			if arg in ['--categories','--root','--command']:
				args['expecting_value'] = arg

	if args['--root']=='':
		panic('root directory not specified')
	if len(args['--command'])==0:
		panic('command not specified')
	if not os.path.isdir(args['--root']):
		panic('requested root directory ('+args['--root']+') does not exist.')

	errCount = 0
	if args['--unit']:
		errCount += UnitTest(args)
	if args['--integration']:
		print()
		errCount += CITest(args)
	if args['--sea-trials']:
		print()
		SeaTrials(args)

	return errCount