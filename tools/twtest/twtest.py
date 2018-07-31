import os
import glob
import sys
import traceback
import subprocess
import datetime
import numpy as np
import re
import twutils.dvdat as dv
import twutils.plot as twplot
import matplotlib as mpl
import matplotlib.pyplot as plt

if len(sys.argv)<2:
	print('Usage: python twtest.py tw_root tw_args')
	print('tw_root = path to turboWAVE root directory.')
	print('tw_args = any arguments to pass to turboWAVE.')
	exit(1)

def cleanup(wildcarded_path):
	cleanstr = glob.glob(wildcarded_path)
	if len(cleanstr)>0:
		cleanstr.insert(0,'rm')
		subprocess.run(cleanstr)

def gen_plot(img_file,plotter,slicing_spec,slices,dyn_range=0.0,my_color_map='viridis',global_contrast=False):
	fig = plt.figure(figsize=(7,5),dpi=75)
	if len(slices)==2:
		data_slice,plot_dict = plotter.falsecolor2d(slicing_spec,slices,dyn_range)
		if global_contrast:
			minval,maxval = plotter.get_global_minmax(dyn_range)
		else:
			minval = plot_dict['vmin']
			maxval = plot_dict['vmax']
		plt.imshow(data_slice,
			origin='lower',
			aspect='auto',
			extent=plot_dict['extent'],
			vmin=minval,
			vmax=maxval,
			cmap=my_color_map)
		if dyn_range==0.0:
			plt.colorbar(label=plotter.name.split('.')[-2])
		else:
			plt.colorbar(label=r'$\log_{10}$'+plotter.name.split('.')[-2])
		plt.xlabel(plot_dict['xlabel'],fontsize=18)
		plt.ylabel(plot_dict['ylabel'],fontsize=18)
	if len(slices)==3:
		abcissa,ordinate,plot_dict = plotter.lineout(slicing_spec,slices,dyn_range)
		plt.plot(abcissa,ordinate)
		plt.xlabel(plot_dict['xlabel'],fontsize=18)
		plt.ylabel(plot_dict['ylabel'],fontsize=18)
	plt.tight_layout()
	plt.savefig(img_file)
	plt.close()

def gen_movie(mov_file,plotter,slicing_spec,slices,dyn_range=0.0,my_color_map='viridis'):
	for i in range(plotter.dims4()[0]):
		img_file = 'frame{:03d}.png'.format(i)
		c = slicing_spec.find('t')
		slices[c-2] = i
		gen_plot(img_file,plotter,slicing_spec,slices,dyn_range,my_color_map,True)
	try:
		subprocess.run(["convert","-delay","30","frame*.png",mov_file])
		cleanup('frame*.png')
	except:
		cleanup('frame*.png')
		raise OSError("The convert program from ImageMagick may not be installed.")

def gen_viz(primitive_file,file_to_plot,slicing_spec,slices,dyn_range=0.0,my_color_map='viridis'):
	'''Figure out whether a still image or movie was requested and call appropriate routine'''
	plotter = twplot.plotter(file_to_plot)
	c = slicing_spec.find('t')
	if (c==2 or c==3) and slices[c-2]==-1:
		gen_movie(primitive_file+'.gif',plotter,slicing_spec,slices,dyn_range,my_color_map)
		return primitive_file+'.gif'
	else:
		gen_plot(primitive_file+'.png',plotter,slicing_spec,slices,dyn_range,my_color_map)
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
				slicing = []
				for s in args[2].split('=')[1].split(','):
					slicing.append(int(s))
				fig_dict['slicing'] = slicing
				fig_dict['data'] = args[3]
				if len(args)>4:
					fig_dict['dyn_range'] = np.float(args[4])
				else:
					fig_dict['dyn_range'] = 0.0
	with open(ex_path) as f:
		data = f.read()
		data = comment_remover(data)
		words = data.replace("="," ").replace(","," ").replace("("," ").replace(")"," ").split()
		start_looking = False
		for i,word in enumerate(words):
			if word=='new' and words[i+1]=='grid':
				start_looking = True
			if start_looking and word=='dimensions':
				input_dict['dims'] = (int(words[i+1]),int(words[i+2]),int(words[i+3]))
				start_looking = False
	return fig_dict,input_dict

def parse_tw_args(arg_list):
	mpi_procs = 1
	omp_threads = 1
	for i,arg in enumerate(arg_list):
		if arg=='-n':
			mpi_procs = int(arg_list[i+1])
		if arg=='-c':
			omp_threads = int(arg_list[i+1])
	return mpi_procs,omp_threads

def form_tw_cmd(num_procs,num_threads,dims):
	num_dims = 0
	if dims[0]>1:
		num_dims += 1
	if dims[1]>1:
		num_dims += 1
	if dims[2]>1:
		num_dims += 1
	if num_dims==1:
		return 'tw3d -nointeractive -n '+str(num_procs*num_threads)
	else:
		return 'tw3d -nointeractive -n '+str(num_procs)+' -c '+str(num_threads)

tw_root = sys.argv[1]
if not os.path.isdir(tw_root):
	print('You specified a root directory ('+tw_root+') which does not exist.')
	exit(1)

mpl.rcParams.update({'text.usetex' : False , 'font.size' : 16})

html_doc = ''
html_doc += '<!DOCTYPE html>'
html_doc += '<html><head><title>TurboWAVE Report</title></head>'
html_doc += '<body>'

html_doc += '<h1>TurboWAVE Test Suite Report</h1>'

html_doc += '<h2>Version Information</h2>'
html_doc += '<p>Date of report : ' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
mpi_procs,omp_threads = parse_tw_args(sys.argv[2:])
html_doc += '<p>Requested MPI processes = ' + str(mpi_procs) + '</p>'
html_doc += '<p>Requested OMP threads = ' + str(omp_threads) + '</p>'

try:
	cleanup('*.gif')
	cleanup('*.png')
	cleanup('*.dvdat')
	cleanup('*.txt')

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
	html_doc += 'ERROR: something went wrong in trying to call git for versioning information.'


try:
	category_path_list = glob.glob(tw_root+'/core/examples/*')
	category_path_list = sorted(category_path_list,key=lambda s: s.split('/')[-1][0])
	category_list = []
	for s in category_path_list:
		category_list.append(s.split('/')[-1])

	for cat in category_list:
		print('=====================================')
		print('Category',cat)
		html_doc += '<h2>"' + cat + '" Subdirectory</h2>'
		ex_path_list = glob.glob(tw_root+'/core/examples/'+cat+'/*')
		ex_path_list = sorted(ex_path_list,key=lambda s: s.split('/')[-1][0])
		ex_list = []
		for s in ex_path_list:
			# Copy all the txt files in case some are required data
			if s[-4:]=='.txt':
				subprocess.run(["cp",s,"."])
			ex_list.append(s.split('/')[-1])
		for i,ex_path in enumerate(ex_path_list):
			print('--------------------------------------')
			print('Example',ex_list[i])
			html_doc += '<h3>"' + ex_list[i] + '" Example</h3>'

			subprocess.run(["cp",ex_path,"stdin"])
			fig_dict,input_dict = parse_input_file(ex_path)

			if len(fig_dict)>0:
				# Run this case
				tw_cmd = form_tw_cmd(mpi_procs,omp_threads,input_dict['dims'])
				print('Executing',tw_cmd,'...')
				html_doc += '<p>TW command line = <samp>'+tw_cmd+'</samp></p>'
				compl = subprocess.run(tw_cmd.split(),stdout=subprocess.PIPE,universal_newlines=True)

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
							html_doc += '<p><samp>' + w + '</samp></p>'
					else:
						print('Completed successfully.')
						html_doc += '<p>Completed successfully.</p>'
					html_doc += '<p>Representative figure:</p>'
					print('Generating visualization...',end='',flush=True)
					primitive_file = cat + '-' + ex_list[i] + '-' + fig_dict['data']
					try:
						viz_file = gen_viz(primitive_file,fig_dict['data'],fig_dict['slicing_spec'],fig_dict['slicing'],fig_dict['dyn_range'])
						html_doc += '<p><img src="' + viz_file + '" alt="Visualization is missing."</p>'
						print('OK')
					except OSError as err:
						print("Trouble")
						print(err)
						html_doc += '<p>There was a problem creating the visualization.  Perhaps ImagMagick is not installed.</p>'
				else:
					print('Did not complete successfully.  Error report from TW follows:')
					html_doc += '<p>Did not complete successfully.  Error report from TW follows:</p>'
					for l in compl.stdout.splitlines():
						if l.find('ERROR')>=0:
							print(l)
							html_doc += '<p><samp>' + l + '</samp></p>'
			else:
				print('  Test not requested, or not a TW input file.')
				html_doc += '<p>Test not requested, or not a TW input file.</p>'

			cleanup('*.dvdat')
			subprocess.run(["rm","-f","stdin"])
			subprocess.run(["rm","-f","twstat"])

		cleanup('*.txt')
except:
	traceback.print_exc()
	print('Unrecoverable error, attempting to close report...')
	html_doc += '<h2>REPORT ENDS PREMATURELY</h2>'

html_doc += '</body></html>'

with open('twreport.html','w') as f:
	print(html_doc,file=f)
