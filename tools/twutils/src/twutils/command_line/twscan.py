import os
import shutil
import sys
import glob
import subprocess
import re
import json
import numpy as np

# Scripted arbitrary dimensional parameter scans.
# Works either as parallel-parallel (typically HPC), or parallel-serial (typically workstation).
# Premise is to change the values of #define constants for each run.
# The values of the constants are determined by a scan file, which is a JSON serialized
# dictionary, whose keys are the identifiers, and whose values are nested lists of floats.

def print_usage():
    print('Usage: twscan <scan_file> [optional arguments...]')
    print('<scan_file> : JSON serialized dictionary of parameters and values')
    print('-n : number of MPI tasks (parallel-serial only)')
    print('-c : number of OpenMP threads per MPI task (parallel-serial only)')
    print('--dry-run : only create files, do not submit jobs')
    print('--input-file,-i : name of input file')
    print('--launch,-l: MPI launcher, if needed')
    print('--submit: job submission command, if needed')
    print('--script: job submission script, if needed')
    print('--help,-h : this message')

def decode_job(n,dims):
    '''Return index space coordinate for job n as a list'''
    stride = [1]
    for i in range(1,len(dims)):
        stride += [stride[i-1]*dims[i-1]]
    idx = list(range(len(dims)))
    idx[len(dims)-1] = int(n / stride[len(dims)-1])
    offset = 0
    for i in range(len(dims)-2,-1,-1):
        offset += idx[i+1]*stride[i+1]
        idx[i] = int((n-offset)/stride[i])
    return idx

def main():
    exec_name = 'tw3d'
    input_name = 'stdin'
    script_name = ''
    submit = ''
    launch = ''
    regex_dict = {   '' : '',
                    'qsub' : r'^#PBS\s+-N\s+',
                    'sbatch' : r'^#SBATCH\s+(--job-name|-J)\s*=\s*' }
    name_dict = {   '' : '',
                    'qsub' : r'#PBS -N ',
                    'sbatch' : r'#SBATCH --job-name=' }

    # Process command line arguments

    if len(sys.argv)<2 or '--help' in sys.argv or '-h' in sys.argv:
        print_usage()
        exit(0)

    with open(sys.argv[1],'r') as f:
        substitutions = json.load(f)
    keys = list(substitutions.keys())
    dims = np.array(substitutions[keys[0]]).shape
    njobs = 1
    for s in dims:
        njobs  = njobs * s
    print('Dimensions of scan =',dims)
    print('Parameters to vary =',keys)

    dry_run = False
    prev = ''
    mpi_tasks = 1
    omp_threads = 1

    for arg in sys.argv[2:]:
        if prev=='' and arg not in ['--dry-run','-h','--help','-n','-c','-i','--input-file','-l','--launch','--submit','--script']:
            raise ValueError('Argument '+arg+' not recognized.')
        if prev=='-n':
            mpi_tasks = int(arg)
        if prev=='-c':
            omp_threads = int(arg)
        if prev=='--input-file' or prev=='-i':
            input_name = arg
        if prev=='--launch' or prev=='-l':
            launch = arg
        if prev=='--submit':
            submit = arg
        if prev=='--script':
            script_name = arg
        if arg=='--dry-run':
            dry_run = True
        if arg in ['-n','-c','-i','--input-file','-l','--launch','--submit','--script']:
            prev = arg
        else:
            prev = ''
    name_re = regex_dict[submit]
    name_cmd = name_dict[submit]

    # Look for needed files

    if submit!='':
        check = glob.glob(exec_name)
        if len(check)==0:
            raise OSError('Required executable "'+exec_name+'" is not in the working directory.')

        check = glob.glob(script_name)
        if len(check)==0:
            raise OSError('Required batch script "'+script_name+'" is not in the working directory.')

        with open(script_name,'r') as f:
            data = f.read()
            if data.find('\ncd ')>=0:
                raise Exception('The batch script changes directories, this is not allowed.')
            if data.find('./'+exec_name)==-1:
                raise Exception('Could not find ./'+exec_name+' in batch script.')
            if data.find('\n'+launch)==-1:
                raise Exception('Could not find [LF]'+launch+' in batch script.')

    check = glob.glob(input_name)
    if len(check)==0:
        raise OSError('Required input file "'+input_name+'" is not in the working directory.')

    # Check existing directories

    rundirs = glob.glob('run*')
    if len(rundirs)>0:
        ans = ''
        while ans!='y' and ans!='n':
            ans = input('Found some files or directories in the form run*, OK to clean (y/n) ?')
        if ans=='n':
            print('STOPPED. Please run script in a directory where there are no important files of the form run*.')
            exit(1)
        for rundir in rundirs:
            shutil.rmtree(rundir)

    # Carry out the scan

    root_dir = os.getcwd()

    for n in range(njobs):
        idx = decode_job(n,dims)
        job_str = ''
        for i in idx:
            job_str += '_{:04d}'.format(i)
        rundir = root_dir+'/run'+job_str
        os.mkdir(rundir)
        shutil.copy(exec_name,rundir)

        # Modify the input file
        with open(root_dir+'/'+input_name,'r') as f:
            data = f.read()
            for key in substitutions:
                esc_key = re.sub(r'\$',r'\$',key)
                val = np.array(substitutions[key])[tuple(idx)]
                macro_def = r'^\s*#define\s+'+esc_key+r'\s+.*\S\s*$'
                if re.search(macro_def,data,flags=re.MULTILINE)==None:
                    raise Exception('no match for parameter ' + key)
                data = re.sub(macro_def,'#define '+key+' '+str(val),data,flags=re.MULTILINE)
            with open(rundir+'/'+input_name,'w') as g:
                g.write(data)

        if submit=='':
            # This is a serial batch
            if dry_run:
                print('dry run: setup',rundir,'(no execution)')
            else:
                print('Executing job',job_str)
                print('-------------------------------\n')
                os.chdir(rundir)
                subprocess.run([exec_name,'-n',str(mpi_tasks),'-c',str(omp_threads),'--no-interactive','--input-file',input_name])
                os.chdir(root_dir)
        else:
            # This is a parallel batch
            # Modify the batch script
            script_path_out = rundir+'/'+script_name
            with open(root_dir+'/'+script_name,'r') as f:
                data = f.read()
                # Add the change of directory
                new_script = re.sub(r'^\s*'+launch,'\ncd '+rundir+'\n'+launch,data,flags=re.MULTILINE)
                # Change the job name
                new_script = re.sub(name_re+r'\w+',name_cmd+'tw'+job_str,new_script,flags=re.MULTILINE)
                with open(script_path_out,'w') as g:
                    g.write(new_script)
            # Submit the job
            if dry_run:
                print('dry run: setup',rundir,'(no job submission)')
            else:
                print('Submitting '+script_path_out)
                subprocess.run([submit,script_path_out])
