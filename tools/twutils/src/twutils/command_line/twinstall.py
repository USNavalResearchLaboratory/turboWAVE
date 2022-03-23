'''This program is used to interactively install turboWAVE core once
the tools have already been installed.  This can be run in two modes:
1. graphical interface based on the tkinter module (default).
2. textual interface based on the curses module (if --terminal argument is given).
Both interfaces are built on the 3 objects config, tasks, cmd.
The config object displays the configuration.
The tasks object displays the list of tasks.
The cmd object accepts inputs, and displays messages.'''

import time
import os
import sys
import shutil
import pathlib
import glob
import subprocess
import re
import copy
import pkg_resources
try:
    import curses
    has_curses = True
except ModuleNotFoundError:
    has_curses = False
try:
    import tkinter as tk
    from tkinter import ttk
    import tkinter.messagebox
    import tkinter.filedialog
    has_tk = True
except ModuleNotFoundError:
    has_tk = False

# Stuff for terminal mode
enter_keys = [ord('\n'),ord('\r')]
panel_width = 80

source_url = 'https://github.com/USNavalResearchLaboratory/turboWAVE.git'
incompatibility_matrix = {
    'Windows' : ['GNU','Cray','CUDA','Apple CL','AMD ROCm','Homebrew','MacPorts'],
    'Linux' : ['Cray','Apple CL','Homebrew','MacPorts'],
    'MacOS' : ['Cray','Intel','CUDA','AMD ROCm','None'],
    'Cray' : ['GNU','LLVM','Apple CL','AMD ROCm','Homebrew','MacPorts'],
    'OpenMPI Cluster' : ['Cray','Apple CL','Homebrew','MacPorts'],
    'IntelMPI Cluster' : ['Cray','Apple CL','Homebrew','MacPorts']
}
help_str = '''This program configures, compiles, and installs turboWAVE core.
Press the buttons in the "Tasks" frame in order from top to bottom.
This will lead you through the necessary steps.'''
config_affirm='''The dropdown menus in the "Current Configuration" frame can be used to edit the configuration.
We have tried to choose the best options by detecting the environment.
Are you satisfied with the selections?'''

###########################################################
### BASE CLASSES FOR EITHER TEXT OR GRAPHICAL INTERFACE ###
###########################################################

def git_get_default_tag(package_path,taglist):
    '''Checks to see if the workspace is identical to any tagged commit.
    If it is return that tag, otherwise return "untagged".'''
    ans = 'untagged'
    save_dir = os.getcwd()
    os.chdir(package_path)
    for tag in taglist:
        compl = subprocess.run(["git","diff",tag],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
        if len(compl.stdout)<2:
            ans = tag
    os.chdir(save_dir)
    return ans

def git_is_clean(package_path):
    save_dir = os.getcwd()
    os.chdir(package_path)
    compl = subprocess.run(["git","status"],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
    os.chdir(save_dir)
    return 'working tree clean' in compl.stdout.lower()

def git_is_detached(package_path):
    save_dir = os.getcwd()
    os.chdir(package_path)
    compl = subprocess.run(["git","status"],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
    os.chdir(save_dir)
    return 'head detached' in compl.stdout.lower()

def git_err_handler(compl,cmd):
    if compl.returncode!=0:
        match = re.search(r'\b(err|ERR|Err).*\n',compl.stdout)
        if match:
            err_str = match.group(0)
        else:
            err_str = '(we cannot find the error string)'
        if cmd.affirm('Git returned an error:\n'+err_str+'\nWould you like to continue anyway?')=='y':
            return False
        else:
            return True
    else:
        return False

class dictionary_view:
    def __init__(self):
        self.highlighted = -1
        self.trouble = []
        self.title = 'Selection'
        self.data = {}
        self.choices = {}
        self.unpretty = {}
        self.widgets = {}
        self.selectable = []
    def highlight(self,i0):
        if i0==-1:
            self.highlighted = -1
            return
        i1 = 0
        for i,key in enumerate(self.data):
            if key in self.selectable:
                if i1==i0:
                    break
                i1 += 1
        self.highlighted = i
    def set_trouble(self,key):
        self.trouble += [key]
    def set_ok(self,key=None):
        if type(key)==type(None):
            self.trouble = []
        else:
            self.trouble = list(set(self.trouble)-{key})
    def get(self,key):
        return self.data[key]
    def verify(self):
        self.set_ok()
        return len(self.trouble)

class base_config(dictionary_view):
    def __init__(self):
        super().__init__()
        self.title = 'Current Configuration'
        self.unpretty['None'] = 'NONE'
        self.unpretty['Linux'] = 'LINUX'
        self.unpretty['Windows'] = 'WIN'
        self.unpretty['MacOS'] = 'OSX'
        self.unpretty['Cray'] = 'CRAY'
        self.unpretty['OpenMPI Cluster'] = 'OPENMPI'
        self.unpretty['IntelMPI Cluster'] = 'INTELMPI'
        self.unpretty['GNU'] = 'GNU'
        self.unpretty['LLVM'] = 'LLVM_CLANG'
        self.unpretty['Intel'] = 'INTEL'
        self.unpretty['Cray'] = 'CRAY'
        self.unpretty['OpenMP'] = 'OMP'
        self.unpretty['CUDA'] = 'CUDA'
        self.unpretty['Apple CL'] = 'APPLE_CL'
        self.unpretty['AMD ROCm'] = 'RADEON_PRO'
        self.unpretty['Homebrew'] = 'HOMEBREW'
        self.unpretty['MacPorts'] = 'MACPORTS'
        self.unpretty['Scalar'] = '64'
        self.unpretty['SSE'] = '128'
        self.unpretty['SSE2'] = '128'
        self.unpretty['AVX'] = '256'
        self.unpretty['AVX2'] = '256'
        self.unpretty['AVX512'] = '512'
        self.unpretty['Little Endian'] = 'LITTLE'
        self.unpretty['Big Endian'] = 'BIG'
        self.choices['Platform'] = ['Linux','Windows','MacOS','Cray','OpenMPI Cluster','IntelMPI Cluster']
        self.choices['Compiler'] = ['GNU','LLVM','Intel','Cray']
        self.choices['Accelerator'] = ['OpenMP','CUDA','Apple CL','AMD ROCm']
        self.choices['Packager'] = ['None','Homebrew','MacPorts']
        # the list of vectors is assumed to be in order of ascending preference
        self.choices['Vectors'] = ['Scalar','SSE','SSE2','AVX','AVX2','AVX512']
        self.choices['Byte Order'] = ['Little Endian','Big Endian']
        self.data['Tools Version'] = pkg_resources.get_distribution('twutils').version
        self.data['Core Version'] = ''
        self.data['Build Path'] = ''
        self.data['Platform'] = self.choices['Platform'][0]
        self.data['Compiler'] = self.choices['Compiler'][0]
        self.data['Accelerator'] = self.choices['Accelerator'][0]
        self.data['Vectors'] = self.choices['Vectors'][0]
        self.data['Packager'] = self.choices['Packager'][0]
        self.data['Byte Order'] = self.choices['Byte Order'][0]
        self.selectable = ['Platform','Compiler','Accelerator','Vectors','Packager','Byte Order']
        # Try to initialize the configuration based on environment
        if 'linux' in sys.platform:
            self.data['Platform'] = 'Linux'
            self.data['Compiler'] = 'LLVM'
            compl = subprocess.run(['lscpu'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
            for v in reversed(self.choices['Vectors']):
                if re.search(r'\b'+v,compl.stdout,flags=re.IGNORECASE)!=None:
                    self.data['Vectors'] = v
                    break
        if 'darwin' in sys.platform:
            self.data['Platform'] = 'MacOS'
            self.data['Compiler'] = 'GNU'
            self.data['Packager'] = 'Homebrew'
            compl = subprocess.run(['sysctl','machdep.cpu'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
            for v in reversed(self.choices['Vectors']):
                if re.search(r'\b'+v,compl.stdout,flags=re.IGNORECASE)!=None:
                    self.data['Vectors'] = v
                    break
        if 'win32' in sys.platform:
            self.data['Platform'] = 'Windows'
            self.data['Compiler'] = 'LLVM'
            self.data['Vectors'] = 'AVX2'
        if sys.byteorder=='big':
            self.data['Byte Order'] = 'Big Endian'
    def connect(self,tasks):
        self.tasks = tasks
    def backup_params(self):
        self.backup = copy.copy(self.data)
    def restore_params(self):
        self.data = copy.copy(self.backup)
    def get_make_val(self,key):
        return self.unpretty[self.data[key]]
    def get_params(self):
        ans = []
        for key in self.data:
            if key in self.choices:
                ans += [key]
        return ans
    def get_choices(self,key):
        return self.choices[key]
    def select(self,key,choice):
        self.set(key,choice)
        self.tasks.dirty_config()
    def verify(self):
        self.set_ok()
        if self.data['Compiler'] in incompatibility_matrix[self.data['Platform']]:
            self.set_trouble('Compiler')
        if self.data['Packager'] in incompatibility_matrix[self.data['Platform']]:
            self.set_trouble('Packager')
        if self.data['Accelerator'] in incompatibility_matrix[self.data['Platform']]:
            self.set_trouble('Accelerator')
        return len(self.trouble)

class base_tasks(dictionary_view):
    def __init__(self):
        super().__init__()
        self.title = 'Tasks'
        self.highlighted = 0
        self.data['Get Components'] = 'not started'
        self.data['Set Version'] = 'not started'
        self.data['Configure Makefile'] = 'not started'
        self.data['Compile Code'] = 'not started'
        self.data['Install Files'] = 'not started'
        self.data['Cleanup'] = 'not started'
        self.selectable = [key for key in self.data]
    def connect(self,cmd,conf):
        self.cmd = cmd
        self.conf = conf
    def dirty_config(self):
        if self.get('Configure Makefile')=='done':
            self.set('Configure Makefile','incomplete')
        if self.get('Compile Code')=='done':
            self.set('Compile Code','incomplete')
        if self.get('Install Files')=='done':
            self.set('Install Files','incomplete')
        if self.get('Cleanup')=='done':
            self.set('Cleanup','incomplete')
    def finished(self):
        for key in self.data:
            if self.data[key]!='done':
                return False
        return True
    def verify(self):
        self.set_ok()
        for key in self.data:
            if self.data[key]=='incomplete':
                self.set_trouble(key)
        return len(self.trouble)
    def GetComponents(self):
        self.set('Get Components','incomplete')
        test = []
        primitive_path = self.AskDirectory()
        if type(primitive_path)==tuple:
            return
        self.build_path = pathlib.Path(primitive_path) / 'turboWAVE' / 'core' / 'build'
        self.package_path = pathlib.Path(primitive_path) / 'turboWAVE'
        verified_path = glob.glob(str(self.build_path))
        if len(verified_path)==0:
            self.build_path = pathlib.Path(primitive_path) / 'core' / 'build'
            self.package_path = pathlib.Path(primitive_path)
            verified_path = glob.glob(str(self.build_path))
        if len(verified_path)==0:
            ans = self.cmd.affirm('Path is '+primitive_path+'\nTW source not found, would you like to retrieve it?')
            if ans=='y':
                dest = str(pathlib.Path(primitive_path) / 'turboWAVE')
                self.cmd.display('Cloning repository...')
                compl = subprocess.run(["git","clone",source_url,dest],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
                self.build_path = pathlib.Path(dest) / 'core' / 'build'
                self.package_path = pathlib.Path(dest)
                verified_path = glob.glob(str(self.build_path))
                if git_err_handler(compl,self.cmd):
                    verified_path = []
                self.cmd.display('')
        else:
            if git_is_clean(self.package_path):
                ans = self.cmd.affirm('TW source found, and working tree clean.\nWould you like to pull the latest?')
                if ans=='y':
                    save_dir = os.getcwd()
                    os.chdir(verified_path[0])
                    self.cmd.display('Pulling from repository...')
                    compl = subprocess.run(["git","pull","origin","master"],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
                    if git_err_handler(compl,self.cmd):
                        verified_path = []
                    compl = subprocess.run(["git","pull","origin","--tag"],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
                    if git_err_handler(compl,self.cmd):
                        verified_path = []
                    self.cmd.display('')
                    os.chdir(save_dir)
            else:
                self.cmd.info('Warning','TW source found, but working tree not clean.\nYou can continue if you know what you are doing.')
        if len(verified_path)>0:
            self.conf.set('Build Path',verified_path[0])
            self.set('Get Components','done')
        else:
            self.conf.set('Build Path','')
            self.cmd.err('Get components did not succeed.')

    def StartSetVersion(self):
        if self.get('Get Components')!='done':
            self.cmd.err('Please complete <Get Components> first.')
            return
        self.set('Set Version','incomplete')
        save_dir = os.getcwd()
        os.chdir(self.package_path)
        compl = subprocess.run(["git","tag","-l"],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
        os.chdir(save_dir)
        if git_err_handler(compl,self.cmd)==False:
            rawlist = compl.stdout.splitlines()
            stable = []
            major = []
            prerelease = []
            for item in rawlist:
                if 'v' in item:
                    major += [item]
                elif 'a' in item or 'b' in item or 'rc' in item:
                    prerelease += [item]
                else:
                    stable += [item]
            taglist = sorted(stable[-5:])[::-1] + sorted(prerelease[-5:])[::-1]
            self.default_tag = 'workspace (' + git_get_default_tag(self.package_path,taglist) + ')'
            taglist = [self.default_tag] + taglist
            self.popup(taglist,self.FinishSetVersion)
    def FinishSetVersion(self,tag):
        if tag=='escape':
            return
        if tag==self.default_tag:
            self.conf.set('Core Version',self.default_tag)
            self.set('Set Version','done')
        else:
            save_dir = os.getcwd()
            os.chdir(self.package_path)
            compl = subprocess.run(["git","checkout",tag],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
            os.chdir(save_dir)
            if git_err_handler(compl,self.cmd)==False:
                self.conf.set('Core Version',tag)
                self.set('Set Version','done')
                if git_is_detached(self.package_path):
                    self.cmd.info('Notice','Local repo in detached state.  This can be fixed during cleanup.')
        toolsv = self.conf.data['Tools Version'].split('.')
        if '(' in tag:
            corev = ((tag.split('(')[1]).split(')')[0]).split('.')
        else:
            corev = tag.split('.')
        if toolsv[0]!=corev[0] or toolsv[1]!=corev[1]:
            self.cmd.info('Notice','Core and tool version mismatch.  For best results match core and tools to within a minor version number.')

    def Configure(self):
        if self.get('Set Version')!='done':
            self.cmd.err('Please complete <Set Version> first.')
            return
        self.set('Configure Makefile','incomplete')
        if self.cmd.affirm(config_affirm)=='n':
            return
        if self.conf.verify()>0:
            self.cmd.err('There are features selected that are incompatible with the platform.')
            return
        maker = 'GNU'
        if self.conf.get('Platform')=='Windows' and self.conf.get('Compiler')=='Intel':
            maker = 'MS'
        if maker=='MS':
            makefile_in = str(self.build_path / 'win.make')
            makefile_out = str(self.build_path / 'win.make.edited')
        if maker=='GNU':
            makefile_in = str(self.build_path / 'makefile')
            makefile_out = str(self.build_path / 'makefile.edited')
        with open(makefile_in,'r') as f:
            makefile = f.read()
        key_on = ['Platform','Compiler','Accelerator','Vectors','Packager','Byte Order']
        unpretty = {'Platform':'PLATFORM',
            'Compiler':'COMPILER_PREF',
            'Accelerator':'HARDWARE_ACCEL',
            'Vectors':'VBITS',
            'Packager':'PACKAGE_PREF',
            'Byte Order':'ENDIANNESS'}
        not_found = []
        for k in key_on:
            var = unpretty[k]
            val = self.conf.get_make_val(k)
            if re.search('\n[#\s]*' + var + '\s*=\s*',makefile)==None:
                not_found += [var]
            else:
                # Comment out all the values
                makefile = re.sub('\n[#\s]*'+var,'\n#'+var,makefile)
                # Don't care about preserving unused values, safest to just overwrite the first value found.
                makefile = re.sub('\n[#\s]*'+var+'\s*=\s*[\w\d]+','\n'+var+' = '+val,makefile,count=1)
        if self.conf.get('Packager')=='Homebrew' and self.conf.get('Compiler')=='GNU':
            gcc_glob_list = ['/usr/local/opt/**/g++-*','/opt/homebrew/**/g++-*']
            gcc_re = r'\/usr\/local\/opt\/gcc\/bin\/g\+\+-\d+'
            brew_gcc = []
            for gcc_glob in gcc_glob_list:
                brew_gcc += glob.glob(gcc_glob,recursive=True)
            if len(brew_gcc)==0:
                self.cmd.err('Did not find GCC of the form '+gcc_glob_list)
                return
            gcc_sel = brew_gcc[0]
            if len(brew_gcc)>1:
                brew_callback = lambda sel : (gcc_sel := sel,)[0]
                self.popup(brew_gcc,brew_callback)
            if re.search(gcc_re,makefile)==None:
                not_found += [gcc_sel]
            else:
                makefile = re.sub(gcc_re,gcc_sel,makefile)
        if len(not_found)==0:
            with open(makefile_out,'w') as f:
                f.write(makefile)
                self.set('Configure Makefile','done')
        else:
            missing = ''
            for s in not_found:
                missing += s + ','
            self.cmd.err('Did not find variables:'+missing[:-1]+'\nThis probably means core and tools are out of sync.\nYou might try a different version.')

    def Compile(self,enter_shell,leave_shell):
        if self.get('Configure Makefile')!='done':
            self.cmd.err('Please complete <Configure Makefile> first')
            return
        self.set('Compile Code','incomplete')
        if self.cmd.affirm('Compiler progress will appear in shell window.\nPlease wait for completion.\nReady to proceed?')=='y':
            self.cmd.display('Compiling...')
            save_dir = os.getcwd()
            os.chdir(self.build_path)
            maker = 'GNU'
            if self.conf.get('Platform')=='Windows' and self.conf.get('Compiler')=='Intel':
                maker = 'MS'
            enter_shell()
            try:
                if maker=='GNU':
                    compl = subprocess.run(['make','-f','makefile.edited','clean'],universal_newlines=True)
                    compl = subprocess.run(['make','-f','makefile.edited','tw3d_release'],universal_newlines=True)
                    if compl.returncode==0:
                        self.set('Compile Code','done')
                if maker=='MS':
                    compl = subprocess.run(['nmake','/F','win.make.edited','clean'],universal_newlines=True)
                    compl = subprocess.run(['nmake','/F','win.make.edited','tw3d'],universal_newlines=True)
                    if compl.returncode==0:
                        self.set('Compile Code','done')
                if self.get('Compile Code')!='done':
                    raise RuntimeError
            except FileNotFoundError:
                leave_shell()
                self.cmd.display('')
                if maker=='MS':
                    self.cmd.err('nmake is not supported in the synced core version.')#self.cmd.err('Could not find nmake. Try running twinstall from Intel shell.')
                else:
                    self.cmd.err('Could not find make, perhaps it is not installed?')
            except RuntimeError:
                leave_shell()
                self.cmd.display('')
                if maker=='MS':
                    self.cmd.err('nmake is not supported in the synced core version.')
                else:
                    self.cmd.err('There was an error while trying to compile.')
            except Exception as ex:
                leave_shell()
                self.cmd.display('')
                self.cmd.err('Unexpected error '+type(ex).__name__+', force exit.')
                exit(1)
            else:
                leave_shell()
                self.cmd.display('')
            os.chdir(save_dir)

    def Install(self):
        if self.get('Compile Code')!='done':
            self.cmd.err('Please complete <Compile Code> first.')
            return
        self.set('Install Files','incomplete')
        if self.conf.get('Platform')=='Windows':
            if self.conf.get('Compiler')=='Intel':
                src_exe = 'tw3d.exe'
            else:
                src_exe = 'tw3d'
            dst_exe = 'tw3d.exe'
            install_path = pathlib.Path(os.environ['CONDA_PREFIX']) / 'Scripts' / dst_exe
        else:
            src_exe = 'tw3d'
            dst_exe = 'tw3d'
            install_path = pathlib.Path(os.environ['CONDA_PREFIX']) / 'bin' / dst_exe
        repo_path = self.build_path / src_exe
        # tw3d executable
        if self.conf.get('Platform')=='Cray':
            dv_dir = self.AskDirectory('Select executable install location',str(pathlib.Path.home()))
            if type(dv_dir)!=tuple:
                shutil.copy(repo_path,dv_dir)
                self.set('Install Files','done')
        else:
            if self.cmd.affirm('We will install the main executable ('+dst_exe+') into '+str(install_path)+'. Proceed?')=='y':
                shutil.copy(repo_path,install_path)
                self.set('Install Files','done')
        # Jupyter DataViewer styles
        if self.cmd.affirm('Can we change Jupyter (increase notebook width)?')=='y':
            repo_path = self.package_path / 'tools' / 'config-files' / 'custom.css'
            dest_path = pathlib.Path.home() / '.jupyter' / 'custom'
            if len(glob.glob(str(dest_path)))==0:
                os.makedirs(dest_path,exist_ok=True)
            shutil.copy(repo_path,dest_path)
        # Jupyter DataViewer
        if self.cmd.affirm('Would you like make a copy of DataViewer?')=='y':
            repo_path = self.package_path / 'tools' / 'DataViewer.ipynb'
            dv_dir = self.AskDirectory('Select DataViewer install location',str(pathlib.Path.home()))
            if type(dv_dir)!=tuple:
                shutil.copy(repo_path,dv_dir)
        # MICRO highlighting
        if self.cmd.affirm('Would you like to install syntax highlights for MICRO?')=='y':
            repo_path = self.package_path / 'tools' / 'config-files' / 'turbowave.micro.yaml'
            dest_path = pathlib.Path.home() / '.config' / 'micro' / 'syntax'
            if len(glob.glob(str(dest_path)))==0:
                os.makedirs(dest_path,exist_ok=True)
            shutil.copy(repo_path,dest_path)
        # NANO highlighting
        if self.cmd.affirm('Would you like to install syntax highlights for NANO?')=='y':
            repo_path = self.package_path / 'tools' / 'config-files' / '.turbowave.nanorc'
            dest_path = pathlib.Path.home()
            shutil.copy(repo_path,dest_path)

            if self.conf.get('Platform')=='Windows':
                dest_path = pathlib.Path.home() / 'nano.rc'
            else:
                dest_path = pathlib.Path.home() / '.nanorc'
            try:
                with open(dest_path,'r') as f:
                    nanorc = f.read()
            except FileNotFoundError:
                nanorc = '# nano settings file, created by twinstall\n\n'
            inc_path = pathlib.Path.home() / '.turbowave.nanorc'
            search_str = r'include.+turbowave\.nanorc.*'
            inc_str = 'include ' + '"' + str(inc_path) + '"'
            if re.search(search_str,nanorc)==None:
                nanorc += '\n' + inc_str
            else:
                nanorc = re.sub(search_str,inc_str,nanorc)
            with open(dest_path,'w') as f:
                f.write(nanorc)
        # VIM highlighting
        if self.cmd.affirm('Would you like to install syntax highlights for VIM? (overwrites other filetype customizations)')=='y':
            repo_path = self.package_path / 'tools' / 'config-files' / 'filetype.vim'
            if self.conf.get('Platform')=='Windows':
                dest_path = pathlib.Path.home() / 'vimfiles'
            else:
                dest_path = pathlib.Path.home() / '.vim'
            if len(glob.glob(str(dest_path)))==0:
                os.makedirs(dest_path,exist_ok=True)
            shutil.copy(repo_path,dest_path)

            repo_path = self.package_path / 'tools' / 'config-files' / 'turbowave.vim'
            if self.conf.get('Platform')=='Windows':
                dest_path = pathlib.Path.home() / 'vimfiles' / 'syntax'
            else:
                dest_path = pathlib.Path.home() / '.vim' / 'syntax'
            if len(glob.glob(str(dest_path)))==0:
                os.makedirs(dest_path,exist_ok=True)
            shutil.copy(repo_path,dest_path)
        # Other highlighting
        self.cmd.info('Notice','Syntax highlights are available for Atom and VS Code, search for "turbowave" in the respective package managers.')

    def Clean(self):
        if self.get('Install Files')!='done':
            self.cmd.err('Please complete <Install Files> first.')
            return
        self.set('Cleanup','incomplete')
        if self.conf.get('Platform')=='Windows' and self.conf.get('Compiler')=='Intel':
            makefile_name = 'win.make.edited'
            obj_ext = '.obj'
            exe_name = 'tw3d.exe'
        else:
            makefile_name = 'makefile.edited'
            obj_ext = '.o'
            exe_name = 'tw3d'
        d = self.build_path / makefile_name
        if self.cmd.affirm('Delete '+str(d)+'?')=='y':
            os.remove(str(d))
        d = self.build_path / ('*' + obj_ext)
        if self.cmd.affirm('Delete '+str(d)+'?')=='y':
            obj_list = glob.glob(str(d))
            for f in obj_list:
                os.remove(f)
        d = self.build_path / exe_name
        if self.cmd.affirm('Delete '+str(d)+'?')=='y':
            os.remove(str(d))
        d = self.build_path / '*.d'
        if self.cmd.affirm('Delete '+str(d)+'?')=='y':
            dep_list = glob.glob(str(d))
            for f in dep_list:
                os.remove(f)
        if git_is_detached(self.package_path):
            if self.cmd.affirm('Restore local repository to master?')=='y':
                save_dir = os.getcwd()
                os.chdir(self.package_path)
                compl = subprocess.run(["git","checkout","master"],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
                os.chdir(save_dir)
                git_err_handler(compl,self.cmd)
        self.set('Cleanup','done')
        if '--terminal' in sys.argv or '-t' in sys.argv:
            self.cmd.display('All tasks are finished, you can press q to quit.')
        else:
            self.cmd.display('All tasks are finished, you can press <Quit>.')


######################
### GRAPHICAL MODE ###
######################

class gui_popup:
    def __init__(self,master_win,master_frame,items,completion_callback):
        self.window = tk.Toplevel()
        self.window.title('Selection')
        self.window.transient(master_frame)
        w = master_frame.winfo_width()
        h = master_frame.winfo_height()
        x = master_win.winfo_x() + master_frame.winfo_x()
        y = master_win.winfo_y() + master_frame.winfo_x()
        self.window.geometry('+'+str(int(w/4+x))+'+'+str(int(h/4+y)))
        self.window.grab_set()
        self.completion_callback = completion_callback
        mb = ttk.Menubutton(self.window,text='Select Version')
        mb.grid(column=0,row=0,sticky='ew',padx=5,pady=5)
        mb.menu = tk.Menu(mb,tearoff=0)
        mb['menu'] = mb.menu
        for choice in items:
            cmd = lambda choice=choice : self.select(choice)
            mb.menu.add_command(label=choice,command=cmd)
    def select(self,choice):
        self.completion_callback(choice)
        self.window.grab_release()
        self.window.destroy()

class gui_dictionary_view(dictionary_view):
    def __init__(self):
        super().__init__()
    def set(self,key,val):
        self.data[key] = val
        self.widgets[key][1].configure(text=val)
        self.verify()
        for key in self.data:
            if key in self.trouble:
                self.widgets[key][1].configure(foreground='red')
            else:
                self.widgets[key][1].configure(foreground='black')
    def popup(self,items,callback):
        return gui_popup(self.window.winfo_toplevel(),self.window,items,callback)

class gui_config(base_config,gui_dictionary_view):
    def __init__(self,master):
        super().__init__()
        self.window = ttk.LabelFrame(master,text=self.title)
        self.window.grid(row=0,column=0,padx=5,pady=5,sticky='ew')
        for i,key in enumerate(self.data):
            if key in self.choices:
                mb = ttk.Menubutton(self.window,text=key)
                mb.grid(column=0,row=i,sticky='ew',padx=5,pady=5)
                mb.menu = tk.Menu(mb,tearoff=0)
                mb['menu'] = mb.menu
                for choice in self.choices[key]:
                    cmd = lambda key=key,choice=choice : self.select(key,choice)
                    mb.menu.add_command(label=choice,command=cmd)
                var_label = ttk.Label(self.window,text=self.data[key],width=60,anchor='w')
                var_label.grid(column=1,row=i,padx=5,pady=5)
                self.widgets[key] = (mb,var_label)
            else:
                static_label = ttk.Label(self.window,text=key+':',anchor='e')
                static_label.grid(column=0,row=i,padx=5,pady=5,sticky='e')
                var_label = ttk.Label(self.window,text=self.data[key],width=60,anchor='w')
                var_label.grid(column=1,row=i,padx=5,pady=5)
                self.widgets[key] = (static_label,var_label)

class gui_tasks(base_tasks,gui_dictionary_view):
    def __init__(self,master):
        super().__init__()
        self.window = ttk.LabelFrame(master,text=self.title)
        self.window.grid(row=1,column=0,padx=5,pady=5,sticky='ew')
        for i,key in enumerate(self.data):
            button = ttk.Button(self.window,text=key)
            status = ttk.Label(self.window,text=self.data[key],width=60)
            button.grid(row=i,column=0,padx=5,pady=5,sticky='ew')
            status.grid(row=i,column=1,padx=5,pady=5)
            self.widgets[key] = (button,status)
        self.widgets['Get Components'][0].configure(command=self.GetComponents)
        self.widgets['Set Version'][0].configure(command=self.StartSetVersion)
        self.widgets['Configure Makefile'][0].configure(command=self.Configure)
        self.widgets['Compile Code'][0].configure(command=lambda : self.Compile(self.enter_shell,self.exit_shell))
        self.widgets['Install Files'][0].configure(command=self.Install)
        self.widgets['Cleanup'][0].configure(command=self.Clean)
    def enter_shell(self):
        print('Start Shell Output.')
    def exit_shell(self):
        print('Reached Completion. Stop Shell Output.')
    def AskDirectory(self,title='Select Working Directory',initialdir=str(pathlib.Path.home())):
        return tkinter.filedialog.askdirectory(title=title,initialdir=initialdir)

class gui_command:
    def __init__(self,master):
        self.window = ttk.Frame(master)
        self.window.grid(row=2,column=0,padx=5,pady=5)
        self.help = ttk.Button(self.window,text='Help',command=self.main_help)
        self.help.grid(row=0,column=0,padx=5,pady=5)
        self.quit = ttk.Button(self.window,text='Quit',command=master.quit)
        self.quit.grid(row=0,column=1,padx=5,pady=5)
        self.mess = ttk.Label(self.window)
        self.mess.grid(row=1,column=0,columnspan=2)
    def err(self,message):
        tk.messagebox.showerror(title='Notice',message=message)
    def affirm(self,message):
        affirmed = tk.messagebox.askyesno('Confirm',message)
        if affirmed:
            return 'y'
        else:
            return 'n'
    def display(self,message=''):
        self.mess.configure(text=message,font=('Times','14','bold'))
        self.mess.update()
    def info(self,title,message):
        tk.messagebox.showinfo(title=title,message=message)
    def main_help(self):
        tk.messagebox.showinfo(title='TurboWAVE Installer Help',message=help_str)

class Application(tk.Frame):
    def __init__(self,master=None):
        super().__init__(master)
        self.master = master
        master.title('TurboWAVE Core Installer')
        self.pack()
        self.cmd = gui_command(self)
        self.conf = gui_config(self)
        self.tasks = gui_tasks(self)
        self.tasks.connect(self.cmd,self.conf)
        self.conf.connect(self.tasks)
        master.update_idletasks()
        ws = master.winfo_screenwidth()
        hs = master.winfo_screenheight()
        w = master.winfo_width()
        h = master.winfo_height()
        master.geometry('+%d+%d' % (ws/2-w/2,hs/2-h/2))
        tst = self.conf.data['Tools Version']
        if 'a' in tst or 'b' in tst or 'rc' in tst:
            confirm = self.cmd.affirm('Tools version is unstable, proceed anyway')
            if confirm=='n':
                self.master.destroy()
    def quit(self):
        if not self.tasks.finished():
            confirm = self.cmd.affirm('There are incomplete tasks, do you want to quit?')
            if confirm=='n':
                return
        self.master.destroy()

#####################
### TERMINAL MODE ###
#####################

class term_titlebar:
    def __init__(self,title,y=0,x=0):
        self.title = title
        self.window = curses.newwin(1,panel_width,y,x)
    def display(self):
        self.window.clear()
        self.window.addstr(0,0,self.title,curses.A_BOLD)
        self.window.refresh()

class term_popup:
    def __init__(self,items,y=1,x=1):
        self.title = 'Select'
        self.highlighted = 0
        lens = [len(s) for s in items]
        self.dims = (2+len(items),2+max(lens+[len(self.title)]))
        self.window = curses.newwin(self.dims[0],self.dims[1],y,x)
        self.choices = items
    def getyx(self):
        return self.dims
    def display(self):
        self.window.clear()
        self.window.border()
        self.window.addstr(0,int(self.window.getmaxyx()[1]/2-len(self.title)/2),self.title,curses.A_BOLD)
        for i,item in enumerate(self.choices):
            if i==self.highlighted:
                self.window.addstr(i+1,1,item,curses.A_REVERSE)
            else:
                self.window.addstr(i+1,1,item)
        self.window.keypad(True)
        self.window.refresh()
    def keypress(self,c):
        N = len(self.choices)
        selection = -1
        if c==ord('q') or c==ord('Q'):
            return 'escape'
        if c==27:
            self.window.nodelay(True)
            if self.window.getch()==-1:
                return 'escape'
            self.window.nodelay(False)
        if self.highlighted==-1:
            self.highlighted = 0
        else:
            if c==curses.KEY_LEFT or c==curses.KEY_UP:
                self.highlighted -= 1
                if self.highlighted<0:
                    self.highlighted = N-1
            if c==curses.KEY_RIGHT or c==curses.KEY_DOWN:
                self.highlighted += 1
                if self.highlighted>=N:
                    self.highlighted = 0
            if c==ord('\n') or c==ord('\r'):
                selection = self.highlighted
        self.display()
        if selection==-1:
            return ''
        else:
            return self.choices[selection]
    def run(self):
        self.display()
        sel = ''
        while sel=='':
            c = self.window.getch()
            sel = self.keypress(c)
        return sel

class term_dictionary_view(dictionary_view):
    def __init__(self):
        super().__init__()
    def display(self):
        self.window.clear()
        self.window.border()
        self.window.addstr(0,int(self.window.getmaxyx()[1]/2-len(self.title)/2),self.title,curses.A_BOLD)
        for i,key in enumerate(self.data):
            line = str(i+1)+') '+key+' = '+self.data[key]
            if i==self.highlighted:
                if key in self.trouble:
                    self.window.addstr(i+1,1,line,curses.color_pair(1) | curses.A_REVERSE)
                else:
                    self.window.addstr(i+1,1,line,curses.A_REVERSE)
            else:
                if key in self.trouble:
                    self.window.addstr(i+1,1,line,curses.color_pair(1))
                else:
                    self.window.addstr(i+1,1,line)
        self.window.keypad(True)
        self.window.refresh()
    def set(self,key,val):
        self.data[key] = val
    def cyclic(self,i,N,delta,indices):
        i += delta
        if i<0:
            i = N-1
        if i>=N:
            i = 0
        while (i not in indices):
            i += delta
            if i<0:
                i = N-1
            if i>=N:
                i = 0
        return i
    def keypress(self,c):
        selection = -1
        max_index = len(self.data)
        indices = []
        for i,key in enumerate(self.data):
            if key in self.selectable:
                indices += [i]
        if self.highlighted==-1:
            self.highlighted = indices[0]
        else:
            if c==curses.KEY_LEFT or c==curses.KEY_UP:
                self.highlighted = self.cyclic(self.highlighted,max_index,-1,indices)
            if c==curses.KEY_RIGHT or c==curses.KEY_DOWN:
                self.highlighted = self.cyclic(self.highlighted,max_index,1,indices)
            if c==ord('\n') or c==ord('\r'):
                selection = self.highlighted
        self.display()
        if selection==-1:
            return ''
        else:
            return list(self.data.keys())[selection]
    def popup(self,items,callback):
        obj = term_popup(items)
        choice = obj.run()
        if choice!='escape':
            callback(choice)
        return choice

class term_config(base_config,term_dictionary_view):
    def __init__(self,y):
        super().__init__()
        self.window = curses.newwin(len(self.data)+2,panel_width,y,0)
    def rows(self):
        return len(self.data)+2
    def display(self):
        self.verify()
        super().display()

class term_tasks(base_tasks,term_dictionary_view):
    def __init__(self,y):
        super().__init__()
        self.window = curses.newwin(len(self.data)+2,panel_width,y,0)
        self.highlighted = 0
    def rows(self):
        return len(self.data)+2
    def display(self):
        self.set_ok()
        for key in self.data:
            if self.data[key]=='incomplete':
                self.set_trouble(key)
        super().display()
    def enter_shell(self):
        curses.endwin()
    def exit_shell(self):
        newscr = curses.initscr()
        newscr.refresh()
        curses.doupdate()
    def AskDirectory(self,title='Select Working Directory',initialdir=str(pathlib.Path.home())):
        self.cmd.display(title+'\nMust give complete path (no tab completion).\nEnvironment variables and ~ OK to use.')
        curses.echo()
        ans = self.cmd.getstr(3,0)
        curses.noecho()
        ans = ans.replace('~',str(pathlib.Path.home()))
        pattern = r'\$(\w+)'
        match = re.search(pattern,ans)
        while match!=None:
            ans = re.sub(pattern,os.environ[match[1]],ans,count=1)
            match = re.search(pattern,ans)
        if len(glob.glob(ans))==0:
            self.cmd.err('Could not find path.')
            ans = ()
        else:
            self.cmd.display('')
        return ans

class term_command:
    def __init__(self,y):
        self.window = curses.newwin(5,panel_width,y,0)
        self.default_mess = 'arrows=highlight, enter=select, q=exit, t=toggle'
        self.last_mess = 'arrows=highlight, enter=select, q=exit, t=toggle'
    def connect(self,conf,tasks,titlebar):
        self.conf = conf
        self.tasks = tasks
        self.titlebar = titlebar
    def display(self,message):
        self.window.clear()
        if message=='':
            self.window.addstr(0,0,self.default_mess)
            self.last_mess = self.default_mess
        elif message=='persist':
            self.window.addstr(0,0,self.last_mess)
        else:
            self.window.addstr(0,0,message)
            self.last_mess = message
        self.window.keypad(True)
        self.window.refresh()
    def err(self,message):
        self.display(message+'\n'+self.default_mess)
    def info(self,title,message):
        self.display(message+' (press any key)')
        self.window.getch()
        self.display('')
    def affirm(self,message):
        self.display(message+' (y/n)?')
        confirm = ord('x')
        while confirm not in [ord('y'),ord('n')]:
            confirm = self.window.getch()
        self.display('')
        return chr(confirm)
    def getstr(self,y,x):
        curses.curs_set(1)
        s = str(self.window.getstr(y,x),encoding='utf-8')
        curses.curs_set(0)
        return s
    def run(self):
        key_focus = self.tasks
        tst = self.conf.data['Tools Version']
        if 'a' in tst or 'b' in tst or 'rc' in tst:
            confirm = self.affirm('Tools version is unstable, proceed anyway')
            if confirm=='n':
                return
        while True:
            self.titlebar.display()
            self.conf.display()
            self.tasks.display()
            self.display('persist')
            c = self.window.getch()
            if c==ord('q'):
                if self.tasks.finished():
                    break
                confirm = self.affirm('There are incomplete tasks, do you want to quit')
                if confirm=='y':
                    break
            if c==ord('t'):
                if key_focus==self.tasks:
                    key_focus = self.conf
                    self.conf.highlight(0)
                    self.tasks.highlight(-1)
                else:
                    key_focus = self.tasks
                    self.tasks.highlight(0)
                    self.conf.highlight(-1)
            sel = key_focus.keypress(c)
            if sel in self.conf.get_params():
                obj = term_popup(self.conf.get_choices(sel))
                res = obj.run()
                if res!='escape':
                    self.tasks.dirty_config()
                    self.conf.set(sel,res)
            if sel=='Get Components':
                self.tasks.GetComponents()
            if sel=='Set Version':
                self.tasks.StartSetVersion()
            if sel=='Configure Makefile':
                self.tasks.Configure()
            if sel=='Compile Code':
                self.tasks.Compile(self.tasks.enter_shell,self.tasks.exit_shell)
            if sel=='Install Files':
                self.tasks.Install()
            if sel=='Cleanup':
                self.tasks.Clean()

def terminal(stdscr):
    curses.init_pair(1,curses.COLOR_RED,curses.COLOR_BLACK)
    stdscr.clear()
    curses.curs_set(0) # make cursor invisible

    titlebar = term_titlebar('TurboWAVE Core Installer')
    conf = term_config(1)
    tasks = term_tasks(1+conf.rows())
    cmd = term_command(1+conf.rows()+tasks.rows())
    conf.connect(tasks)
    tasks.connect(cmd,conf)
    cmd.connect(conf,tasks,titlebar)
    cmd.run()

def main():
    if '--help' in sys.argv or '-h' in sys.argv:
        print('TurboWAVE Interactive Core Installer.')
        print('Usage: twinstall [--help] [--terminal]')
        print('Arguments: --help,-h :  displays this message')
        print('           --version,-v : display version number')
        print('           --terminal,-t : use textual interface')
        exit(0)

    if '--version' in sys.argv or '-v' in sys.argv:
        print('twinstall is provided by twutils, version '+pkg_resources.get_distribution('twutils').version)
        exit(0)

    if '--terminal' in sys.argv or '-t' in sys.argv:
        os.environ.setdefault('ESCDELAY','25')
        curses.wrapper(terminal)

    else:
        root = tk.Tk()
        ttk.Style().theme_use('default')
        app = Application(master=root)
        app.mainloop()
