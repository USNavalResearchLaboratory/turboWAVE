import sys
import os
import glob
import json
import warnings
import numpy as np
import h5py
from scipy import constants as C

if len(sys.argv)<2:
    print('Usage: os-convert.py cmd component[,modifier] [n1=density] [dt=frame_time] [mode=m,part] [frame=:] [root=MS]')
    print('Extracts data from OSIRIS directory tree, consolidates into npy file, and creates metadata file.')
    print('1D coordinate label mapping: 1->z')
    print('2D coordinate label mapping: 1->z, 2->x')
    print('3D coordinate label mapping: 1->z, 2->y, 3->x')
    print('-------------Examples----------------')
    print('Cartesian field example: os-convert.py e1')
    print('Cartesian charge example: os-convert.py charge,electrons')
    print('Q3D field example: os-convert.py e3 mode=1,re')
    print('Q3D charge example: convert.py charge,electrons mode=0,re')
    print('Phase space example: os-convert.py p1x1,electrons')
    print('NOTE: you do NOT need to include _cyl_m')
    print('-------------Arguments----------------')
    print('cmd = init or append, if init restart the metadata, otherwise append to it.')
    print('component = e1,e2,e3,b1,b2,b3,etc..')
    print('modifier = typically a species name, but could include other things.')
    print('n1 = a normalizing density in cgs units, needed to support TW unit conversions.')
    print('dt = time between frames in plasma units.')
    print('m = 0,1,2... (only include for q3d).')
    print('part = re,im (only include for q3d).')
    print('frame = python style range indicating frames to keep, e.g., 2:11:4, etc.')
    print('root = OSIRIS root data directory relative to working directory')
    exit()

def range_str_to_tuple(rng_str):
    '''Convert a python style range string to a tuple'''
    rng = rng_str.split(':')
    tup = ()
    for i,el in enumerate(rng):
        if el=='' and i==0:
        	el = '0'
        if el=='' and i==1:
        	el = '-1'
        if el=='' and i==2:
        	el = '1'
        tup += (int(el),)
    return tup

def get_data_ext(ax):
    '''Convert OSIRIS bounds data, involving some re-arrangement and padding.

    :param dict ax: result of AXIS key in OSIRIS HDF5 file'''
    data_ext = []
    data_ext += [ax['AXIS1'][0],ax['AXIS1'][1]]
    try:
        data_ext += [ax['AXIS2'][0],ax['AXIS2'][1]]
        try:
            data_ext += [ax['AXIS3'][0],ax['AXIS3'][1]]
            data_ext = data_ext[4:6] + data_ext[2:4] + data_ext[:2]
        except:
            data_ext = data_ext[2:4] + [-0.5,0.5] + data_ext[:2]
    except:
        data_ext = [-0.5,0.5,-0.5,0.5] + data_ext
    return data_ext

# Process command line
# Required arguments
cmd = sys.argv[1]
if cmd not in ['init','append']:
    print('ERROR: not a valid command.')
    exit(1)
id = sys.argv[2].split(',')
component = id[0]
if len(id)>1:
    modifier = id[-1]
else:
    modifier = ''
# Defaults for optional arguments
mode = 0
part = ''
frame_range = (0,-1,1)
root_dir = 'MS'
n1_mks = 1e19
dt = 1.0
# Process optional arguments
for keyval in sys.argv[3:]:
    key = keyval.split('=')[0]
    val = keyval.split('=')[1]
    if key=='frame':
        frame_range = range_str_to_tuple(val)
    if key=='root':
        root_dir = val
    if key=='mode':
        mode = int(val.split(',')[0])
        part = val.split(',')[1]
    if key=='n1':
        n1_mks = 1e6*np.float(val)
    if key=='dt':
        dt = np.float(val)

# Construct OSIRIS filename prefix from the information
if len(part)>0 and component[-6:]!='_cyl_m':
    component = component + '_cyl_m'
prefix = ''
if len(modifier)>0:
    if part=='':
        prefix = component + '-' + modifier + '-'
    else:
        prefix = component + '-' + modifier + str(mode) + '-' + part + '-'
else:
    if part=='':
        prefix = component + '-'
    else:
        prefix = component + '-' + str(mode) + '-' + part + '-'

# See if we can find the specified file prefix
glob_str = root_dir + '/**/' + prefix + '*.h5'
print('Searching for',glob_str)
file_list = glob.glob(glob_str,recursive=True)
num_frames = len(file_list)
print('Found {:04d} h5 files with expected prefix.'.format(num_frames))
if num_frames==0:
    print('No files, so nothing to do.')
    exit(1)

# Get the axis and dimension information based on first frame

try:
    f = h5py.File(file_list[0],'r')
except:
    print('Could not open '+file_list[0])
    exit(1)

print('Found HDF5 groups',list(f.keys()))
ax = f['AXIS']
dset = f[component]

print('Dimensions',dset.shape)
print('Initial data bounds',get_data_ext(ax))
print('Initial value range',np.min(dset),np.max(dset))

# Set up the frames to grab
first = frame_range[0]
if frame_range[1]<0:
    last = num_frames + frame_range[1] + 1
else:
    last = frame_range[1]
step = frame_range[2]
total_frames = int((last-first)/step)

# name of the file to create
new_file = prefix[:-1]+'.npy'
# path of the file to read, not counting frame and extension
full_path_prefix = file_list[0][:-9]

def main_metadata(name,meta):
    meta[name] = {}
    meta[name]['native units'] = 'plasma'
    meta[name]['grid'] = name + '_grid_warp.txt'
    meta[name]['axes'] = {}
    return meta

def axis_metadata(ax,fullpath):
    path = os.path.split(fullpath)[0]
    file = os.path.split(fullpath)[1]
    wp = np.sqrt(n1_mks*C.e**2/C.epsilon_0/C.m_e)
    xa = (C.c/wp) / (C.hbar/C.m_e/C.c/C.alpha)
    ta = (1/wp) / (C.hbar/C.m_e/C.c**2/C.alpha**2)
    xn = (C.c/wp) / (C.hbar/C.m_e/C.c)
    tn = (1/wp) / (C.hbar/C.m_e/C.c**2)
    emks = C.m_e*C.c*wp/C.e
    ecgs = emks / (1e-4*C.c)
    ea = emks / (C.m_e**2*C.c**3*C.alpha**3 / C.hbar / C.e)
    en = emks / (C.m_e**2*C.c**3 / C.hbar / C.e)
    if 'FLD' not in path and 'PHA' not in path:
        print('ERROR: we only support converting fields or phase space.')
        exit(1)
    if 'FLD' in path:
        comp = file[0].upper()
        if comp!='E' and comp!='B':
            warnings.warn('Metadata is specialized for E or B')
        m = { 'E' : (emks,ecgs,1.0,ea,en) , 'B' : (emks/C.c,ecgs,1.0,ea,en) }
        u = { 'E' : ('V$/$m','SV$/$cm') , 'B' : ('T','G') }
        lab_key = ['mks','cgs','plasma','atomic','natural']
        label = ['$t$','$x$','$y$','$z$','$'+comp+'$']
        mult = [[1/wp]+[C.c/wp]*3+[m[comp][0]] , [1/wp]+[100*C.c/wp]*3+[m[comp][1]] , [1.0,1.0,1.0,1.0,1.0] , [ta,xa,xa,xa,m[comp][3]] , [tn,xn,xn,xn,m[comp][4]]]
        ulab = []
        ulab += [['s']+['m']*3+[u[comp][0]]]
        ulab += [['s']+['cm']*3+[u[comp][1]]]
        ulab += [[r'$\omega$']+[r'$\omega/c$']*3+[r'$e/mc\omega$']]
    if 'PHA' in path:
        warnings.warn('Metadata is specialized for p1x1')
        comp = file[0].upper()
        lab_key = ['mks','cgs','plasma','atomic','natural']
        label = ['$t$','$p_z$','$y$','$z$','$f$']
        mult = [[1/wp]+[C.m_e*C.c]+[C.c/wp]*2+[1.0] , [1/wp]+[1e5*C.m_e*C.c]+[100*C.c/wp]*2+[1.0] , [1.0,1.0,1.0,1.0,1.0] , [ta,0.0,xa,xa,1.0] , [tn,0.0,xn,xn,1.0]]
        ulab = []
        ulab += [['s']+['kg$\cdot$m$/$s']+['m']*2+['']]
        ulab += [['s']+['g$\cdot$cm$/$s']+['cm']*2+['']]
        ulab += [[r'$\omega$']+[r'$/mc$']+[r'$\omega/c$']*2 + ['']]
    ulab += [['a.u.']*5]
    ulab += [['n.u.']*5]
    a = {}
    a['label'] = label[ax]
    for i in range(len(lab_key)):
        a[lab_key[i]+' label'] = ulab[i][ax]
        a[lab_key[i]+' multiplier'] = mult[i][ax]
    return a

def grid_string(frame,bounds,shp):
    s = 't = ' + str(frame*dt) + '\n'
    pts = lambda x0,x1,n : np.linspace(x0+0.5*(x1-x0)/n , x1-0.5*(x1-x0)/n , n)
    for ax in range(1,4):
        s += 'axis' + str(ax) + ' ='
        for val in pts(bounds[int(2*ax-2)],bounds[int(2*ax-1)],shp[ax-1]):
            s += ' ' + str(val)
        s += '\n'
    return s

# TW shapes are always 4 dimensional, so add necessary axes
if len(dset.shape)==1:
	new_shape = (1,1) + dset.shape
if len(dset.shape)==2:
	new_shape = (dset.shape[0],1,dset.shape[1])
if len(dset.shape)==3:
	new_shape = dset.shape
full_array = np.zeros((total_frames,)+new_shape)
grid_str = ''
for n in range(first,last,step):
    fullpath = full_path_prefix + '{:06d}.h5'.format(n)
    try:
        f = h5py.File(fullpath,'r')
    except:
        print('Could not open HDF5 file.')
        exit(1)
    full_array[n,...] = np.array(f[component]).reshape(new_shape)
    grid_str += grid_string(n,get_data_ext(f['AXIS']),new_shape)
    print('.',end='',flush=True)

# Save the data into the .npy file
np.save(new_file,full_array)
print()

# create the metadata

if cmd=='init':
    meta = {}
else:
    with open('tw_metadata.json','r') as f:
        meta = json.load(f)
meta = main_metadata(new_file,meta)
for ax in range(5):
    meta[new_file]['axes'][str(ax)] = axis_metadata(ax,fullpath)
with open('tw_metadata.json','w') as f:
    json.dump(meta,f,indent=4)

with open(new_file+'_grid_warp.txt','w') as f:
    f.write(grid_str)
