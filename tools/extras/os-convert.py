import sys
import glob
import subprocess
import numpy as np
import h5py

if len(sys.argv)<2:
    print('Usage: os-convert.py component[,modifier] [mode=m,part] [frame=:] [root=MS]')
    print('Extracts data from OSIRIS directory tree and consolidates into npy file.')
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
    print('component = e1,e2,e3,b1,b2,b3,etc..')
    print('modifier = typically a species name, but could include other things.')
    print('m = 0,1,2... (only include for q3d).')
    print('part = re,im (only include for q3d).')
    print('frame = python style range indicating frames to keep, e.g., 2:11:4, etc.')
    print('root = OSIRIS root data directory relative to working directory')
    exit()

# Process command line

def range_str_to_tuple(rng_str):
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

id = sys.argv[1].split(',')
component = id[0]
if len(id)>1:
    modifier = id[-1]
else:
    modifier = ''

mode = 0
part = ''
frame_range = (0,-1,1)
root_dir = 'MS'
for keyval in sys.argv[2:]:
    key = keyval.split('=')[0]
    val = keyval.split('=')[1]
    if key=='frame':
        frame_range = range_str_to_tuple(val)
    if key=='root':
        root_dir = val
    if key=='mode':
        mode = int(val.split(',')[0])
        part = val.split(',')[1]

if len(part)>0:
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
    f = h5py.File(file_list[0])
except:
    print('Could not open '+file_list[0])
    exit(1)

print('Found HDF5 groups',list(f.keys()))
ax = f['AXIS']
dset = f[component]

# Set up the coordinate range along each axis.
# At present we do not account for changing coordinate range over time.
# Hence there will be a disparity in cases of moving window.
data_ext = []
data_ext += [ax['AXIS1'][0],ax['AXIS1'][1]]
try:
    data_ext += [ax['AXIS2'][0],ax['AXIS2'][1]]
    try:
        data_ext += [ax['AXIS3'][0],ax['AXIS3'][1]]
        data_ext = data_ext[4:6] + data_ext[2:4] + data_ext[:2]
    except:
        if 'PHA' in file_list[0]:
            data_ext = data_ext[2:4] + data_ext[:2] + [-0.5,0.5]
        else:
            data_ext = data_ext[2:4] + [-0.5,0.5] + data_ext[:2]
except:
    data_ext = [-0.5,0.5,-0.5,0.5] + data_ext

print('Data bounds',data_ext)
print('Dimensions',dset.shape)
print('Value range',np.min(dset),np.max(dset))

first = frame_range[0]
if frame_range[1]<0:
    last = num_frames + frame_range[1] + 1
else:
    last = frame_range[1]
step = frame_range[2]
total_frames = int((last-first)/step)

if len(part)>0:
    new_file = component+'-'+str(mode)+'-'+part+'.npy'
else:
    new_file = component+'.npy'
full_path_prefix = file_list[0][:-9]


for n in range(first,last,step):
    fullpath = full_path_prefix + '{:06d}.h5'.format(n)
    try:
        f = h5py.File(fullpath)
    except:
        print('Could not open HDF5 file.')
        exit(1)
    data_slice = np.array(f[component])
    if len(data_slice.shape)==1:
        data_slice = data_slice.reshape(1,1,data_slice.shape[0])
    if len(data_slice.shape)==2:
        if 'PHA' in file_list[0]:
            data_slice = data_slice.reshape(data_slice.shape[0],data_slice.shape[1],1)
        else:
            data_slice = data_slice.reshape(data_slice.shape[0],1,data_slice.shape[1])
    if n==frame_range[0]:
        print('Create',new_file)
        print('New shape',data_slice.shape)
        #dv.Create(new_file,data_slice,data_ext,True)
    else:
        print('Append',n)
        #dv.Append(new_file,data_slice,True)
