import sys
import pathlib
import glob
import json
import warnings
import numpy as np
import h5py
from scipy import constants as C

def print_usage():
    print('Usage: os2tw <cmd> <component>[,<modifier>] [n1=<density>] [dt=<frame_time>] [frame=<::>] [root=<MS>]')
    print('Extracts data from OSIRIS directory tree, consolidates into npy file, and creates metadata file.')
    print('1D coordinate label mapping: 1->z')
    print('2D coordinate label mapping: 1->z, 2->x')
    print('3D coordinate label mapping: 1->z, 2->y, 3->x')
    print('-------------Examples----------------')
    print('Field example: os-convert.py init e1')
    print('Charge example: os-convert.py init charge,electrons')
    print('Phase space example: os-convert.py init p1x1,electrons')
    print('-------------Arguments----------------')
    print('cmd = init or append, if init restart the metadata, otherwise append to it.')
    print('component = e1,e2,e3,b1,b2,b3,etc..')
    print('modifier = used to narrow the file search, typically a species name, but could include other things.')
    print('n1 = a normalizing density in cgs units, needed to support TW unit conversions.')
    print('dt = time between frames in plasma units.')
    print('frame = python style range indicating frames to keep, e.g., 2:11:4, etc.')
    print('root = OSIRIS root data directory relative to working directory')

def str_to_range(rng_str,num):
    '''Convert a numpy-style slice string to a range'''
    rng = rng_str.split(':')
    tup = ()
    for i,el in enumerate(rng):
        if el=='' and i==0:
        	el = '0'
        if el=='' and i==1:
        	el = str(num)
        if el=='' and i==2:
        	el = '1'
        tup += (int(el),)
    return range(tup[0],tup[1],tup[2])

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

def main_metadata(name,meta):
    meta[name] = {}
    meta[name]['native units'] = 'plasma'
    meta[name]['grid'] = name + '_grid_warp.txt'
    meta[name]['axes'] = {}
    return meta

def axis_metadata(ax,component):
    wp = np.sqrt(n1_mks*C.e**2/C.epsilon_0/C.m_e)
    xa = (C.c/wp) / (C.hbar/C.m_e/C.c/C.alpha)
    ta = (1/wp) / (C.hbar/C.m_e/C.c**2/C.alpha**2)
    xn = (C.c/wp) / (C.hbar/C.m_e/C.c)
    tn = (1/wp) / (C.hbar/C.m_e/C.c**2)
    emks = C.m_e*C.c*wp/C.e
    ecgs = emks / (1e-4*C.c)
    ea = emks / (C.m_e**2*C.c**3*C.alpha**3 / C.hbar / C.e)
    en = emks / (C.m_e**2*C.c**3 / C.hbar / C.e)
    units = ['mks','cgs','plasma','atomic','natural']
    label = ['$t$','$x$','$y$','$z$','$'+component+'$']
    fields = ['e1','e2','e3','b1','b2','b3']
    phase_space_axes = ['x1','x2','x3','p1','p2','p3']
    phase_spaces = []
    for i in range(6):
        for j in range(6):
            if i!=j:
                phase_spaces += [phase_space_axes[i]+phase_space_axes[j]]
    l = {}
    l['none'] = {'mks':'AU' , 'cgs':'AU' , 'plasma':'', 'atomic':'AU', 'natural':'AU'}
    l['t'] = {'mks':'s' , 'cgs':'s' , 'plasma':r'$\omega$', 'atomic':'a.u.', 'natural':'n.u.'}
    l['x'] = {'mks':'m' , 'cgs':'cm' , 'plasma':r'$\omega/c$', 'atomic':'a.u.', 'natural':'n.u.'}
    l['p'] = {'mks':'kg$\cdot$m$/$s' , 'cgs':'g$\cdot$cm$/$s' , 'plasma':r'$/mc$', 'atomic':'a.u.', 'natural':'n.u.'}
    l['e'] = {'mks':'V$/$m' , 'cgs':'SV$/$cm' , 'plasma':r'$e/mc\omega$', 'atomic':'a.u.', 'natural':'n.u.'}
    l['b'] = {'mks':'T' , 'cgs':'G' , 'plasma':r'$e/mc\omega$', 'atomic':'a.u.', 'natural':'n.u.'}
    m = {}
    m['none'] = {'mks':1.0 , 'cgs':1.0 , 'plasma':1.0, 'atomic':1.0, 'natural':1.0}
    m['t'] = {'mks':1/wp , 'cgs':1/wp , 'plasma':1.0, 'atomic':xa, 'natural':xn}
    m['x'] = {'mks':C.c/wp , 'cgs':100*C.c/wp , 'plasma':1.0, 'atomic':xa, 'natural':xn}
    m['p'] = {'mks':C.m_e*C.c , 'cgs':1e5*C.m_e*C.c , 'plasma':1.0, 'atomic':C.m_e*C.c*C.alpha, 'natural':C.m_e*C.c}
    m['e'] = {'mks':emks , 'cgs':ecgs , 'plasma':1.0, 'atomic':ea, 'natural':en}
    m['b'] = {'mks':emks/C.c , 'cgs':ecgs , 'plasma':1.0, 'atomic':ea, 'natural':en}
    mult = [[]]*5
    ulab = [[]]*5
    if component not in fields+phase_spaces:
        warnings.warn('Requested component not known, units will be arbitrary.')
        for i,s in enumerate(units):
            mult[i] = [m['none'][s]]*5
            ulab[i] = [l['none'][s]]*5
    if component in fields:
        label[4] = '$'+component[0].upper()+'_'+component[1]+'$'
        for i,s in enumerate(units):
            mult[i] = [m['t'][s]] + [m['x'][s]]*3 + [m[component[0]][s]]
            ulab[i] = [l['t'][s]] + [l['x'][s]]*3 + [l[component[0]][s]]
    if component in phase_spaces:
        h = component[0]+'_'+component[1]
        v = component[2]+'_'+component[3]
        label[1] = '$'+h+'$'
        label[3] = '$'+v+'$'
        label[4] = '$f('+h+','+v+')$'
        for i,s in enumerate(units):
            mult[i] = [1.0,m[component[0]][s],1.0,m[component[2]][s],1.0]
            ulab[i] = [l['none'][s],l[component[0]][s],l['none'][s],l[component[2]][s],l['none'][s]]
    a = {}
    a['label'] = label[ax]
    for i in range(len(units)):
        a[units[i]+' label'] = ulab[i][ax]
        a[units[i]+' multiplier'] = mult[i][ax]
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

def main():
    # Process command line
    # Required arguments
    if len(sys.argv)<2:
        print_usage()
        exit(0)
    cmd = sys.argv[1]
    if cmd not in ['init','append']:
        raise ValueError('First argument <'+cmd+'> not recognized.')
        print('ERROR: not a valid command.')
        exit(1)
    id = sys.argv[2].split(',')
    component = id[0]
    if len(id)>1:
        modifier = id[-1]
    else:
        modifier = ''
    # Defaults for optional arguments
    range_str = '::'
    root_dir = 'MS'
    n1_mks = 1e19
    dt = 1.0
    keys_found = []
    # Process optional arguments
    for keyval in sys.argv[3:]:
        key = keyval.split('=')[0]
        val = keyval.split('=')[1]
        keys_found += [key]
        if key not in ['frame','root','n1','dt']:
            raise KeyError('Command line argument <'+key+'> not recognized.')
        if key=='frame':
            range_str = val
        if key=='root':
            root_dir = val
        if key=='n1':
            n1_mks = 1e6*np.float(val)
        if key=='dt':
            dt = np.float(val)
    if 'n1' not in keys_found:
        warnings.warn('No density given, physical units will be arbitrary.')
    if 'dt' not in keys_found:
        warnings.warn('No time step given, time axis will be arbitrary.')

    # Search for the desired OSIRIS data file
    if modifier=='':
        glob_str = root_dir + '/**/' + component + '*[0-9][0-9][0-9][0-9][0-9][0-9].h5'
    else:
        glob_str = root_dir + '/**/' + component + '*' + modifier + '*[0-9][0-9][0-9][0-9][0-9][0-9].h5'
    file_list = glob.glob(glob_str,recursive=True)
    num_frames = len(file_list)
    if num_frames==0:
        print('ERROR: no files matched the requested pattern.')
        exit(1)
    prefix_list = list(set([s[:-9] for s in file_list]))
    if len(prefix_list)>1:
        print('ERROR: found more than one match, use modifier to narrow search.')
        for i,res in enumerate(prefix_list):
            print('    Match '+str(i)+': '+res)
        exit(1)
    full_path_prefix = prefix_list[0]
    prefix = pathlib.Path(full_path_prefix).name
    print('Found {:04d} h5 files with prefix {}.'.format(num_frames,prefix))

    # Set up the frames to grab
    frame_list = []
    for frame in str_to_range(range_str,num_frames):
        frame_list += [frame]

    # Get the axis and dimension information based on first frame

    try:
        f = h5py.File(file_list[0],'r')
    except:
        print('Could not open '+file_list[0])
        exit(1)

    print('Found HDF5 groups',list(f.keys()))
    ax = f['AXIS']
    data_key = None
    for k in f.keys():
        if component in k:
            data_key = k
    if data_key==None:
        print('ERROR: could not match to an HDF5 key.')
        exit(1)

    dset = f[data_key]
    print('Dimensions',dset.shape)
    print('Initial data bounds',get_data_ext(ax))
    print('Initial value range',np.min(dset),np.max(dset))

    # name of the file to create
    new_file = prefix[:-1]+'.npy'

    # TW shapes are always 4 dimensional, so add necessary axes
    if len(dset.shape)==1:
    	new_shape = (1,1) + dset.shape
    if len(dset.shape)==2:
    	new_shape = (dset.shape[0],1,dset.shape[1])
    if len(dset.shape)==3:
    	new_shape = dset.shape
    full_array = np.zeros((len(frame_list),)+new_shape)
    grid_str = ''
    count = 0
    for n in frame_list:
        fullpath = full_path_prefix + '{:06d}.h5'.format(n)
        try:
            f = h5py.File(fullpath,'r')
        except:
            print('Could not open HDF5 file <' + fullpath + '>.')
            exit(1)
        full_array[count,...] = np.array(f[data_key]).reshape(new_shape)
        grid_str += grid_string(n,get_data_ext(f['AXIS']),new_shape)
        count += 1
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
        meta[new_file]['axes'][str(ax)] = axis_metadata(ax,component)
    with open('tw_metadata.json','w') as f:
        json.dump(meta,f,indent=4)

    with open(new_file+'_grid_warp.txt','w') as f:
        f.write(grid_str)
