"""Module to setup plots based on npy and accompanying metadata.
Primary interaction is through the plotter class."""

import sys
import os
import json
import warnings
import importlib
import numpy as np
from scipy import constants as C

def axis_num(axis_letter):
    if axis_letter=='t':
        return 0
    if axis_letter=='x':
        return 1
    if axis_letter=='y':
        return 2
    if axis_letter=='z':
        return 3

def axis_lett(axis_num):
    if axis_num==0:
        return 't'
    if axis_num==1:
        return 'x'
    if axis_num==2:
        return 'y'
    if axis_num==3:
        return 'z'

def expand_slice_spec(primitive_spec):
    """Given a 2-character plane of interest add the letters for the remaining 2 axes.
    The order is such that time will never be last."""
    h = axis_num(primitive_spec[0])
    v = axis_num(primitive_spec[1])
    s1 = 3
    while s1==h or s1==v:
        s1 = s1 - 1
    s2 = 3
    while s2==h or s2==v or s2==s1:
        s2 = s2 - 1
    return primitive_spec + axis_lett(s2) + axis_lett(s1)

def get_axis_info(slicing_spec):
    """Get an ordered array of axis indices from a slicing string.
    Argument is a 4 character string giving the axis order (field component is always last).
    Return a 5 element axis array of integers, with encoding t=0,x=1,y=2,z=3,component=4.
    Example: 'zxyt' would give [3,1,2,0,4]."""
    axes = np.array([0,0,0,0,0]).astype(np.int)
    for i in range(4):
        axes[i] = axis_num(slicing_spec[i])
    axes[4] = 4
    return axes

def get_mesh_pts(grid_file_path,dims):
    """Try to find a grid warp file matching the npy file.
    If the file is found return [tmesh,xmesh,ymesh,zmesh].
    Each element is a list of the mesh points along the given axis.
    If the file is not found return unit length uniform mesh."""
    try:
        ans = [[],[],[],[]]
        with open(grid_file_path) as f:
            for line in f:
                l = line.split(' ')
                if l[0]=='t':
                    ans[0] += [np.float(l[-1])]
                if l[0]=='axis1' and len(ans[1])==0:
                    ans[1] = [np.float(x) for x in l[2:]]
                if l[0]=='axis2' and len(ans[2])==0:
                    ans[2] = [np.float(x) for x in l[2:]]
                if l[0]=='axis3' and len(ans[3])==0:
                    ans[3] = [np.float(x) for x in l[2:]]
        if len(ans[0])==dims[0] and len(ans[1])==dims[1] and len(ans[2])==dims[2] and len(ans[3])==dims[3]:
            return [np.array(ans[0]),np.array(ans[1]),np.array(ans[2]),np.array(ans[3])]
        else:
            warnings.warn('Grid file found but wrong dimensions ('+grid_file_path+').')
            return [np.linspace(0.0,1.0,dims[i]) for i in range(4)]
    except:
        warnings.warn('No grid file found ('+grid_file_path+').')
        return [np.linspace(0.0,1.0,dims[i]) for i in range(4)]

class plotter:
    """Class that retains information about turboWAVE data,
    and generates objects suitable for direct use in matplotlib."""
    def __init__(self,file_to_plot,units='plasma',buffered=True,small_pos=1e-25):
        self.basepath = os.path.split(file_to_plot)[0]
        self.name = os.path.split(file_to_plot)[1] # filename only
        self.path = file_to_plot # includes path
        self.units = units
        self.buffered = buffered
        try:
            with open(os.path.join(self.basepath,'tw_metadata.json'),'r') as f:
                meta_str = f.read()
        except:
            warnings.warn('Could not read file '+os.path.join(self.basepath,'tw_metadata.json'))
        try:
            self.meta = json.loads(meta_str)
        except:
            warnings.warn('Could not decode tw_metadata.json.')
        try:
            self.axes = self.meta[self.name]['axes']
            for i in range(5):
                self.axes[i] = self.axes[str(i)]
        except:
            warnings.warn('Could not retrieve axes from the metadata for '+self.name+', using arbitrary scales.')
            self.axes = { 0 : {} , 1 : {} , 2 : {} , 3 : {} , 4 : {} }
            self.axes[0]['label'] = 't'
            self.axes[1]['label'] = 'x'
            self.axes[2]['label'] = 'y'
            self.axes[3]['label'] = 'z'
            self.axes[4]['label'] = 'data'
            for i in range(5):
                self.axes[i]['mks label'] = 'None'
                self.axes[i]['mks multiplier'] = 1.0
                self.axes[i]['cgs label'] = 'None'
                self.axes[i]['cgs multiplier'] = 1.0
                self.axes[i]['plasma label'] = 'None'
                self.axes[i]['plasma multiplier'] = 1.0
                self.axes[i]['atomic label'] = 'None'
                self.axes[i]['atomic multiplier'] = 1.0
                self.axes[i]['natural label'] = 'None'
                self.axes[i]['natural multiplier'] = 1.0
        try:
            grid_file_path = os.path.join(self.basepath,self.meta[self.name]['grid'])
        except:
            warnings.warn('Could not find a reference to a grid file in the metadata.')
            grid_file_path = ''
        if buffered:
            self.data = np.load(self.path)
        else:
            self.data = np.load(self.path,mmap_mode='c')
        self.small_pos = small_pos
        self.mesh = get_mesh_pts(grid_file_path,self.data.shape)
        self.data *= self.axes[4][self.units+' multiplier']
        for a,mesh_pts in enumerate(self.mesh):
            mesh_pts *= self.axes[a][self.units+' multiplier']

    def name(self):
        return self.name
    def get(self):
        return self.data
    def set(self,new_data):
        self.data = new_data
    def scale(self,factor):
        self.data *= factor
    def max_frame(self):
        return self.data.shape[0]-1
    def phys_dims4(self):
        if self.mesh==None:
            span = [0.0,1.0]*4
        else:
            span = []
            for ax in range(4):
                m = self.mesh[ax]
                if len(m)<2:
                    span += [0.0,1.0]
                else:
                    span += [1.5*m[0]-0.5*m[1] , 1.5*m[-1]-0.5*m[-2]]
        return np.array(span)
    def dims4(self):
        return self.data.shape
    def display_info(self):
        print('-------------------------')
        print('grid cells =',self.data.shape)
        print('physical size =',self.phys_dims4())
        print('-------------------------')
    def plot_coord(self,dict,ax,idx):
        # N.b. ax is relative, i.e., 0=hor., 1=ver.
        ds = (dict['extent'][2*ax+1] - dict['extent'][2*ax])/dict['dims'][ax]
        return dict['extent'][2*ax] + ds*idx + ds/2
    def get_global_minmax(self,dyn_range):
        if dyn_range==0.0:
            return np.min(self.data),np.max(self.data)
        else:
            logmax = np.log10(np.max(np.abs(self.data))+self.small_pos)
            return logmax-dyn_range , logmax
    def update_dyn_range(self,plot_dict,data_slice,dyn_range):
        lab = plot_dict['blabel']
        if dyn_range==0.0:
            plot_dict['vmin'] = np.min(data_slice)
            plot_dict['vmax'] = np.max(data_slice)
            if len(lab)>11:
                if lab[:11]=='log$_{10}$ ':
                    plot_dict['blabel'] = lab[11:]
        else:
            maxval = np.max(np.abs(data_slice))
            logmax = np.log10(maxval+self.small_pos)
            plot_dict['vmin'] = logmax-dyn_range
            plot_dict['vmax'] = logmax
            data_slice[...] = np.log10(np.abs(data_slice)+self.small_pos)
            if len(lab)>11:
                if lab[:11]!='log$_{10}$ ':
                    plot_dict['blabel'] = 'log$_{10}$ ' + lab
            else:
                plot_dict['blabel'] = 'log$_{10}$ ' + lab

    def start_plot(self,slicing_spec='zxty'):
        '''Not directly called.
        Returns an initial plotting dictionary for use with any kind of plot'''
        if len(slicing_spec)!=4:
            raise TypeError('The slicing_spec was not 4 characters.')
        axes = get_axis_info(slicing_spec)
        phys_dims = self.phys_dims4()
        plot_dict = {}

        lab = [self.axes[axes[i]]['label'] for i in range(5)]
        ulab = [self.axes[axes[i]][self.units+' label'] for i in range(5)]
        dkey = ['xlabel','ylabel','zlabel','none','blabel']
        if self.units=='plasma':
            for i in [0,1,2,4]:
                if ulab[i]=='None':
                    plot_dict[dkey[i]] = lab[i]
                else:
                    s = ulab[i].find('/')
                    if s==-1:
                        plot_dict[dkey[i]] = ulab[i] + lab[i]
                    else:
                        plot_dict[dkey[i]] = (ulab[i][:s] + ' $ ' + lab[i] + ' $ ' + ulab[i][s:]).replace('$ $','')
        else:
            for i in [0,1,2,4]:
                if ulab[i]=='None':
                    plot_dict[dkey[i]] = lab[i]
                else:
                    plot_dict[dkey[i]] = lab[i] + ' (' + ulab[i] + ')'
        plot_dict['dims'] = [self.dims4()[axes[i]] for i in range(3)]
        plot_dict['xpts'] = self.mesh[axes[0]]
        plot_dict['ypts'] = self.mesh[axes[1]]
        plot_dict['zpts'] = self.mesh[axes[2]]
        return plot_dict

    def volume3d(self,slicing_spec='xyzt',slice_to_plot=0,dyn_range=0.0):
        '''Returns data,dict
        data = the numpy array with the 3D slice data
        dict = a plotting dictionary to be used externally with mayavi
        The dictionary keys often have the same names as mayavi function arguments.  Available keys:
        "extent", "xlabel", "ylabel", "zlabel", "blabel", "dims", "vmin", "vmax", "mesh"
        Inputs:
        slicing_spec is a 4 character code ordering the axes, e.g., zxty
        The first 3 characters are the plotting volume, the last is the slicing axis
        slice_to_plot is an integer selecting a point on the slicing axis'''
        if len(slicing_spec)!=4:
            raise TypeError('The slicing_spec was not 4 characters.')
        if not slicing_spec=='xyzt':
            raise ValueError('The plotting axes have to be in order at present.')
        axes = get_axis_info(slicing_spec)
        phys_dims = self.phys_dims4()
        plot_dict = self.start_plot(slicing_spec)
        data_slice = np.take(self.data,[slice_to_plot],axis=axes[3])
        data_slice = np.squeeze(data_slice,axis=(axes[3],))

        ext = [[phys_dims[axes[i]*2],phys_dims[axes[i]*2+1]] for i in range(3)]
        plot_dict['extent'] = ext[0]+ext[1]+ext[2]
        self.update_dyn_range(plot_dict,data_slice,dyn_range)
        return data_slice,plot_dict

    def falsecolor2d(self,slicing_spec='zxty',slice_to_plot=(0,0),dyn_range=0.0,lib='matplotlib'):
        '''Returns data,dict
        data = the numpy array with the 2D slice data
        dict = a plotting dictionary to be used externally with MPL
        The dictionary keys have the same names as MPL function arguments.  Available keys:
        "extent", "xlabel", "ylabel", "dims", "aspect", "vmin", "vmax"
        Inputs:
        slicing_spec is a 4 character code ordering the axes, e.g., zxty
        The first two characters are the plotting plane, the next two are the slicing plane
        slice_to_plot is a tuple (s3,s4) selecting a point in the slicing plane'''
        if len(slice_to_plot)!=2:
            raise TypeError('The slice for falsecolor2d is supposed to be a 2-tuple.')
        #my_aspect = 'equal'
        my_aspect = 'auto'
        axes = get_axis_info(slicing_spec)
        phys_dims = self.phys_dims4()
        plot_dict = self.start_plot(slicing_spec)
        data_slice = np.take(self.data,[slice_to_plot[0]],axis=axes[2])
        data_slice = np.take(data_slice,[slice_to_plot[1]],axis=axes[3])
        data_slice = np.squeeze(data_slice,axis=(axes[2],axes[3]))

        # matplotlib plot data is accessed as [v,h]
        # mayavi plot data is accessed as [h,v]
        # twutils is supposed to put horizontal as the axis listed first
        if axes[0]<axes[1] and lib=='matplotlib':
            data_slice = data_slice.swapaxes(0,1)
        if axes[1]<axes[0] and lib=='mayavi':
            data_slice = data_slice.swapaxes(0,1)

        ext = [[phys_dims[axes[i]*2],phys_dims[axes[i]*2+1]] for i in range(2)]
        plot_dict['extent'] = ext[0]+ext[1]
        plot_dict['aspect'] = my_aspect
        self.update_dyn_range(plot_dict,data_slice,dyn_range)
        return data_slice,plot_dict

    def lineout(self,slicing_spec='xyzt',slice_to_plot=(0,0,0),dyn_range=0.0):
        '''Returns abcissa,ordinate,dict
        abcissa,ordinate = the numpy arrays with the 1D axis values and data
        dict = a plotting dictionary to be used externally with MPL
        The dictionary keys have the same names as MPL function arguments. Available keys:
        "extent", "xlabel", "ylabel", "dims"
        Inputs:
        slicing_spec is a 4 character code ordering the axes, e.g., zxty
        The first character is the plotting axis, the next three are the slicing volume
        slice_to_plot is a tuple (s2,s3,s4) selecting a point in the slicing volume'''
        if len(slice_to_plot)!=3:
            raise TypeError('The slice for lineout is supposed to be a 3-tuple.')
        axes = get_axis_info(slicing_spec)
        phys_dims = self.phys_dims4()
        plot_dict = self.start_plot(slicing_spec)
        data_slice = np.take(self.data,[slice_to_plot[0]],axis=axes[1])
        data_slice = np.take(data_slice,[slice_to_plot[1]],axis=axes[2])
        data_slice = np.take(data_slice,[slice_to_plot[2]],axis=axes[3])
        data_slice = np.squeeze(data_slice,axis=(axes[1],axes[2],axes[3]))
        if self.mesh!=None:
            abcissa = self.mesh[axes[0]]
        else:
            abcissa = np.linspace(0.0,1.0,len(data_slice))
        self.update_dyn_range(plot_dict,data_slice,dyn_range)
        plot_dict['ylabel'] = plot_dict['blabel']
        return abcissa,data_slice,plot_dict
