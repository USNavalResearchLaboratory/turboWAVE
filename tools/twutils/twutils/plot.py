"""Module to setup plots based on dvdat files.
Primary interaction is through the plotter class."""

import numpy as np
from scipy import constants as C
import twutils.dvdat as dv
import twutils.pre as twpre
import glob

def axis_num(axis_letter):
    # DVDAT ordering is t,z,y,x (C) or x,y,z,t (FORTRAN).
    # While reading into numpy array this is transformed to array[t,x,y,z]
    if axis_letter=='t':
        return 0
    if axis_letter=='x':
        return 1
    if axis_letter=='y':
        return 2
    if axis_letter=='z':
        return 3

def axis_lett(axis_num):
    # DVDAT ordering is t,z,y,x (C) or x,y,z,t (FORTRAN).
    # While reading into numpy array this is transformed to array[t,x,y,z]
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
    """Argument is a 4 character string giving the axis order.
    Return a 4 element axis array of integers, with encoding t=0,x=1,y=2,z=3"""
    axes = np.array([0,0,0,0]).astype(np.int)
    for i in range(4):
        axes[i] = axis_num(slicing_spec[i])
    return axes

def get_mesh_pts(data_file_name,dims):
    """Try to find a grid warp file matching the dvdat file.
    If the file is found return [tmesh,xmesh,ymesh,zmesh].
    Each element is a list of the mesh points along the given axis.
    If the file is not found return None."""
    if data_file_name[-6:]==".dvdat":
        basename = data_file_name[:-6]
    else:
        basename = data_file_name
    frags = basename.split('_')
    # Strategy is to try matching the most complicated names first
    for i in range(len(frags),-1,-1):
        test = ''
        for s in frags[:i]:
            test += s + '_'
        res = glob.glob(test+'grid_warp.txt')
        if len(res)>0:
            break
    else:
        return None
    ans = [[],[],[],[]]
    with open(res[0]) as f:
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
        return None

class plotter:
    """Class that retains information about a DVDAT file,
    and generates objects suitable for direct use in matplotlib.
    Will also try to get grid warp and time data."""
    def __init__(self,file_to_plot,buffered=True,small_pos=1e-25):
        # Following 2 are redundant for now, but may not be later
        self.name = file_to_plot
        self.file_to_plot = file_to_plot
        self.params = dv.GetParameters(file_to_plot)
        self.buffered = buffered
        if buffered:
            self.data = dv.GetFullArray(file_to_plot,c_order=True)
        else:
            self.data = dv.GetFrameArray(file_to_plot,self.params.frames-1,c_order=True)
        self.vmin_global = np.min(self.data)
        self.vmax_global = np.max(self.data)
        self.vmax_global_log10 = np.log10(np.max(np.abs(self.data)))
        self.small_pos = small_pos
        self.mesh = get_mesh_pts(file_to_plot,[self.params.frames,self.params.dims[0],self.params.dims[1],self.params.dims[2]])

    def name(self):
        return self.name
    def get(self):
        return self.data
    def set(self,new_data):
        self.data = new_data
    def scale(self,factor):
        self.data *= factor
    def max_frame(self):
        return self.params.frames-1
    def phys_dims4(self):
        if self.mesh==None:
            span = [0.0,1.0] + list(self.params.phys_dims)
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
        return np.concatenate(([self.params.frames],self.params.dims))
    def display_info(self):
        dv.PrintParameters(self.params)
        if self.mesh==None:
            print('Note: compatible warp file could not be found.')
    def plot_coord(self,dict,ax,idx):
        # N.b. ax is relative, i.e., 0=hor., 1=ver.
        ds = (dict['extent'][2*ax+1] - dict['extent'][2*ax])/dict['dims'][ax]
        return dict['extent'][2*ax] + ds*idx + ds/2
    def get_global_minmax(self,dyn_range):
        if dyn_range==0.0:
            return self.vmin_global,self.vmax_global
        else:
            return self.vmax_global_log10-dyn_range , self.vmax_global_log10

    def start_plot(self,slicing_spec='zxty',slice_to_plot=(0,0),lib='matplotlib'):
        '''Not directly called.
        Returns data,dict
        data = the numpy array with the 2D slice data
        dict = an initial plotting dictionary for use with any kind of plot'''
        if len(slicing_spec)!=4:
            raise TypeError('The slicing_spec was not 4 characters.')
        if len(slice_to_plot)!=2:
            raise TypeError('From within start_plot slice is supposed to be 2-tuple.')
        phys_dims = self.phys_dims4()
        axes = get_axis_info(slicing_spec)

        if not self.buffered:
            tpos = slicing_spec.find('t')
            if tpos<2:
                raise ValueError('Time cannot be a plotting axis unless the plotter is buffered.')
            spos = 5-tpos # evalutes to 2 or 3 if tpos does
            frame = slice_to_plot[tpos-2]
            slc = slice_to_plot[spos-2]
            self.data = dv.GetFrameArray(self.file_to_plot,frame,c_order=True)
            data_slice = np.take(self.data,[slc],axis=axes[spos]-1)
            data_slice = np.squeeze(data_slice,axis=axes[spos]-1)
        else:
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

        h0 = phys_dims[axes[0]*2]
        h1 = phys_dims[axes[0]*2+1]
        v0 = phys_dims[axes[1]*2]
        v1 = phys_dims[axes[1]*2+1]

        plot_dict = { 'extent' : [h0,h1,v0,v1] }
        plot_dict['xlabel'] = 'Axis'+str(axes[0])
        plot_dict['ylabel'] = 'Axis'+str(axes[1])
        plot_dict['dims'] = [self.dims4()[axes[0]],self.dims4()[axes[1]]]
        return data_slice,plot_dict

    def volume3d(self,slicing_spec='xyzt',slice_to_plot=0,dyn_range=0.0):
        '''Returns data,dict
        data = the numpy array with the 3D slice data
        dict = a plotting dictionary to be used externally with mayavi
        The dictionary keys have the same names as mayavi function arguments.  Available keys:
        "extent", "xlabel", "ylabel", "zlabel", "dims", "vmin", "vmax"
        Inputs:
        slicing_spec is a 4 character code ordering the axes, e.g., zxty
        The first 3 characters are the plotting volume, the last is the slicing axis
        slice_to_plot is an integer selecting a point on the slicing axis'''
        if len(slicing_spec)!=4:
            raise TypeError('The slicing_spec was not 4 characters.')
        phys_dims = self.phys_dims4()
        axes = get_axis_info(slicing_spec)
        if axes[0]>axes[1] or axes[0]>axes[2] or axes[1]>axes[2]:
            raise ValueError('The plotting axes have to be in order at present.')
        if not self.buffered:
            tpos = slicing_spec.find('t')
            if tpos<3:
                raise ValueError('Time cannot be a plotting axis unless the plotter is buffered.')
            self.data = dv.GetFrameArray(self.file_to_plot,slice_to_plot,c_order=True)
            data_slice = self.data
        else:
            data_slice = np.take(self.data,[slice_to_plot],axis=axes[3])
            data_slice = np.squeeze(data_slice,axis=(axes[3],))

        x0 = phys_dims[axes[0]*2]
        x1 = phys_dims[axes[0]*2+1]
        y0 = phys_dims[axes[1]*2]
        y1 = phys_dims[axes[1]*2+1]
        z0 = phys_dims[axes[2]*2]
        z1 = phys_dims[axes[2]*2+1]

        plot_dict = { 'extent' : [x0,x1,y0,y1,z0,z1] }
        plot_dict['xlabel'] = 'Axis'+str(axes[0])
        plot_dict['ylabel'] = 'Axis'+str(axes[1])
        plot_dict['zlabel'] = 'Axis'+str(axes[2])
        plot_dict['dims'] = [self.dims4()[axes[0]],self.dims4()[axes[1]],self.dims4()[axes[2]]]
        if dyn_range==0.0:
            plot_dict['vmin'] = np.min(data_slice)
            plot_dict['vmax'] = np.max(data_slice)
        else:
            maxval = np.max(np.abs(data_slice))
            logmax = np.log10(maxval+self.small_pos)
            plot_dict['vmin'] = logmax-dyn_range
            plot_dict['vmax'] = logmax
            data_slice = np.log10(np.abs(data_slice)+self.small_pos)
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
        #my_aspect = 'equal'
        my_aspect = 'auto'
        data_slice,plot_dict = self.start_plot(slicing_spec,slice_to_plot,lib)
        plot_dict['aspect'] = my_aspect
        if dyn_range==0.0:
            plot_dict['vmin'] = np.min(data_slice)
            plot_dict['vmax'] = np.max(data_slice)
        else:
            maxval = np.max(np.abs(data_slice))
            logmax = np.log10(maxval+self.small_pos)
            plot_dict['vmin'] = logmax-dyn_range
            plot_dict['vmax'] = logmax
            data_slice = np.log10(np.abs(data_slice)+self.small_pos)
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
        data_slice,plot_dict = self.start_plot(slicing_spec,slice_to_plot[1:])
        axes = get_axis_info(slicing_spec)
        # At this point we have 2 axes internally ordered as v,h, we are plotting h
        data_slice = np.take(data_slice,[slice_to_plot[0]],axis=0)
        data_slice = np.squeeze(data_slice,axis=(0,))
        if self.mesh!=None:
            abcissa = self.mesh[axes[0]]
        else:
            abcissa = np.linspace(plot_dict['extent'][0],plot_dict['extent'][1],len(data_slice))
        plot_dict['ylabel'] = 'Value'
        return abcissa,data_slice,plot_dict
