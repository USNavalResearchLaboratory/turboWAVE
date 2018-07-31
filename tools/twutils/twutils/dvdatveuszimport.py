"""
Veusz Plugin for turboWAVE dvdat files
"""
from veusz.plugins import *
import numpy as np

## really should import this from dvdat.py...

# define the data types occuring in a dvdat file
dt1 = np.dtype(np.int32)
dt1 = dt1.newbyteorder('>')
dt2 = np.dtype(np.float32)
dt2 = dt2.newbyteorder('>')

class dvdatParams:
	"""Contains parameters describing a dvdat file"""
	header = 'DataViewer 2.0.0'
	file_size = 0
	frames = 0
	floats_per_frame = 0
	dims = np.array([1,1,1],dtype=dt1)
	phys_dims = np.array([0.,1.,0.,1.,0.,1.],dtype=dt2)
	
def GetParameters(file_name):
	"""Returns the parameters of the dvdat file
	inputs:
	file_name (string)"""
	p = dvdatParams()
	with open(file_name,"rb") as f:
		f.seek(16,0)
		p.dims = np.fromfile(f,dt1,3)
		p.phys_dims = np.fromfile(f,dt2,6)
		f.seek(0,0)
		# for some reason, in python3, the following cannot come before np.fromfile(...)
		p.header = f.read(16) # read header string "DataViewer 2.0.0"
		f.seek(0,2) # go to end of file
		p.file_size = f.tell() # get length of file
		p.floats_per_frame = p.dims[0]*p.dims[1]*p.dims[2]
		p.frames = (p.file_size-52)/(p.floats_per_frame*4)
	return p

def PrintParameters(p):
	"""Prints the parameters of the dvdat file
	inputs:
	p (dvdatParams)"""
	text = ""
	text +="------------------------------------------------------------------\n"
	text +="File properties:\n"
	text +="------------------------------------------------------------------\n"
	text +="header string: " + str(p.header)
	text +="length: " + str(p.file_size)
	text +="frames: " + str(p.frames)
	text +="index dimensions: " + str(p.dims)
	text +='X Bounds: {0:.3f} , {1:.3f}'.format(float(p.phys_dims[0]),float(p.phys_dims[1])) + "\n"
	text +='Y Bounds: {0:.3f} , {1:.3f}'.format(float(p.phys_dims[2]),float(p.phys_dims[3])) + "\n"
	text +='Z Bounds: {0:.3f} , {1:.3f}'.format(float(p.phys_dims[4]),float(p.phys_dims[5])) + "\n"
	text +="------------------------------------------------------------------\n"
	return text

def GetFrameArray(file_name,frame_to_plot):
	"""Returns numpy array containing dvdat frame
	inputs:
	file_name (string)
	frame_to_plot (int)"""
	p = GetParameters(file_name)
	with open(file_name,"rb") as f:
		f.seek(52+p.floats_per_frame*4*frame_to_plot,0) # go to requested frame
		the_data = np.fromfile(f,dt2,p.floats_per_frame,"") # read the frame
		the_data = the_data.reshape(p.dims[2],p.dims[1],p.dims[0]) # dimension the resulting numpy array
	return the_data
	
def GetTimeArray(file_name,idx):
	"""Returns numpy array containing time series from dvdat file
	inputs:
	file_name (string)
	idx (int) = index relative to each frame"""
	p = GetParameters(file_name)
	with open(file_name,"rb") as f:
		the_data = np.zeros([p.frames],dt2)
		for j in range(0,p.frames):
			f.seek(52 + j*p.floats_per_frame*4 + idx*4,0)
			the_data[j] = np.fromfile(f,dt2,1,"")[0]
	return the_data

class dvdatImportPlugin(ImportPlugin):
    """A plugin for reading turboWAVE dvdat files."""

    name = "turboWAVE dvdat plugin"
    author = "Steve Richardson"
    description = "Reads data from turboWAVE dvdat files"

    # Uncomment this line for the plugin to get its own tab
    #promote_tab='turboWAVE'

    file_extensions = set(['.dvdat'])

    def __init__(self):
        ImportPlugin.__init__(self)
        self.fields = [ImportFieldInt("frameNum", descr="Index of frame to import", default=0),]
 
    def doImport(self, params):
        """Actually import data
        params is a ImportPluginParams object.
        Return a list of ImportDataset1D, ImportDataset2D objects
        """
        f = params.filename

        frameNum = params.field_results["frameNum"]
        frameData = GetFrameArray(f, frameNum)
        data = []
        data.append(ImportDataset2D("dvdat", np.squeeze(frameData) ))

        return data
        
    def getPreview(self, params):
        """Import data preview
        params is a ImportPluginParams object.
        Return (text, ok) where text is the text to display and ok is a 
            boolean saying whether the import should be allowed.
        """
        f = params.filename

        p = GetParameters(f)
        text = PrintParameters(p)
        ok = True
        
        return (text, ok)




importpluginregistry.append(dvdatImportPlugin())




