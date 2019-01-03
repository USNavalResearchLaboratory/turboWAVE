import sys
import numpy as np
import re
import OCC.Display.SimpleGui as occgui
import OCC.BRepPrimAPI as occprim
import OCC.BRepBuilderAPI as occbuild
import OCC.BRepAlgoAPI as occalgo
from OCC.gp import *

# Things needed for STEP or STL export
from OCC.STEPControl import STEPControl_Writer, STEPControl_AsIs
from OCC.StlAPI import StlAPI_Writer
from OCC.Interface import Interface_Static_SetCVal
from OCC.IFSelect import IFSelect_RetDone

# The master shape list
tw_shape = []

# Function to strip comments from TW input file
def comment_remover(text):
	# From stackoverflow.com, author Markus Jarderot
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

# Set up python classes to parallel turboWAVE classes

def RotateZ(v,q):
	st = np.sin(q);
	ct = np.cos(q);
	xOld = v[0];
	v[0] = xOld*ct - v[1]*st;
	v[1] = xOld*st + v[1]*ct;

def RotateY(v,q):
	st = np.sin(q);
	ct = np.cos(q);
	zOld =  v[2];
	v[2] = zOld*ct - v[0]*st;
	v[0] = zOld*st + v[0]*ct;

def RotateX(v,q):
	st = np.sin(q);
	ct = np.cos(q);
	yOld = v[1];
	v[1] = yOld*ct - v[2]*st;
	v[2] = yOld*st + v[2]*ct;

class basis:
	def __init__(self):
		self.u = np.array([1,0,0]).astype(np.float)
		self.v = np.array([0,1,0]).astype(np.float)
		self.w = np.array([0,0,1]).astype(np.float)
	def AboutX(self,q):
		RotateX(self.u,q)
		RotateX(self.v,q)
		RotateX(self.w,q)
	def AboutY(self,q):
		RotateY(self.u,q)
		RotateY(self.v,q)
		RotateY(self.w,q)
	def AboutZ(self,q):
		RotateZ(self.u,q)
		RotateZ(self.v,q)
		RotateZ(self.w,q)

class tw_rgn:
	type = 'none'
	def __init__(self,name):
		self.name = name
		self.center = np.array([0,0,0]).astype(np.float)
		self.rbox = np.array([0,0,0]).astype(np.float)
		self.complement = False
		self.orientation = basis()
		self.occ_shape = None
		self.transparency = 0.0
		self.elements = []
	def ParseItem(self,words,i):
		if words[i]=='bounds':
			p1 = np.array([np.float(words[i+1]),np.float(words[i+3]),np.float(words[i+5])])
			p2 = np.array([np.float(words[i+2]),np.float(words[i+4]),np.float(words[i+6])])
			self.center = 0.5*(p1 + p2)
			self.rbox = 0.5*(p2 - p1)
			return 7
		if words[i]=='center':
			print('ERROR: Direct assignment to center illegal.  Use translation instead.')
			exit(1)
		if words[i]=='translation':
			self.center[0] += np.float(words[i+1])
			self.center[1] += np.float(words[i+2])
			self.center[2] += np.float(words[i+3])
			return 4
		if words[i]=='rotation':
			if words[i+2]=='x':
				self.orientation.AboutX(np.float(words[i+3])*np.pi/180)
				RotateX(self.center,np.float(words[i+3])*np.pi/180)
			if words[i+2]=='y':
				self.orientation.AboutY(np.float(words[i+3])*np.pi/180)
				RotateY(self.center,np.float(words[i+3])*np.pi/180)
			if words[i+2]=='z':
				self.orientation.AboutZ(np.float(words[i+3])*np.pi/180)
				RotateZ(self.center,np.float(words[i+3])*np.pi/180)
			return 4
		if words[i]=='complement':
			self.complement = words[i+1]=='true' or words[i+1]=='yes' or words[i+1]=='on'
			return 2
		return 0
	def ParseElementList(self,words,i):
		if words[i]!='{':
			print('Format error in input file (element list initiation).')
			exit(1)
		i += 1
		while words[i]!='}':
			self.elements.append(next(x for x in tw_shape if x.name==words[i]))
			i += 1
	def Parse(self,words,i):
		if words[i]!='{':
			print('Format error in input file (block initiation).')
			exit(1)
		i += 1
		while words[i]!='}':
			i0 = i
			i += self.ParseItem(words,i)
			if i==i0:
				print('Format error in input file (block termination):',self.name)
				exit(1)
	def Name(self):
		return self.name
	def Get_center(self):
		return gp_Pnt(self.center[0],self.center[1],self.center[2])
	def OrientShape(self):
		# Intended to be called at end of CreateShape
		# Shape is expected to be in default position and orientation
		T = gp_Trsf()
		u = self.orientation.u
		v = self.orientation.v
		w = self.orientation.w
		#T.SetValues(u[0],u[1],u[2],0,v[0],v[1],v[2],0,w[0],w[1],w[2],0)
		T.SetValues(u[0],v[0],w[0],0,u[1],v[1],w[1],0,u[2],v[2],w[2],0)
		self.occ_shape = occbuild.BRepBuilderAPI_Transform(self.occ_shape,T,False).Shape()
		T.SetTranslation(gp_Vec(self.center[0],self.center[1],self.center[2]))
		self.occ_shape = occbuild.BRepBuilderAPI_Transform(self.occ_shape,T,False).Shape()
	def CreateShape(self):
		p1 = gp_Pnt(-self.rbox[0],-self.rbox[1],-self.rbox[2])
		p2 = gp_Pnt(self.rbox[0],self.rbox[1],self.rbox[2])
		self.occ_shape = occprim.BRepPrimAPI_MakeBox(p1,p2).Shape()
		self.OrientShape()

class tw_domain(tw_rgn):
	type = 'domain'
	def ParseItem(self,words,i):
		if words[i]=='corner':
			self.corner = np.array([0,0,0]).astype(np.float)
			self.corner[0] = np.float(words[i+1])
			self.corner[1] = np.float(words[i+2])
			self.corner[2] = np.float(words[i+3])
			return 4
		if words[i]=='dimensions':
			self.dim = np.array([np.float(words[i+1]),np.float(words[i+2]),np.float(words[i+3])])
			return 4
		if words[i]=='cell':
			self.dr = np.array([np.float(words[i+2]),np.float(words[i+3]),np.float(words[i+4])])
			return 5
		if words[i]=='geometry':
			return 2
		if words[i]=='decomposition':
			return 4
		if words[i]=='region':
			return 8
		if words[i]=='radial':
			return 4
		if words[i]=='adaptive':
			return 3
		return 0
	def Parse(self,words,i):
		tw_rgn.Parse(self,words,i)
		p1 = self.corner
		p2 = self.corner + self.dim*self.dr
		self.center = 0.5*(p1 + p2)
		self.rbox = 0.5*(p2 - p1)
	def CreateShape(self):
		tw_rgn.CreateShape(self)
		self.transparency = 0.7

class tw_rect(tw_rgn):
	type = 'rect'

class tw_circ(tw_rgn):
	type = 'circ'
	def ParseItem(self,words,i):
		if words[i]=='radius':
			self.radius = np.float(words[i+1])
			return 2
		return tw_rgn.ParseItem(self,words,i)
	def CreateShape(self):
		self.occ_shape = occprim.BRepPrimAPI_MakeSphere(gp_Pnt(0,0,0),self.radius).Shape()
		self.OrientShape()

class tw_cyl(tw_rgn):
	type = 'cylinder'
	def ParseItem(self,words,i):
		if words[i]=='radius':
			self.radius = np.float(words[i+1])
			return 2
		if words[i]=='length':
			self.length = np.float(words[i+1])
			return 2
		return tw_rgn.ParseItem(self,words,i)
	def CreateShape(self):
		self.occ_shape = occprim.BRepPrimAPI_MakeCylinder(self.radius,self.length).Shape()
		T = gp_Trsf()
		T.SetTranslation(gp_Vec(0,0,-self.length/2))
		self.occ_shape = occbuild.BRepBuilderAPI_Transform(self.occ_shape,T,False).Shape()
		self.OrientShape()

class tw_rounded_cyl(tw_rgn):
	type = 'rounded_cylinder'
	def ParseItem(self,words,i):
		if words[i]=='radius':
			self.radius = np.float(words[i+1])
			return 2
		if words[i]=='length':
			self.length = np.float(words[i+1])
			return 2
		return tw_rgn.ParseItem(self,words,i)
	def CreateShape(self):
		endcap1 = occprim.BRepPrimAPI_MakeSphere(gp_Pnt(0,0,0),self.radius).Shape()
		endcap2 = occprim.BRepPrimAPI_MakeSphere(gp_Pnt(0,0,self.length),self.radius).Shape()
		body = occprim.BRepPrimAPI_MakeCylinder(self.radius,self.length).Shape()
		self.occ_shape = occalgo.BRepAlgoAPI_Fuse(endcap1,body).Shape()
		self.occ_shape = occalgo.BRepAlgoAPI_Fuse(self.occ_shape,endcap2).Shape()
		T = gp_Trsf()
		T.SetTranslation(gp_Vec(0,0,-self.length/2))
		self.occ_shape = occbuild.BRepBuilderAPI_Transform(self.occ_shape,T,False).Shape()
		self.OrientShape()

class tw_cyl_shell(tw_rgn):
	type = 'cylindrical_shell'
	def ParseItem(self,words,i):
		if words[i]=='inner':
			self.inner_radius = np.float(words[i+2])
			return 3
		if words[i]=='outer':
			self.outer_radius = np.float(words[i+2])
			return 3
		if words[i]=='length':
			self.length = np.float(words[i+1])
			return 2
		return tw_rgn.ParseItem(self,words,i)
	def CreateShape(self):
		outer = occprim.BRepPrimAPI_MakeCylinder(self.outer_radius,self.length).Shape()
		inner = occprim.BRepPrimAPI_MakeCylinder(self.inner_radius,self.length).Shape()
		self.occ_shape = occalgo.BRepAlgoAPI_Cut(outer,inner).Shape()
		T = gp_Trsf()
		T.SetTranslation(gp_Vec(0,0,-self.length/2))
		self.occ_shape = occbuild.BRepBuilderAPI_Transform(self.occ_shape,T,False).Shape()
		self.OrientShape()

class tw_torus(tw_rgn):
	type = 'torus'
	def ParseItem(self,words,i):
		if words[i]=='minor':
			self.minor_radius = np.float(words[i+2])
			return 3
		if words[i]=='major':
			self.major_radius = np.float(words[i+2])
			return 3
		return tw_rgn.ParseItem(self,words,i)
	def CreateShape(self):
		self.occ_shape = occprim.BRepPrimAPI_MakeTorus(self.major_radius,self.minor_radius).Shape()
		self.OrientShape()

class tw_cone(tw_rgn):
	type = 'cone'
	def ParseItem(self,words,i):
		if words[i]=='length':
			self.length = np.float(words[i+1])
			return 2
		if words[i]=='tip':
			self.tip_radius = np.float(words[i+2])
			return 3
		if words[i]=='base':
			self.base_radius = np.float(words[i+2])
			return 3
		return tw_rgn.ParseItem(self,words,i)
	def CreateShape(self):
		self.occ_shape = occprim.BRepPrimAPI_MakeCone(self.base_radius,self.tip_radius,self.length).Shape()
		T = gp_Trsf()
		T.SetTranslation(gp_Vec(0,0,-self.length/2))
		self.occ_shape = occbuild.BRepBuilderAPI_Transform(self.occ_shape,T,False).Shape()
		self.OrientShape()

class tw_union(tw_rgn):
	type = 'union'
	def ParseItem(self,words,i):
		if words[i]=='elements':
			self.ParseElementList(words,i+1)
			return len(self.elements)+3
		return tw_rgn.ParseItem(self,words,i)
	def CreateShape(self):
		self.occ_shape = self.elements[0].occ_shape
		for s in self.elements:
			if s.occ_shape!=self.occ_shape:
				self.occ_shape = occalgo.BRepAlgoAPI_Fuse(self.occ_shape,s.occ_shape).Shape()
		self.OrientShape()

class tw_intersection(tw_rgn):
	type = 'union'
	def ParseItem(self,words,i):
		if words[i]=='elements':
			self.ParseElementList(words,i+1)
			return len(self.elements)+3
		return tw_rgn.ParseItem(self,words,i)
	def CreateShape(self):
		self.occ_shape = self.elements[0].occ_shape
		if self.elements[0].complement:
			print('ERROR: complementary shape cannot be first in list')
			exit(1)
		for s in self.elements:
			if s.occ_shape!=self.occ_shape:
				if s.complement:
					self.occ_shape = occalgo.BRepAlgoAPI_Cut(self.occ_shape,s.occ_shape).Shape()
				else:
					self.occ_shape = occalgo.BRepAlgoAPI_Common(self.occ_shape,s.occ_shape).Shape()
		self.OrientShape()

# Read the input file into a list of words
with open("stdin","r") as f:
	data = f.read()
	data = comment_remover(data)
	words = data.replace("="," ").replace(","," ").replace("("," ").replace(")"," ").split()

# Process unit conversion macros
# First get the unit of density
i = 0
# Default to ~1 um unit
unitDensityCGS = 2.8e19
while i<len(words):
	if words[i]=='unit':
		if words[i+1]=='density':
			unitDensityCGS = np.double(words[i+2])
	i += 1
omegaP2 = 4.0*np.pi*unitDensityCGS*4.803e-10**2 / 9.109e-28
lengthUnitMKS = 2.998e8/np.sqrt(omegaP2)
# Apply the macro substitution for length units only
i = 0
while i<len(words):
	if words[i][0]=='%':
		old_word = words[i]
		temp = words[i].split('%')[1]
		if temp[-1]=='m':
			if temp[-2]=='u':
				words[i] = str(1e-6*np.float(temp[:-2])/lengthUnitMKS)
			if temp[-2]=='m':
				words[i] = str(1e-3*np.float(temp[:-2])/lengthUnitMKS)
			if temp[-2]=='c':
				words[i] = str(1e-2*np.float(temp[:-2])/lengthUnitMKS)
			if temp[-2].isdecimal():
				words[i] = str(np.float(temp[:-1])/lengthUnitMKS)
			print('Substitution:',old_word,'replaced with',words[i])
	i += 1

# Parse the input file to populate the shape list
i = 0
while i<len(words):
	if words[i]=='new':
		if words[i+1]=='grid':
			tw_shape.append(tw_domain('domain'))
			tw_shape[-1].Parse(words,i+2)
		if words[i+1]=='region':
			rgnType = words[i+2]
			if rgnType=='rect':
				tw_shape.append(tw_rect(words[i+3]))
				tw_shape[-1].Parse(words,i+4)
			if rgnType=='circ':
				tw_shape.append(tw_circ(words[i+3]))
				tw_shape[-1].Parse(words,i+4)
# 			if rgnType=='prism':
# 				tw_shape.append(tw_rect(words[i+3]))
# 				tw_shape[-1].Parse(words,i+4)
# 			if rgnType=='ellipsoid':
# 				tw_shape.append(tw_rect(words[i+3]))
# 				tw_shape[-1].Parse(words,i+4)
			if rgnType=='cylinder':
				tw_shape.append(tw_cyl(words[i+3]))
				tw_shape[-1].Parse(words,i+4)
			if rgnType=='rounded_cylinder':
				tw_shape.append(tw_rounded_cyl(words[i+3]))
				tw_shape[-1].Parse(words,i+4)
			if rgnType=='true_sphere':
				tw_shape.append(tw_circ(words[i+3]))
				tw_shape[-1].Parse(words,i+4)
# 			if rgnType=='box_array':
# 				tw_shape.append(tw_rect(words[i+3]))
# 				tw_shape[-1].Parse(words,i+4)
			if rgnType=='torus':
				tw_shape.append(tw_torus(words[i+3]))
				tw_shape[-1].Parse(words,i+4)
			if rgnType=='cone':
				tw_shape.append(tw_cone(words[i+3]))
				tw_shape[-1].Parse(words,i+4)
# 			if rgnType=='tangent_ogive':
# 				tw_shape.append(tw_rect(words[i+3]))
# 				tw_shape[-1].Parse(words,i+4)
			if rgnType=='cylindrical_shell':
				tw_shape.append(tw_cyl_shell(words[i+3]))
				tw_shape[-1].Parse(words,i+4)
			if rgnType=='union':
				tw_shape.append(tw_union(words[i+3]))
				tw_shape[-1].Parse(words,i+4)
			if rgnType=='intersection':
				tw_shape.append(tw_intersection(words[i+3]))
				tw_shape[-1].Parse(words,i+4)
	i += 1

# Use the turboWAVE shapes to create CAD shapes
for shape in tw_shape:
	print(shape.name,[x.name for x in shape.elements])
	shape.CreateShape()

# Create the display window
display,start_display,add_menu,add_function_to_menu=occgui.init_display()
for shape in tw_shape:
	if len(sys.argv)==1 or shape.Name() in sys.argv or shape.Name()=='domain':
		display.DisplayShape(shape.occ_shape,update=True,transparency=shape.transparency)
start_display()

# Export geometry to STEP or STL if requested
lastArg = sys.argv[-1].split('=')
cmd = lastArg[0]
if cmd=='stepfile' or cmd=='stlfile':
	if len(lastArg)==1:
		print('ERROR: no name was given for the requested output file.')
		exit(1)
	filename = lastArg[1]
	if cmd=='stepfile':
		step_writer = STEPControl_Writer()
		Interface_Static_SetCVal("write.step.schema","AP203")
		for shape in tw_shape:
			if shape.Name() in sys.argv:
				step_writer.Transfer(shape.occ_shape,STEPControl_AsIs)
				status = step_writer.Write(filename+'.step')
				assert(status==IFSelect_RetDone)
				exit(0)
	if cmd=='stlfile':
		stl_writer = StlAPI_Writer()
		stl_writer.SetASCIIMode(True)
		for shape in tw_shape:
			if shape.Name() in sys.argv:
				status = stl_writer.Write(shape.occ_shape,filename+'.stl')
				exit(0)
	print('ERROR: failed to write requested ouptut file.')
	print('Did you list a valid region on the command line?')
