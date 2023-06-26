'''
Module :samp:`twutils.QO.stationary`
====================================

This module computes stationary states of numerical atoms.
It can be used to generate state files for use with quantum optics modules.

Classes and functions
----------------------
'''

import numpy as np
import scipy as sp
from scipy.sparse.linalg import eigs
from scipy import constants as C
from scipy.optimize import newton

def DiracHypergeometricArgument(w,cylindrical,Qnuc,qorb,morb,nr,jam,jzam,Bz):
	"""Auxilliary function used to find weak field energy levels.  At present cylindrical only.
	Energy is the solution of f(w)==0.

	:parameter float w: energy of the level, depending on usage possibly guess
	:parameter boolean cylindrical: is this a cylincrical atom (currently must be true)
	:parameter float Qnuc: charge of the nucleus
	:parameter float qorb: charge of the orbiting particle
	:parameter float morb: mass of the orbiting particle
	:parameter float nr: radial quantum number
	:parameter float jam: total angular momentum
	:parameter float jzam: total angular momentum projection
	:parameter float Bz: external magnetic field (must be weak)"""
	# write the square roots so we dont get complex numbers
	# we are just trying to get the root finder to not bail
	l0 = np.sqrt(np.abs(morb**2 - w**2))
	Zw = qorb*Qnuc*w/l0
	Zm = qorb*Qnuc*morb/l0
	Omega = qorb*Bz/2
	sigma = np.sqrt(np.abs(jzam**2-(qorb*Qnuc)**2))
	beta = np.sqrt(np.abs(l0**2 - 2*jzam*Omega))
	sgn = -1.0
#	  if jzam>0.0:
#		beta = np.sqrt(np.abs(l0**2 - 2*Omega*(jzam+(sigma+Zw)/(jzam-Zm))))
#		sgn = -1.0
#	  else:
#		beta = np.sqrt(np.abs(l0**2 - 2*Omega*(jzam+(sigma-Zw)/(jzam+Zm))))
#		sgn = -1.0
	if cylindrical:
		return ((sgn+2*Zw)*l0+(1+2*sigma)*beta)/(2*beta) + nr
	else:
		return 0.0

def SchroedingerEnergyLevel(Qnuc,qorb,morb,nr,lam):
	"""Returns Schroedinger energy levels indexed by nr and lam.  The principle quantum number is nr+lam+1.

	:parameter float Qnuc: charge of the nucleus
	:parameter float qorb: charge of the orbiting particle
	:parameter float morb: mass of the orbiting particle
	:parameter int nr: radial quantum number [0,1,2,...]
	:parameter in lam: orbital angular momentum [0,1,2,...]"""
	return morb*(1.0 - (Qnuc*qorb)**2/(2.0*(lam+nr+1.0)**2))

def KGEnergyLevel(cylindrical,Qnuc,qorb,morb,nr,lam,lzam,Bz):
	"""Returns KG energy level indexed by nr,lam,lzam.
	Spherical stability of an electronic atom requires Z<=68.
	Cylindrical stability requires lam>0.

	:parameter boolean cylindrical: is this a cylincrical atom
	:parameter float Qnuc: charge of the nucleus
	:parameter float qorb: charge of the orbiting particle
	:parameter float morb: mass of the orbiting particle
	:parameter int nr: radial quantum number [0,1,2,...]
	:parameter int lam: orbital angular momentum [0,1,2,...]
	:parameter int lzam: orbital angular momentum projection [-lam,...,lam]
	:parameter float Bz: magnetic field (must be weak)"""
	if cylindrical:
		if lzam**2 - (Qnuc*qorb)**2 < 0.0:
			return 0.0
		e2i = 1.0 + (Qnuc*qorb)**2/(np.sqrt(lzam**2 - (Qnuc*qorb)**2) + nr + 0.5)**2
	else:
		if (lam+0.5)**2 - (Qnuc*qorb)**2 < 0.0:
			return 0.0
		e2i = 1.0 + (Qnuc*qorb)**2/(np.sqrt((lam+0.5)**2 - (Qnuc*qorb)**2) + nr + 0.5)**2
	return morb*np.sqrt(1.0-lzam*qorb*Bz/morb**2)/np.sqrt(e2i)

def DiracEnergyLevel(cylindrical,Qnuc,qorb,morb,nr,jam,jzam,Bz,energy_estimate=0.99):
	"""Returns Dirac energy levels indexed by nr,jam,jzam.
	Spherical stability of an electronic atom requires Z<=137.
	Cylindrical stability requires Z<=68.

	:parameter boolean cylindrical: is this a cylincrical atom
	:parameter float Qnuc: charge of the nucleus
	:parameter float qorb: charge of the orbiting particle
	:parameter float morb: mass of the orbiting particle
	:parameter int nr: radial quantum number [0,1,2,...]
	:parameter float jam: total angular momentum [0.5,1.5,2.5,...]
	:parameter float jzam: total angular momentum projection [-jam,...,jam]
	:parameter float Bz: magnetic field (must be weak)"""
	if Bz==0:
		x2 = (Qnuc*qorb)**2
		if cylindrical:
			if jzam**2 - x2 < 0.0:
				return 0.0
			e2i = 1.0 + x2/(np.sqrt(jzam**2-x2) + nr)**2
		else:
			if (jam+0.5)**2 - x2 < 0.0:
				return 0.0
			e2i = 1.0 + x2/(np.sqrt((jam+0.5)**2 - x2) + nr)**2
		return morb/np.sqrt(e2i)
	else:
		return newton(DiracHypergeometricArgument,energy_estimate,args=(cylindrical,Qnuc,qorb,morb,nr,jam,jzam,Bz))

def nu_to_au(energy):
	'''take energy from natural to atomic units'''
	return (energy - 1.0)*C.m_e*C.c**2/(2.0*C.Rydberg*C.c*C.h)

def nu_to_eV(energy):
	'''take energy from natural units to eV'''
	return (energy - 1.0)*C.m_e*C.c**2/C.e

def nu_to_keV(energy):
	'''take energy from natural units to keV'''
	return 1e-3*(energy - 1.0)*C.m_e*C.c**2/C.e

def nu_to_MeV(energy):
	'''take energy from natural units to MeV'''
	return 1e-6*(energy - 1.0)*C.m_e*C.c**2/C.e

def SoftCorePotential(Qnuc,r0,num_pts,dr):
	"""Create scalar potential on a radial grid

	:parameter float Qnuc: charge of the nucleus
	:parameter float r0: soft core radius
	:parameter int num_pts: number of interior grid points, total is num_pts+2
	:parameter float dr: radial grid spacing
	:returns: an array of radial points and an array of potentials, including ghost cells
	:returnType: tuple with the two arrays (grid,phi)"""
	grid = np.linspace(-0.5*dr,num_pts*dr + 0.5*dr,num_pts+2)
	phi = Qnuc/np.sqrt(r0**2 + grid**2)
	return grid , phi

def KGMatrix(cylindrical,grid,phi,qorb,morb,lam,lzam,Bz):
	"""Returns matrix of the discretized time-independent KG operator.
	This uses the leapfrog friendly Hamiltonian representation.
	For N grid points we obtain a 2Nx2N matrix corresponding to 2 components.
	The second component is :math:`\chi = (\omega-q\phi)\psi/m`

	:parameter boolean cylindrical: is this a cylincrical atom
	:parameter array grid: grid array, can get from SoftCorePotential function
	:parameter array phi: potential array, can get from SoftCorePotential function
	:parameter float qorb: charge of the orbiting particle
	:parameter float morb: mass of the orbiting particle
	:parameter int lam: orbital angular momentum [0,1,2,...], ignored in cylindrical case
	:parameter int lzam: orbital angular momentum projection [-lam,...,lam], or any integer in cylindrical case
	:parameter float Bz: magnetic field (must be weak in spherical case)"""
	N = grid.shape[0]
	num_pts = N-2
	dr = grid[1]-grid[0]
	r1_array = grid[1:N-1] - 0.5*dr
	r2_array = grid[1:N-1]
	r3_array = grid[1:N-1] + 0.5*dr
	phi_array = phi[1:N-1]
	if cylindrical:
		T1 = -(1.0/dr**2)*r1_array[1:]/r2_array[1:]
		T2 = morb*morb + lzam**2/r2_array**2 + (1.0/dr**2)*(r1_array + r3_array)/r2_array
		T2 -= lzam*qorb*Bz - (0.5*qorb*Bz*r2_array)**2
		T3 = -(1.0/dr**2)*r3_array[:num_pts-1]/r2_array[:num_pts-1]
	else:
		T1 = -(1.0/dr)*3.0*r1_array[1:]**2/(r3_array[1:]**3-r1_array[1:]**3)
		T2 = morb*morb + lam*(lam+1)/r2_array**2 + (1.0/dr)*3.0*(r1_array**2 + r3_array**2)/(r3_array**3-r1_array**3)
		T2 -= lzam*qorb*Bz # have to neglect B^2 or we would need to solve in r-theta plane
		T3 = -(1.0/dr)*3.0*r3_array[:num_pts-1]**2/(r3_array[:num_pts-1]**3-r1_array[:num_pts-1]**3)
	# may have to use "spdiags" with older scipy
	mA = sp.sparse.diags([qorb*phi_array],[0])
	mB = sp.sparse.diags([np.ones(num_pts)*morb],[0])
	mC = sp.sparse.diags([T1/morb,T2/morb,T3/morb],[-1,0,1])
	mD = sp.sparse.diags([qorb*phi_array],[0])
	return sp.sparse.bmat([[mA,mB],[mC,mD]])

def DiracMatrix(cylindrical,grid,phi,qorb,morb,jam,lam,jzam,Bz):
	"""Returns matrix of the discretized time-independent Dirac operator.
	For N grid points we obtain a 2Nx2N matrix corresponding to 2 components.
	The two components are :math:`f` and :math:`g` which in turn form the bispinor.
	See, e.g., Landau and Lifshitz QED section 35.

	:parameter boolean cylindrical: is this a cylincrical atom
	:parameter array grid: grid array, can get from SoftCorePotential function
	:parameter array phi: potential array, can get from SoftCorePotential function
	:parameter float qorb: charge of the orbiting particle
	:parameter float morb: mass of the orbiting particle
	:parameter float jam: total angular momentum [0.5,1.5,...], ignored in cylindrical case
	:parameter int lam: parity number, must be jam :math:`\pm 1/2`, or in cylindrical case, jzam :math:`\pm 1/2`
	:parameter float jzam: total angular momentum projection [-jam,...,jam], or in cylindrical case, any half integer
	:parameter float Bz: magnetic field (must be weak in spherical case)"""
	N = grid.shape[0]
	num_pts = N-2
	# uniform grid is assumed
	# may have to use "spdiags" with older scipy
	dr = grid[1:N-1]-grid[0:N-2]
	rc = grid[1:N-1]
	phi_array = phi[1:N-1]
	if cylindrical:
		kappa = 2*(lam-jzam)*(jzam - 0.5*qorb*Bz*rc*rc)
		D1 = np.complex(-0.5)/dr[1:]
		D2f = np.complex(0.5)/rc
		D2g = np.complex(0.5)/rc
		D3 = np.complex(0.5)/dr[:num_pts-1]
		# Boundary condition based on parity number
		D2f[0] += D1[0]*(-1)**lam
		D2g[0] -= D1[0]*(-1)**lam
		mA = sp.sparse.diags([np.complex(qorb)*phi_array+morb],[0])
		mB = sp.sparse.diags([-1j*D1,-1j*(D2g-kappa/rc),-1j*D3],[-1,0,1])
		mC = sp.sparse.diags([-1j*D1,-1j*(D2f+kappa/rc),-1j*D3],[-1,0,1])
		mD = sp.sparse.diags([np.complex(qorb)*phi_array-morb],[0])
	else:
		if Bz!=0.0:
			exit(1)
		# spherical case ignores B-field
		kappa = 2*(lam-jam)*(jam+0.5)
		D1 = -0.5/dr[1:]
		D2f = 1.0/rc
		D2g = 1.0/rc
		D3 = 0.5/dr[:num_pts-1]
		# Boundary condition based on parity number
		# But what are the parities of f & g ??
		D2f[0] += D1[0]*(-1)**lam
		D2g[0] -= D1[0]*(-1)**lam
		mA = sp.sparse.diags([np.complex(qorb)*phi_array+morb],[0])
		mB = sp.sparse.diags([-D1,-D2g+kappa/rc,-D3],[-1,0,1])
		mC = sp.sparse.diags([D1,D2f+kappa/rc,D3],[-1,0,1])
		mD = sp.sparse.diags([np.complex(qorb)*phi_array-morb],[0])
	return sp.sparse.bmat([[mA,mB],[mC,mD]])

class StationaryStateTool:
	'''Class for managing stationary states, and creating state files'''

	def __init__(self,cylindrical,grid,phi,nr,jam,lam,jzam):
		'''Create a stationary state tool.

		:parameter bool cylindrical: is this a cylindrical atom
		:parameter array grid: radial grid, can be output from SoftCorePotential()
		:parameter array phi: potential on the radial grid, can be output from SoftCorePotential()
		:parameter int nr: radial quantum number [0,1,2,...] (cannot be 0 in certain cases)
		:parameter float jam: see lower level functions for possible values
		:parameter int lam: see lower level functions for possible values
		:parameter float jzam: see lower level functions for possible values.  This triggers Dirac or KG treatment depending on whether integer or half-integer.'''
		self.cylindrical = cylindrical
		self.grid = grid
		self.phi = phi
		self.nr = nr
		self.jam = jam
		self.lam = lam
		self.jzam = jzam
		self.qorb = -np.sqrt(C.alpha)
		self.morb = 1.0
		self.Bz = 0.0

	def SetSpecialConditions(self,qorb,morb,Bz):
		self.qorb = qorb
		self.morb = morb
		self.Bz = Bz

	def RescalePhi(self,scale_factor):
		self.phi *= scale_factor

	def CheckQuantumNumbers(self):
		nr = self.nr
		jam = self.jam
		lam = self.lam
		jzam = self.jzam
		if np.abs(np.modf(jzam)[0])!=0.5 and np.abs(np.modf(jzam)[0])!=0.0:
			return 'jzam is neither integer nor half-integer'

		# Spin 1/2 Case
		if np.abs(np.modf(jzam)[0])==0.5:
			if self.cylindrical:
				if np.abs(jzam-lam)!=0.5:
					return 'abs(jzam-lam) is not 1/2.'
				if nr==0.0 and (lam-jzam)*jzam>0.0:
					return 'nr cannot be 0 with requested angular mode.'
			else:
				if np.abs(np.modf(jam)[0])!=0.5 or jam<0.5:
					return 'jam is not a positive half-integer.'
				if np.abs(jam-lam)!=0.5:
					return 'abs(jam-lam) is not 1/2.'
				if nr==0.0 and (lam-jam)*(jam+0.5)>0.0:
					return 'nr cannot be 0 with requested angular mode.'
				if jzam<-jam or jzam>jam:
					return 'violation of -jam<=jzam<=jam'

		# Spin 0 Case
		if np.abs(np.modf(jzam)[0])==0.0:
			if not self.cylindrical:
				if np.abs(np.modf(lam)[0])!=0.0 or lam<0.0:
					return 'lam is not a positive integer or zero.'
				if jzam<-lam or jzam>lam:
					return 'violation of -lam<=jzam<=lam'

		return 'OK'

	def spin0(self):
		return np.abs(np.modf(self.jzam)[0])==0.0

	def GetEnergy(self,energy_guess,scale_factor):
		nr_max = int(self.nr + 1)
		# Find requested state; sorting energies identifies nr
		if self.spin0():
			mat = KGMatrix(self.cylindrical,self.grid,scale_factor*self.phi,self.qorb,self.morb,self.lam,self.jzam,self.Bz)
		else:
			mat = DiracMatrix(self.cylindrical,self.grid,scale_factor*self.phi,self.qorb,self.morb,self.jam,self.lam,self.jzam,self.Bz)
		vals,vecs = eigs(mat.tocsc(),k=nr_max,sigma=energy_guess,which='LM')
		sorted_indices = np.argsort(vals)
		idx = sorted_indices[int(self.nr)]
		En = np.real(vals[idx])
		return En

	def GetBaseComponents(self,energy_guess):
		'''Gets 2 components used to form the wavefunction, call them f and g.
		For KG returns a Hamiltonian decomposition, with f the usual scalar, and g an auxiliary function.
		The Feshbach-Villars decomposition is then :math:`(f+g)/\sqrt{2}` , :math:`(f-g)/\sqrt{2}`.
		For Dirac returns f,g from L&L, which may be used to form the bispinor.

		:parameter float energy_guess: initial guess of the energy
		:returns: En,f,g,sorted (En is the selected eigenvalue, sorted are the first few eigenvalues)'''

		num_pts = self.grid.shape[0]-2
		nr_max = int(self.nr + 4)

		# Find requested state; sorting energies identifies nr

		if self.spin0():
			mat = KGMatrix(self.cylindrical,self.grid,self.phi,self.qorb,self.morb,self.lam,self.jzam,self.Bz)
		else:
			mat = DiracMatrix(self.cylindrical,self.grid,self.phi,self.qorb,self.morb,self.jam,self.lam,self.jzam,self.Bz)
		vals,vecs = eigs(mat.tocsc(),k=nr_max,sigma=energy_guess,which='LM')

		sorted_indices = np.argsort(vals)
		idx = sorted_indices[int(self.nr)]
		f = vecs[:num_pts,idx]
		g = vecs[num_pts:,idx]
		En = np.real(vals[idx])
		sorted = np.sort(np.real(vals))
		return En,f,g,sorted

	def GetWavefunction(self,f,g):
		'''Return 4 components (0,1,2,3) derived from the base components f and g.
		For KG only components 1 and 2 are relevant, components 0 and 3 may be ignored.
		For KG we return the simple decomposition with psi1 the usual scalar wavefunction.
		The Feshbach-Villars decomposition is then :math:`(\psi_1+\psi_2)/\sqrt{2}` , :math:`(\psi_1-\psi_2)/\sqrt{2}`.
		For Dirac the 4 components are the radial part of the bispinor.'''
		nr = self.nr
		jam = self.jam
		lam = self.lam
		jzam = self.jzam
		num_pts = self.grid.shape[0]-2
		dr = self.grid[1]-self.grid[0] # assume uniform grid
		if self.cylindrical:
			jacobian = 2*np.pi*self.grid[1:num_pts+1]
		else:
			jacobian = self.grid[1:num_pts+1]**2
		if self.spin0():
			psi0 = np.zeros(num_pts)
			psi1 = f
			psi2 = g
			psi3 = np.zeros(num_pts)
			rho = 0.5*(np.abs(f+g)**2 - np.abs(f-g)**2)
		else:
			if self.cylindrical:
				if jzam>lam:
					psi0 = f
					psi1 = np.zeros(num_pts)
					psi2 = np.zeros(num_pts)
					psi3 = g
				else:
					psi0 = np.zeros(num_pts)
					psi1 = f
					psi2 = g
					psi3 = np.zeros(num_pts)
			else:
				psi0 = f
				psi1 = f
				psi2 = (-1.0)**((1+2*lam-2*jam)/2)*g
				psi3 = (-1.0)**((1+2*lam-2*jam)/2)*g
			rho = np.abs(psi0)**2 + np.abs(psi1)**2 + np.abs(psi2)**2 + np.abs(psi3)**2

		CN = np.sqrt(1/np.sum(dr*rho*jacobian))
		psi0 = CN*psi0
		psi1 = CN*psi1
		psi2 = CN*psi2
		psi3 = CN*psi3

		return psi0,psi1,psi2,psi3

	def StandardFilename(self,Z):
		'''Return a filename appropriate for the current quantum state

		:parameter float Z: atomic number, purely for labeling
		:returns: string with the filename'''
		file_name = 'z'+str(int(Z))
		if self.cylindrical:
			file_name += '_c'
			if self.spin0():
				angular_num = self.jzam
			else:
				angular_num = self.lam
		else:
			file_name += '_s'
			angular_num = self.lam
		file_name += str(int(self.nr))
		file_name += str(int(angular_num))
		if not self.spin0():
			if self.cylindrical:
				if self.jzam>self.lam:
					file_name += '+'
				else:
					file_name += '-'
			else:
				if self.jam>self.lam:
					file_name += '+'
				else:
					file_name += '-'
		file_name += '.txt'
		return file_name

	def WriteFile(self,file_name,psi0,psi1,psi2,psi3,En,softCoreRadius,Qnuc):
		'''Write the state file given the current quantum state'''
		num_pts = self.grid.shape[0]-2
		dr = self.grid[1]-self.grid[0] # assume uniform grid
		outFile = open(file_name,'w')
		outFile.write('energy = ')
		outFile.write(str(En))
		outFile.write('\n')
		outFile.write('pts = ')
		outFile.write(str(num_pts))
		outFile.write('\n')
		if self.spin0():
			outFile.write('components = 2')
		else:
			outFile.write('components = 4')
		outFile.write('\n')
		outFile.write('cell_width = ')
		outFile.write(str(dr))
		outFile.write('\n')
		outFile.write('soft_core_radius = ')
		outFile.write(str(softCoreRadius))
		outFile.write('\n')
		outFile.write('nuclear_charge = ')
		outFile.write(str(Qnuc))
		outFile.write('\n')
		outFile.write('nr_j_l_m = ')
		outFile.write(str(self.nr)+' ')
		outFile.write(str(self.jam)+' ')
		outFile.write(str(self.lam)+' ')
		outFile.write(str(self.jzam)+ ' ')
		outFile.write('\n')
		outFile.write('Bz = ')
		outFile.write(str(self.Bz))
		outFile.write('\n')
		if self.cylindrical:
			outFile.write('cylindrical = true')
		else:
			outFile.write('cylindrical = false')
		outFile.write('\n\n')
		outFile.write('start_data\n')
		if self.spin0():
			np.column_stack((np.real(psi1),np.imag(psi1))).tofile(outFile,'\n')
			outFile.write('\n')
			np.column_stack((np.real(psi2),np.imag(psi2))).tofile(outFile,'\n')
		else:
			np.column_stack((np.real(psi0),np.imag(psi0))).tofile(outFile,'\n')
			outFile.write('\n')
			np.column_stack((np.real(psi1),np.imag(psi1))).tofile(outFile,'\n')
			outFile.write('\n')
			np.column_stack((np.real(psi2),np.imag(psi2))).tofile(outFile,'\n')
			outFile.write('\n')
			np.column_stack((np.real(psi3),np.imag(psi3))).tofile(outFile,'\n')
		outFile.close()
