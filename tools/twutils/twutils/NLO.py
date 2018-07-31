"""Functions to evaluate numerical dispersion in
fully explicit PIC nonlinear optics model"""

import numpy as np
import scipy as sp
from scipy import constants as C

class vacuum:
	"""Contains parameters describing vacuum"""
	num = 0
	Omega = []
	wp = []
	gamma = []
	
class fused_silica:
	"""Contains parameters describing 3 oscillator model for fused silica"""
	num = 3
	Omega = [2.346e14 , 1.627e16 , 2.776e16]
	wp = [1.794e14 , 1.055e16 , 2.296e16]
	gamma = [0.0,0.0,0.0]

class mag_fluoride:
	"""Contains parameters describing 3 oscillator model for magnesium fluoride"""
	num = 3
	Omega = [4.342e16 , 1.991e16 , 7.917e13]
	wp = [3.032e16 , 1.257e16 , 1.204e14]
	gamma = [0.0,0.0,0.0]

class silicon_carbide:
	"""Contains parameters describing 2 oscillator model for silicon carbide"""
	num = 2
	Omega = [1.828e14 , 1.885e15]
	wp = [2.733e14 , 4.524e15]
	gamma = [8.954e11*0.5 , 0.0]
	
def analyticWavenumber(w,oscParams):
	"""Compute k(w) from analytical dispersion relation
	inputs:
	w = angular frequency [mks]
	oscParams = instance of an oscillator class defined in this module"""
	ans = w**2 + 0j
	for i in range(oscParams.num):
		ans -= (oscParams.wp[i]**2 * w**2) / (w**2 + 2j*oscParams.gamma[i]*w - oscParams.Omega[i]**2)
	return np.sqrt(ans)/C.c

def numericalWavenumber(w,oscParams,dt,dz):
	"""Compute k(w) from numerical dispersion relation
	inputs:
	w = angular frequency [mks]
	oscParams = instance of an oscillator class defined in this module
	dt = timestep [mks]
	dz = cell size [mks]"""
	Pol = 0.0
	for i in range(oscParams.num):
		wp = oscParams.wp[i]
		W = oscParams.Omega[i]
		g = oscParams.gamma[i]
		denom = (1.0+g*dt)*np.sin(0.5*w*dt)**2 - 0.25*(W*dt)**2 + 0.5j*g*dt*np.sin(w*dt)
		Pol += (0.5*(wp*dt)**2 * np.sin(0.5*w*dt)**2)/denom
	Vac_deduced = np.cos(w*dt) + Pol
	sin_factor = np.sqrt(0.5*(1.0-Vac_deduced)*dz**2/(C.c*dt)**2)
	return np.arcsin(sin_factor)*2.0/dz

def analyticDF(wa,kz,oscParams):
	"""Returns analytical dispersion function (intended as argument to fsolve(...))
	inputs:
	wa = two-element array with real & imaginary part of angular frequency
	kz = real wavenumber
	oscParams = instance of an oscillator class defined in this module"""
	w = wa[0] + 1j*wa[1]
	ans = w**2 - (C.c*kz)**2
	for i in range(oscParams.num):
		ans -= (oscParams.wp[i]**2 * w**2) / (w**2 + 2j*oscParams.gamma[i]*w - oscParams.Omega[i]**2)
	return [ans.real,ans.imag]

def numericalDF(wa,kz,oscParams,dt,dz):
	"""Returns numerical dispersion function (intended as argument to fsolve(...))
	inputs:
	wa = two-element array with real & imaginary part of angular frequency
	kz = real wavenumber
	oscParams = instance of an oscillator class defined in this module
	dt = timestep [mks]
	dz = cell size [mks]"""
	w = wa[0] + 1j*wa[1]
	Vac = 1.0 - 2.0 * (C.c*dt)**2 * (np.sin(0.5*kz*dz)**2 / dz**2)
	Pol = 0.0
	for i in range(oscParams.num):
		wp = oscParams.wp[i]
		W = oscParams.Omega[i]
		g = oscParams.gamma[i]
		denom = (1.0+g*dt)*np.sin(0.5*w*dt)**2 - 0.25*(W*dt)**2 + 0.5j*g*dt*np.sin(w*dt)
		Pol += (0.5*(wp*dt)**2 * np.sin(0.5*w*dt)**2)/denom
	ans = np.cos(w*dt) - Vac + Pol
	return [ans.real,ans.imag]
	