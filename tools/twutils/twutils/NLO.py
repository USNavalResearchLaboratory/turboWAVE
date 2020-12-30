'''
Module :samp:`twutils.NLO`
====================================

This module contains parameters for dielectric materials that can be used in turboWAVE, as well as functions to evaluate their numerical dispersion properties.  This module uses mks units.

Classes and functions
----------------------
'''

import numpy as np
import scipy as sp
from scipy import constants as C

class vacuum:
	"""Contains parameters describing vacuum"""
	num = 0
	Omega = []
	wp = []
	gamma = []
	f = []
	a = []
	b = []
	d = []

class fused_silica:
	"""Contains parameters describing 3 oscillator model for fused silica"""
	num = 3
	Omega = [[2.346e14]*3 , [1.627e16]*3 , [2.776e16]*3]
	wp = [1.794e14 , 1.055e16 , 2.296e16]
	gamma = [[0.0]*3]*3
	f = [[1.0]*3]*3
	a = [[[0.0]*6]*3]*3
	b = [0.0]*3
	d = [0.0]*3

class mag_fluoride:
	"""Contains parameters describing 3 oscillator model for magnesium fluoride"""
	num = 3
	Omega = [4.342e16 , 1.991e16 , 7.917e13]
	wp = [3.032e16 , 1.257e16 , 1.204e14]
	gamma = [[0.0]*3]*3
	f = [[1.0]*3]*3
	a = [[[0.0]*6]*3]*3
	b = [0.0]*3
	d = [0.0]*3

class silicon_carbide:
	"""Contains parameters describing 2 oscillator model for silicon carbide"""
	num = 2
	Omega = [[1.828e14]*3 , [1.885e15]*3]
	wp = [2.733e14 , 4.524e15]
	gamma = [[8.954e11*0.5]*3 , [0.0]*3]
	f = [[1.0]*3]*2
	a = [[[0.0]*6]*3]*2
	b = [0.0]*2
	d = [0.0]*2

class BBO:
	"""Contains parameters describing 2 oscillator model for BBO crystals"""
	num = 2
	Omega = [[1.714e16,1.714e16,1.880e16],[1.590e15,1.590e15,1.356e15]]
	wp = [1e16,1e14]
	gamma = [[0.0]*3]*2
	f = [[5.075,5.075,4.833],[0.6823,0.6823,1.653]]
	a = [
		[[-0.00e+00, -0.00e+00, -0.00e+00, -0.00e+00, -4.89e+40, -1.11e+42],
		[-1.11e+42,  1.11e+42, -0.00e+00, -4.89e+40, -0.00e+00, -0.00e+00],
		[-2.59e+40, -2.59e+40, -3.77e+40, -0.00e+00, -0.00e+00, -0.00e+00]],
		[[0.0]*6]*3
		]
	b = [1.6e54,0.0]
	d = [-7.2e74,0.0]


def analyticWavenumber(w,oscParams,n=1):
	"""Compute k(w) from analytical dispersion relation

	:param float w: angular frequency [mks]
	:param class oscParams: instance of an oscillator class defined in this module
	:param int n: polarization axis numbered from 1-3"""
	ans = w**2 + 0j
	for i in range(oscParams.num):
		ans -= (oscParams.f[i][n]*oscParams.wp[i]**2 * w**2) / (w**2 + 2j*oscParams.gamma[i][n]*w - oscParams.Omega[i][n]**2)
	return np.sqrt(ans)/C.c

def numericalWavenumber(w,oscParams,dt,dz,n=1):
	"""Compute k(w) from numerical dispersion relation

	:param float w: angular frequency [mks]
	:param class oscParams: instance of an oscillator class defined in this module
	:param float dt: timestep [mks]
	:param float dz: cell size [mks]
	:param int n: polarization axis numbered from 1-3"""
	Pol = 0.0
	for i in range(oscParams.num):
		f = oscParams.f[i][n]
		wp = oscParams.wp[i]
		W = oscParams.Omega[i][n]
		g = oscParams.gamma[i][n]
		denom = (1.0+g*dt)*np.sin(0.5*w*dt)**2 - 0.25*(W*dt)**2 + 0.5j*g*dt*np.sin(w*dt)
		Pol += (0.5*f*(wp*dt)**2 * np.sin(0.5*w*dt)**2)/denom
	Vac_deduced = np.cos(w*dt) + Pol
	sin_factor = np.sqrt(0.5*(1.0-Vac_deduced)*dz**2/(C.c*dt)**2)
	return np.arcsin(sin_factor)*2.0/dz

def analyticDF(wa,kz,oscParams,n=1):
	"""Returns analytical dispersion function (intended as argument to ``fsolve``)

	:param array-like wa: two-element array with real & imaginary part of angular frequency
	:param float kz: real wavenumber
	:param class oscParams: instance of an oscillator class defined in this module
	:param int n: polarization axis numbered from 1-3"""
	w = wa[0] + 1j*wa[1]
	ans = w**2 - (C.c*kz)**2
	for i in range(oscParams.num):
		ans -= (oscParams.f[i][n]*oscParams.wp[i]**2 * w**2) / (w**2 + 2j*oscParams.gamma[i][n]*w - oscParams.Omega[i][n]**2)
	return [ans.real,ans.imag]

def numericalDF(wa,kz,oscParams,dt,dz,n=1):
	"""Returns numerical dispersion function (intended as argument to ``fsolve``)

	:param array-like wa: two-element array with real & imaginary part of angular frequency
	:param float kz: real wavenumber
	:param class oscParams: instance of an oscillator class defined in this module
	:param float dt: timestep [mks]
	:param float dz: cell size [mks]
	:param int n: polarization axis numbered from 1-3"""
	w = wa[0] + 1j*wa[1]
	Vac = 1.0 - 2.0 * (C.c*dt)**2 * (np.sin(0.5*kz*dz)**2 / dz**2)
	Pol = 0.0
	for i in range(oscParams.num):
		f = oscParams.f[i][n]
		wp = oscParams.wp[i]
		W = oscParams.Omega[i][n]
		g = oscParams.gamma[i][n]
		denom = (1.0+g*dt)*np.sin(0.5*w*dt)**2 - 0.25*(W*dt)**2 + 0.5j*g*dt*np.sin(w*dt)
		Pol += (0.5*f*(wp*dt)**2 * np.sin(0.5*w*dt)**2)/denom
	ans = np.cos(w*dt) - Vac + Pol
	return [ans.real,ans.imag]
