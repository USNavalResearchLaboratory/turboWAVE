'''
Module :samp:`twutils.pre`
--------------------------

This module contains functions and classes to help with setting parameters for turboWAVE input files.'''

import numpy as np
import scipy.constants as C

class SimUnits:
	"""Contains constants characterizing plasma units.
	Class is created with the unit density in mks units.

	t1 = unit of time in seconds

	x1 = unit of length in meters

	E1 = unit of electric field in volts per meter"""
	t1 = 0.0
	x1 = 0.0
	E1 = 0.0
	def __init__(self,n1):
		wp = np.sqrt(n1*C.e*C.e/(C.epsilon_0*C.m_e))
		self.t1 = 1/wp
		self.x1 = C.c/wp
		self.E1 = C.m_e * C.c * wp / C.e

class QuinticPulseData:
	"""Contains constants characterizing a quintic pulse with unit amplitude and base-to-base duration.

	energy = integration over square of quintic pulse

	FWHM = FWHM of square of quintic pulse"""
	energy = 0.3918
	FWHM = 0.3856

def ncrit(wavelength):
	"""Returns critical plasma density given a wavelength (mks units)"""
	w0 = 2*np.pi*C.c/wavelength
	return w0**2 * C.epsilon_0 * C.m_e / C.e**2

def courant(dx,dy,dz,c=1.0):
	"""Returns maximum time step allowed by the vacuum Courant condition.

	:param float dx: cell size in x
	:param float dy: cell size in y
	:param float dz: cell size in z
	:param float c: speed of light in the preferred units (defaults to plasma units)

	:return: maximum time step"""
	return 1/np.sqrt(1/dx/dx + 1/dy/dy + 1/dz/dz)/c

def getexp(a0,w0,r0,risetime,n1):
	"""Get experimental pulse parameters from simulation parameters.
	Inputs are in plasma units, outputs in mks.

	:param float a0: vector potential
	:param float w0: angular frequency
	:param float r0: 1/e of amplitude radius
	:param float risetime: quintic rise parameter
	:param float n1: unit density in *mks* units

	:return: (intensity,power,energy,FWHM)"""
	su = SimUnits(n1)
	qpd = QuinticPulseData()
	eta0 = np.sqrt(C.mu_0/C.epsilon_0)
	rmks = r0*su.x1
	base_to_base = 2.0*risetime*su.t1
	A0 = C.m_e * C.c * a0 / C.e
	wp = np.sqrt(n1*C.e*C.e/(C.epsilon_0*C.m_e))
	E0 = (w0/su.t1) * A0
	I0 = 0.5*E0*E0/eta0
	P0 = 0.5*np.pi*rmks*rmks*I0
	U0 = P0 * qpd.energy * base_to_base
	FWHM = qpd.FWHM * base_to_base
	return I0,P0,U0,FWHM

def getsim(U0,l0,r0,FWHM,n1):
	"""Get simulation pulse parameters from experimental parameters.
	Inputs are in mks, outputs in plasma units.

	:param float U0: pulse energy
	:param float l0: central wavelength
	:param float r0: 1/e of amplitude radius
	:param float FWHM: pulse length, full-width-half-maximum of intensity
	:param float n1: unit density

	:return: (vector potential,angular frequency,quintic risetime)"""
	su = SimUnits(n1)
	qpd = QuinticPulseData()
	eta0 = np.sqrt(C.mu_0/C.epsilon_0)
	w0_sim = su.t1 * 2.0*np.pi*C.c / l0
	r_sim = r0/su.x1
	base_to_base = FWHM / qpd.FWHM
	P0 = U0 / (qpd.energy * base_to_base)
	I0 = 2.0*P0 / (np.pi * r0 * r0)
	E0 = np.sqrt(2.0*I0*eta0)
	A0 = su.t1 * E0 / w0_sim
	a0 = C.e * A0 / (C.m_e * C.c)
	risetime = 0.5 * base_to_base / su.t1
	return a0,w0_sim,risetime
