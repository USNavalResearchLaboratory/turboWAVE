"""Functions to help set up turboWAVE input files.
Contents:
class SimUnits
class QuinticPulseData
def courant
def getexp
def getsim"""

import numpy as np
import scipy.constants as C

class SimUnits:
	"""Contains simulation units t1,x1,E1 expressed in mks
	Pass unit density in cm^-3 upon instantiation"""
	t1 = 0.0
	x1 = 0.0
	E1 = 0.0
	def __init__(self,n1):
		wp = np.sqrt(n1*1e6*C.e*C.e/(C.epsilon_0*C.m_e))
		self.t1 = 1/wp
		self.x1 = C.c/wp
		self.E1 = C.m_e * C.c * wp / C.e

class QuinticPulseData:
	"""energy = integration over square of quintic pulse
	FWHM = FWHM of square of quintic pulse
	Parameters are for base to base length unity, and amplitude unity"""
	energy = 0.3918
	FWHM = 0.3856

def courant(dx,dy,dz):
	"""Returns maximum time step
	inputs:
	dx , dy , dz = cell sizes"""
	return 1/np.sqrt(1/dx/dx + 1/dy/dy + 1/dz/dz)
	
def getexp(a0,w0,r0,risetime,n1):
	"""Returns intensity (W/cm^2), Power (TW), energy (J), FWHM (fs)
	inputs:
	a0 = normalized vector potential
	w0 = angular frequency in wp^-1
	r0 = 1/e of field in c/wp
	risetime = quintic rise parameter in wp^-1
	n1 = unit density in cm^-3"""
	su = SimUnits(n1)
	qpd = QuinticPulseData()
	rmks = r0*su.x1
	base_to_base = 2.0*risetime*su.t1
	A0 = C.m_e * C.c * a0 / C.e
	wp = np.sqrt(n1*1e6*C.e*C.e/(C.epsilon_0*C.m_e))
	E0 = (w0/su.t1) * A0
	I0 = 0.5*E0*E0/377.0
	P0 = 0.5*np.pi*rmks*rmks*I0
	U0 = P0 * qpd.energy * base_to_base
	FWHM = qpd.FWHM * base_to_base
	print('I0[W/cm^2] = {0:.3g}'.format(I0*1e-4))
	print('P0[TW] = {0:.3g}'.format(P0*1e-12))
	print('U0[J] = {0:.3g}'.format(U0))
	print('FWHM[fs] = {0:.3g}'.format(FWHM*1e15))
	return I0*1e-4,P0*1e-12,U0,FWHM*1e15
	
def getsim(U0,l0,r0,FWHM,n1):
	"""Returns a, w, risetime
	inputs:
	U0 = pulse energy (J)
	l0 = wavelength (um)
	r0 = 1/e of field (um)
	FWHM = pulse length (fs)
	n1 = unit density in cm^-3"""
	su = SimUnits(n1)
	qpd = QuinticPulseData()
	w0_sim = su.t1 * 2.0*np.pi*C.c / (l0*1e-6)
	r_sim = r0*1e-6/su.x1
	base_to_base = FWHM*1e-15 / qpd.FWHM
	P0 = U0 / (qpd.energy * base_to_base)
	I0 = 2.0*P0 / (np.pi * r0 * r0 * 1e-12)
	E0 = np.sqrt(2.0*I0*377.0)
	A0 = su.t1 * E0 / w0_sim
	a0 = C.e * A0 / (C.m_e * C.c)
	risetime = 0.5 * base_to_base / su.t1
	print('a = {0:.3g}'.format(a0))
	print('w = {0:.3g}'.format(w0_sim))
	print('risetime = {0:.3g}'.format(risetime))
	return a0,w0_sim,risetime
