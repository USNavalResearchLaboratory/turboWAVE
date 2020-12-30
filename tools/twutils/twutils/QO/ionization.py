'''
Module :samp:`twutils.QO.ionization`
====================================

This module can be used to compute ionization rates that may be needed in pre or post processing.  It contains dictionaries of atomic data, and classes for computing ionization rates based on various models appearing in the literature.  The ionization classes use atomic units, but include routines for converting to and from mks units.

References
----------

The ionization models herein are detailed in the following references.

#. L.V. Keldysh, Sov. Phys. JETP 20 (5), 1307 (1965)

#. A.M. Perelemov, V.S. Popov, M.V. Terent'ev, Sov. Phys. JETP 23 (5), 924 (1966)

#. M.V. Ammosov, N.B. Delone, V.P. Krainov, Sov. Phys. JETP 64 (6), 1191 (1986)

#. G.L. Yudin and M.Y. Ivanov, Phys. Rev. A 64, 013409 (2001)

#. S.V. Popruzhenko, V.D. Mur, V.S. Popov, D. Bauer, Phys. Rev. Lett. 101, 193003 (2008)

#. M\. Klaiber, E. Yakaboylu, K.Z. Hatsagortsyan, Phys. Rev. A 87, 023417 and 023418 (2013)

Constants
----------

* Dictionary ``ip``: keys are atomic symbols (e.g. 'H','He',etc.), values are lists of ionization potentials in eV.

* Dictionary ``atomic_mass``: keys are atomic symbols, values are atomic masses in proton mass units.

* Dictionary ``atomic_number`` : keys are atomic symbols, values are the number of protons.

An example of usage is as follows::

	import twutils.QO.ionization as iz

	print('argon ionization potential =',iz.ip['Ar'][0])
	print('argon[3+] ionization potential =',iz.ip['Ar'][3])
	print('atomic mass of argon is',iz.atomic_mass['Ar'])

Classes and functions
----------------------
'''

import numpy as np
from scipy import special as spsf
from scipy import constants as C
import scipy.integrate


# H,He,Li,Be,B,C,N,O,Ne,Ar data from photocopy of three pages (10-205--10-207) of unknown book
# Krypton I-XXVI, J. Sugar, A. Musgrove, J. Phys. Chem. Ref. Data 20, 859 (1991)
# Krypton V-XXVI, T. Shirai, K. Okazaki, and J. Sugar, J. Phys. Chem. Ref. Data 24, 1577 (1995)
# Xenon, E.B. Saloman, J. Phys. Chem. Ref. Data 33, 765 (2004)

ip = {'H' : [13.59844],
	'He' : [24.58741,54.41778],
	'Li' : [5.39172,75.64018,122.45429],
	'Be' : [9.32263,18.21116,153.89661,217.71865],
	'B' : [8.29803,25.15484,37.93064,259.37521,340.22580],
	'C' : [11.26030,24.38332,47.8878,64.4939,392.087,489.9934],
	'N' : [14.53414,29.6013,47.44924,77.4735,97.8902,552.0718,667.046],
	'O' : [13.61806,35.11730,54.9355,77.41353,113.8990,138.1197,739.29,871.4101],
	'Ne' : [21.56454,40.96328,63.45,97.12,126.21,157.93,207.2759,239.0989,1195.8286,1362.1995],
	'Ar' : [15.75962,27.62967,40.74,59.81,75.02,91.009,124.323,143.460,422.45,478.69,
		538.96,618.26,686.10,755.74,854.77,918.03,4120.8857,4426.2296],
	'Kr' : [13.99961,24.35985,36.950,52.5,64.7,78.5,111.0,125.802,230.85,268.2,
		308.2,350.1,390.9,446.6,491.8,540.7,591.5,641,785.9,833.0,
		883.9,936.7,997.7,1050.9,1151.4,1205.3,2928,3070,3227,3381,
		3594,3760,3966,4108.54,17297.16,17936.21],
	'Xe' : [12.130,20.975,31.05,40.9,54.14,66.703,91.6,105.976,179.85,201.7,
		229.02,263.5,294.4,325.3,358.3,389.6,420.9,452.2,572.5,607.7,
		642.9,678.1,726.0,762.4,852.7,857.0,1490,1491,1587,1684,
		1781,1877,1987,2085,2183,2281,2548,2637,2726,2814,
		3001,3093,3296,3334,7224,7491,7758,8024,8617,8899,
		9607,9813,40273,41300]}

atomic_mass = { 'H' : 1.00794 , 'He' : 4.002602 , 'Li' : 6.941 ,
	'Be' : 9.012182 , 'B' : 10.811 , 'C' : 12.0107 , 'N' : 14.0067 ,
	'O' : 15.9994 , 'F' : 18.9984 , 'Ne' : 20.1797 , 'Na' : 22.98977 ,
	'Mg' : 24.3050 , 'Al' : 26.9815 , 'Si' : 28.0855 , 'P' : 30.9738 ,
	'S' : 32.065 , 'Cl' : 35.453 , 'K' : 39.0983 , 'Ar' : 39.948 ,
	'Ca' : 40.078 , 'Sc' : 44.956 , 'Ti' : 47.867 , 'V' : 50.9415 ,
	'Cr' : 51.9961 , 'Mn' : 54.938 , 'Fe' : 55.845 , 'Co' : 58.933 ,
	'Ni' : 58.6934 , 'Cu' : 63.546 , 'Zn' : 65.409 , 'Ga' : 69.723 ,
	'Ge' : 72.64 , 'As' : 74.92160 , 'Se' : 78.96 , 'Br' : 79.904 ,
	'Kr' : 83.798 , 'Rb' : 85.4678 , 'Sr' : 87.62 , 'Y' : 88.90585 ,
	'Zr' : 91.224 , 'Nb' : 92.906 , 'Mo' : 95.94 , 'Tc' : 98 ,
	'Ru' : 101.07 , 'Rh' : 102.905 , 'Pd' : 106.42 , 'Ag' : 107.8682 ,
	'Cd' : 112.411 , 'In' : 114.818 , 'Sn' : 118.710 , 'Sb' : 121.760 ,
	'Te' : 127.60 , 'I' : 126.904 , 'Xe' : 131.293 }

atomic_number = { 'H' : 1 , 'He' : 2 , 'Li' : 3 ,
	'Be' : 4 , 'B' : 5 , 'C' : 6 , 'N' : 7 ,
	'O' : 8 , 'F' : 9 , 'Ne' : 10 , 'Na' : 11 ,
	'Mg' : 12 , 'Al' : 13 , 'Si' : 14 , 'P' : 15 ,
	'S' : 16 , 'Cl' : 17 , 'K' : 19 , 'Ar' : 18 , # n.b. ordered by mass
	'Ca' : 20 , 'Sc' : 21 , 'Ti' : 22 , 'V' : 23 ,
	'Cr' : 24 , 'Mn' : 25 , 'Fe' : 26 , 'Co' : 27 ,
	'Ni' : 28 , 'Cu' : 29 , 'Zn' : 30 , 'Ga' : 31 ,
	'Ge' : 32 , 'As' : 33 , 'Se' : 34 , 'Br' : 35 ,
	'Kr' : 36 , 'Rb' : 37 , 'Sr' : 38 , 'Y' : 39 ,
	'Zr' : 40 , 'Nb' : 41 , 'Mo' : 42 , 'Tc' : 43 ,
	'Ru' : 44 , 'Rh' : 45 , 'Pd' : 46 , 'Ag' : 47 ,
	'Cd' : 48 , 'In' : 49 , 'Sn' : 50 , 'Sb' : 51 ,
	'Te' : 52 , 'I' : 53 , 'Xe' : 54 }

def Shifted_Uion_H(E):
	"""Returns the Stark shifted ionization energy of Hydrogen"""
	return 0.5 + (9.0/4.0)*E**2 + 55.55*E**4 + 4908.0*E**6 + 794200*E**8

def Keldysh_Parameter(omega,Uion,E):
	'''Compute the adiabaticity parameter (Keldysh parameter)

	:param double omega: radiation frequency in atomic units
	:param double Uion: ionization potential in atomic units
	:param double E: electric field in atomic units

	:return: dimensionless adiabaticity parameter'''
	return omega*np.sqrt(2.0*Uion)/E

def wfunc(x):
	"""Exponential integral from PPT theory (Dawson's integral)"""
	return np.real(0.5 * np.exp(-x**2) * np.sqrt(np.pi) * spsf.erf(x*1j) / 1j)

def CycleAveragingFactor(Uion,E):
	'''Returns the PPT cycle averaging factor, multiply static rate by this factor to get average rate.

	:param double Uion: ionization potential in atomic units
	:param double E: electric field in atomic units'''
	return np.sqrt(3.0/np.pi) * np.sqrt(E) / (2*Uion)**0.75

class AtomicUnits:
	'''Class for converting to and from atomic units. Ionization classes derive from this.'''
	def __init__(self):
		'''Create an AtomicUnits object'''
		# atomic units in mks units
		self.xa = C.hbar / (C.m_e*C.c*C.alpha)
		self.ta = C.hbar / (C.m_e*C.c**2*C.alpha**2)
		self.ua = C.m_e*C.c**2*C.alpha**2
		self.ea = self.ua/self.xa/C.e
	def energy_to_eV(self,u):
		'''Convert an energy from atomic units to eV'''
		return u*self.ua/C.e
	def length_to_mks(self,x):
		'''Convert a length from atomic to mks units'''
		return x*self.xa
	def rate_to_mks(self,rate):
		'''Convert a time from atomic to mks units'''
		return rate/self.ta
	def time_to_mks(self,t):
		'''Convert a time from atomic to mks units'''
		return t*self.ta
	def field_to_mks(self,field):
		'''Convert electric field from atomic to mks units'''
		return field*self.ea
	def energy_to_au(self,u):
		'''Convert an energy in eV to atomic units'''
		return u*C.e/self.ua
	def length_to_au(self,x):
		'''Convert a length from mks to atomic units'''
		return x/self.xa
	def rate_to_au(self,rate):
		'''Convert a time from mks to atomic units'''
		return rate*self.ta
	def time_to_au(self,t):
		'''Convert a time from mks to atomic units'''
		return t/self.ta
	def field_to_au(self,field):
		'''Convert electric field from mks to atomic units'''
		return field/self.ea

class Ionization(AtomicUnits):
	'''Base class for computing ionization rates'''
	def __init__(self,w0,Uion,Z,terms=1):
		'''Create an ionization object.  Base class should be treated as abstract.

		:param double w0: radiation angular frequency in atomic units
		:param double Uion: ionization potential in atomic units
		:param double Z: residual charge in atomic units
		:param int terms: terms to keep in PPT expansion'''
		super().__init__()
		self.w0 = w0
		self.Uion = Uion
		self.Z = Z
		self.terms = terms
		self.cutoff_field = 1e-3
	def ExtractEikonalForm(self,E,dt,w00=0.0,bandwidth=1.0):
		"""Extract amplitude, phase, and center frequency from a carrier resolved field E.
		The assumed form is :math:`E = \Re\{Ae^{i\phi}\}`.

		:param numpy.array E: Electric field, any shape, axis 0 is time.
		:param double dt: time step.
		:param double w00: Center frequency, if zero deduce from the data (default).
		:param double bandwidth: relative bandwidth to keep, if unity no filtering (default)."""
		Nt = E.shape[0]
		ndims = len(E.shape)
		dw = 2*np.pi/(Nt*dt)
		# Get the complexified spectrum and carrier frequency
		Ew = np.fft.fft(E,axis=0)
		Ew[np.int(Nt/2):,...] = 0.0 # eliminate negative frequencies before extracting carrier
		if w00==0.0:
			# Compute a suitable carrier frequency based on intensity weighting
			if ndims==1:
				indices = np.arange(0,Nt)
			else:
				indices = np.einsum('i,jk',np.arange(0,Nt),np.ones(Ew.shape[1:]))
			carrier_idx = np.int(np.average(indices,weights=np.abs(Ew)**2))
			w0 = carrier_idx*dw
		else:
			# Carrier frequency is set by caller
			w0 = w00
			carrier_idx = np.int(w0/dw)
		# Bandpass filtering
		if bandwidth!=1.0:
			low_idx = np.int((w0-0.5*bandwidth*w0)/dw)
			high_idx = np.int((w0+0.5*bandwidth*w0)/dw)
			Ew[:low_idx] = 0.0
			Ew[high_idx:] = 0.0
		# Form the complex envelope in time
		Ew = np.roll(Ew,-carrier_idx,axis=0)
		Et = 2*np.fft.ifft(Ew,axis=0)
		# Get amplitude and phase
		amp = np.abs(Et)
		phase = (np.angle(Et).swapaxes(0,ndims-1) + w0*np.linspace(0,(Nt-1)*dt,Nt)).swapaxes(0,ndims-1)
		return amp,phase,w0
	def Coeff(self,averaging):
		# Precompute coefficients for ionization formulas
		return 0.0,1.0,-1.0
	def Rate(self,Es,averaging):
		'''Get the ionization rate, returned value is in atomic units.

		:param np.array Es: Time domain electric field in atomic units, any shape.
		:param bool averaging: Interpret Es a complex envelope and get cycle averaged rate.  In this case :math:`|Es|` is expected to give the crest.'''
		E = np.abs(Es) + self.cutoff_field
		C1,C2,C3 = self.Coeff(averaging)
		return C1*E**C2*np.exp(C3/E)

class Hydrogen(Ionization):
	"""Landau & Lifshitz formula for hydrogen, using PPT averaging"""
	def Coeff(self,averaging):
		C_pre = 4.0
		C_pow = -1.0
		C_exp = -2.0/3.0
		if averaging:
			C_pre *= np.sqrt(3.0/np.pi) / (2*self.Uion)**0.75
			C_pow += 0.5
		return C_pre,C_pow,C_exp

class ADK(Ionization):
	'''Standard ADK tunneling for s orbitals'''
	def Coeff(self,averaging):
		Uion = self.Uion
		Z = self.Z
		nstar = Z/np.sqrt(2*Uion)
		D = ((4.0*np.exp(1.0)*Z**3)/nstar**4)**nstar
		C_pre = D*D/(8.0*np.pi*Z)
		C_pow = 1.0-2*nstar
		C_exp = -2*Z**3/(3.0*nstar**3)
		if averaging:
			C_pre *= np.sqrt(3.0/np.pi) / (2*Uion)**0.75
			C_pow += 0.5
		return C_pre,C_pow,C_exp

class Klaiber(Ionization):
	'''Klaiber dressed coulomb corrected SFA (Eq. 97)'''
	def Coeff(self,averaging):
		Uion = self.Uion
		Ua2 = Uion*C.alpha**2
		Ea = (2*Uion)**1.5 # definition of Ea applies to Eq. 97 but not Fig. 4 caption
		C_pre = 2**(3-4*Ua2)*np.sqrt(3/np.pi)*(1-7*Ua2/72)*np.exp(4*Ua2)*(2*Uion)**(1.75-3*Ua2)
		C_pre /= spsf.gamma(3-2*Ua2)
		C_pow = 2*Ua2 - 0.5
		C_exp = -(2.0/3.0)*Ea*(1-Ua2/12)
		if not averaging:
			C_pre /= np.sqrt(3.0/np.pi) / (2*Uion)**0.75
			C_pow -= 0.5
		return C_pre,C_pow,C_exp

class PPT_Tunneling(ADK):
	'''Evaluating the PPT theory in the tunneling limit, with the infinite sum solved analytically.'''
	def Coeff(self,averaging):
		# Derive from ADK so that in the tunneling limit,
		# we can simply use Popov's note to undo application of the Stirling formula
		nstar = self.Z/np.sqrt(2*self.Uion)
		NPPT = (2**(2*nstar-1)/spsf.gamma(nstar+1))**2
		NADK = (1/(8*np.pi*nstar))*(4*np.exp(1)/nstar)**(2*nstar)
		C1,C2,C3 = super().Coeff(averaging)
		C1 *= NPPT/NADK
		return C1,C2,C3

class PPT(Ionization):
	'''Full PPT theory accounting for multiphoton and tunneling limits.
	If instantaneous rate is requested, the cycle averaging factor is divided out.
	However, it may be better to use the ``PPT_Tunneling`` class in this case.'''
	def Rate(self,Es,averaging):
		E = np.abs(Es) + self.cutoff_field
		w = self.w0
		F0 = (2*self.Uion)**1.5
		nstar = self.Z/np.sqrt(2*self.Uion)
		lstar = nstar - 1
		C2 = 2**(2*nstar) / (nstar*spsf.gamma(nstar+lstar+1)*spsf.gamma(nstar-lstar))
		gam = np.sqrt(2.0*self.Uion)*w/E
		alpha = 2 * (np.arcsinh(gam)-gam/np.sqrt(1+gam**2))
		beta = 2*gam/np.sqrt(1+gam**2)
		g = (3/(2*gam))*((1+1/(2*gam**2))*np.arcsinh(gam)-np.sqrt(1+gam**2)/(2*gam))
		nu = (self.Uion/w) * (1 + 1/(2*gam**2))
		A0 = 0.0*Es # if Es is an array we want A0 to be one also
		dnu = np.ceil(nu) - nu
		for n in range(self.terms):
			A0 += np.exp(-alpha*(n+dnu))*wfunc(np.sqrt(beta*(n+dnu)))
		A0 *= (4/np.sqrt(3*np.pi)) * (gam**2/(1+gam**2))
		ans = A0*(E*np.sqrt(1+gam**2)/(2*F0))**1.5
		ans *= self.Uion*C2*np.sqrt(6/np.pi) * np.exp(-2.0*F0*g/(3*E))
		# Following is the Coulomb correction factor
		ans *= (2*F0/E)**(2*nstar)
		if averaging:
			return ans
		else:
			# If explicit, simply unwind cycle averaging.
			# This does not account for changes in sub-cycle structure with adiabaticity, but gives right average behavior.
			return ans/CycleAveragingFactor(self.Uion,E)

class PMPB(Ionization):
	'''PMBP theory, similar to PPT but with different Coulomb factor and bound state coefficient.
	This routine takes lstar=0, unlike others where lstar=nstar-1.'''
	def Rate(self,Es,averaging):
		# Mappings from our PPT notation to PMPB notation:
		# beta -> beta , gam -> gam , Uion -> I , Uion/w -> K0
		# E/F0 -> F , alpha -> 2*c1 , nu -> nth , g -> g , wfunc -> script-F
		E = np.abs(Es) + self.cutoff_field
		w = self.w0
		F = E / (2*self.Uion)**1.5
		nstar = self.Z/np.sqrt(2*self.Uion)
		C2 = 2**(2*nstar-2.0) / spsf.gamma(nstar+1)**2
		gam = np.sqrt(2.0*self.Uion)*w/E
		c1 = np.arcsinh(gam)-gam/np.sqrt(1+gam**2)
		beta = 2*gam/np.sqrt(1+gam**2)
		g = (3/(2*gam))*((1+1/(2*gam**2))*np.arcsinh(gam)-np.sqrt(1+gam**2)/(2*gam))
		nu = (self.Uion/w) * (1 + 1/(2*gam**2))
		A0 = 0.0*Es # if Es is an array we want A0 to be one also
		dnu = np.ceil(nu) - nu
		for n in range(self.terms):
			A0 += np.exp(-2*c1*(n+dnu))*wfunc(np.sqrt(beta*(n+dnu)))
		ans = A0 * (2/np.pi) * C2 * self.Uion * (self.Uion/w)**(-1.5) * np.sqrt(beta)
		ans *= np.exp(-2*g/(3*F))
		# Following is the improved Coulomb correction factor
		ans *= (2/F)**(2*nstar) * (1+gam*2/np.exp(1))**(-2*nstar)
		if averaging:
			return ans
		else:
			# If explicit, simply unwind cycle averaging.
			# This does not account for changes in sub-cycle structure with adiabaticity, but gives right average behavior.
			return ans/CycleAveragingFactor(self.Uion,E)

class YI(Ionization):
	"""Yudin-Ivanov phase dependent ionization rate."""
	def Rate(self,amp,phase):
		"""Get the phase dependent rate.  The Phase can be extracted using the
		Ionization class method ExtractEikonalForm.

		:param numpy.array amp: Amplitude of the field, any shape
		:param numpy.array phase: Phase of the field, any shape"""
		w = self.w0
		amp = amp + self.cutoff_field
		theta = (phase - 0.5*np.pi)%np.pi - 0.5*np.pi
		nstar = self.Z/np.sqrt(2*self.Uion)
		lstar = nstar - 1
		Anl = 2**(2*nstar) / (nstar*spsf.gamma(nstar+lstar+1)*spsf.gamma(nstar-lstar))
		gam = np.sqrt(2.0*self.Uion)*w/amp
		a = 1+gam*gam-np.sin(theta)**2
		b = np.sqrt(a*a+4*gam*gam*np.sin(theta)**2)
		c = np.sqrt((np.sqrt((b+a)/2)+gam)**2 + (np.sqrt((b-a)/2)+np.sin(np.abs(theta)))**2)
		Phi = (gam**2 + np.sin(theta)**2 + 0.5)*np.log(c)
		Phi -= 3*(np.sqrt(b-a)/(2*np.sqrt(2)))*np.sin(np.abs(theta))
		Phi -= (np.sqrt(b+a)/(2*np.sqrt(2)))*gam
		kappa = np.log(gam+np.sqrt(gam**2+1)) - gam/np.sqrt(1+gam**2)
		alpha = 2 * (np.arcsinh(gam)-gam/np.sqrt(1+gam**2))
		beta = 2*gam/np.sqrt(1+gam**2)
		nu = (self.Uion/w) * (1 + 1/(2*gam**2))
		A0 = 0.0*amp
		dnu = np.ceil(nu) - nu
		for n in range(self.terms):
			A0 += np.exp(-alpha*(n+dnu))*wfunc(np.sqrt(beta*(n+dnu)))
		A0 *= (4/np.sqrt(3*np.pi)) * (gam**2/(1+gam**2))
		pre = Anl * np.sqrt(3*kappa/gam**3)*(1+gam**2)**0.75 * A0 * self.Uion
		pre *= (2*(2*self.Uion)**1.5 / amp)**(2*nstar-1)
		return pre * np.exp(-amp**2 * Phi / w**3)
