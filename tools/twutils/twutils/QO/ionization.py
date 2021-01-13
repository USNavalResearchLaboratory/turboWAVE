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
	def __init__(self,averaging,Uion,Z,lstar=0,l=0,m=0,w0=0,terms=1):
		'''Create an ionization object.  Base class should be treated as abstract.

		:param bool averaging: whether this object will cycle average the rate, or try to return static rate
		:param double Uion: ionization potential in atomic units
		:param double Z: residual charge in atomic units
		:parameter float lstar: effective orbital angular momentum
		:parameter int l: orbital angular momentum quantum number
		:parameter int m: orbital angular momentum projection in direction of field
		:parameter float w0: angular frequency of radiation in atomic units
		:parameter int terms: terms to keep in expansion over multiphoton orders'''
		super().__init__()
		self.avg = averaging
		self.Uion = Uion
		self.Z = Z
		self.lstar = lstar
		self.l = l
		self.m = m
		self.w0 = w0
		self.terms = terms
		self.cutoff_field = 1e-3
		# coefficients for standard form ionization formulas
		self.C_pre = 0.0
		self.C_pow = 0.0
		self.C_exp = 0.0
	def ToAverage(self):
		'''Switch to averaging.  This must not be called more than once.'''
		self.C_pre *= np.sqrt(3.0/np.pi) / (2*self.Uion)**0.75
		self.C_pow += 0.5
	def ToStatic(self):
		'''Switch to static rate.  This must not be called more than once.'''
		self.C_pre *= (2*self.Uion)**0.75 / np.sqrt(3.0/np.pi)
		self.C_pow -= 0.5
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
	def Rate(self,Es):
		'''Get the ionization rate, returned value is in atomic units.
		The object's internal state decides whether this is a static or cycle averaged rate.

		:param np.array Es: Time domain electric field in atomic units, any shape.  If a complex envelope is given :math:`|Es|` is expected to give the crest.'''
		E = np.abs(Es) + self.cutoff_field
		return self.C_pre * E**self.C_pow * np.exp(self.C_exp/E)

class Hydrogen(Ionization):
	"""Landau & Lifshitz formula for hydrogen, using PPT averaging"""
	def __init__(self,averaging,Uion,Z):
		super().__init__(averaging,Uion,Z)
		self.C_pre = 4.0
		self.C_pow = -1.0
		self.C_exp = -2.0/3.0
		if averaging:
			self.ToAverage()

class ADK(Ionization):
	'''Full ADK for arbitrary states'''
	def __init__(self,averaging,Uion,Z,lstar=0,l=0,m=0):
		super().__init__(averaging,Uion,Z,lstar,l,m)
		nstar = Z/np.sqrt(2*Uion)
		m = abs(self.m)
		e = np.exp(1.0)
		sn2l2 = np.sqrt(nstar**2-lstar**2)
		self.C_pre = np.sqrt(3/np.pi**3) * (2*l+1) * spsf.factorial(l+m) / spsf.factorial(m) / spsf.factorial(l-m)
		self.C_pre *= (e/sn2l2)**(m+1.5) # exponent is of poor print quality in JETP
		self.C_pre *= ((nstar+lstar)/(nstar-lstar))**(lstar+0.5)
		self.C_pre *= (Z**2/nstar**3)
		self.C_pre *= (4*e*Z**3/nstar**3/sn2l2)**(2*nstar-m-1.5)
		self.C_pow = m+1.5-2*nstar
		self.C_exp = -2*Z**3/(3*nstar**3) # ADK 1986 has nstar**4 in Eq. 21, must be a typo?
		if not averaging:
			self.ToStatic()

class ADK00(Ionization):
	'''ADK theory simplified for lstar=l=m=0'''
	def __init__(self,averaging,Uion,Z):
		super().__init__(averaging,Uion,Z)
		nstar = Z/np.sqrt(2*Uion)
		D = ((4.0*np.exp(1.0)*Z**3)/nstar**4)**nstar
		self.C_pre = D*D/(8.0*np.pi*Z)
		self.C_pow = 1.0-2*nstar # N.b. this is already for the static rate
		self.C_exp = -2*Z**3/(3.0*nstar**3)
		if averaging:
			self.ToAverage()

class KYH(Ionization):
	'''KYH theory for relativistic tunneling, using the dressed coulomb corrected SFA (Eq. 97)'''
	def __init__(self,averaging,Uion,Z):
		super().__init__(averaging,Uion,Z)
		Ua2 = Uion*C.alpha**2
		Ea = (2*Uion)**1.5 # definition of Ea applies to Eq. 97 but not Fig. 4 caption
		self.C_pre = 2**(3-4*Ua2)*np.sqrt(3/np.pi)*(1-7*Ua2/72)*np.exp(4*Ua2)*(2*Uion)**(1.75-3*Ua2)
		self.C_pre /= spsf.gamma(3-2*Ua2)
		self.C_pow = 2*Ua2 - 0.5
		self.C_exp = -(2.0/3.0)*Ea*(1-Ua2/12)
		if not averaging:
			self.ToStatic()

class PPT_Tunneling(Ionization):
	'''Evaluating the PPT theory in the tunneling limit, with the infinite sum solved analytically.
	This is PPT Eq. 59 multiplied by the standard Coulomb factor, and taking gamma=0.
	For the bound state coefficient we use ADK Eq. 19.'''
	def __init__(self,averaging,Uion,Z,lstar=0,l=0,m=0):
		super().__init__(averaging,Uion,Z,lstar,l,m)
		nstar = Z/np.sqrt(2*Uion)
		m = abs(self.m)
		F0 = Z**3/nstar**3
		# First without the Coulomb factor
		self.C_pre = Uion
		self.C_pre *= 2**(2*nstar) / (nstar*spsf.gamma(nstar+lstar+1)*spsf.gamma(nstar-lstar)) # |C|^2
		self.C_pre *= np.sqrt(6/np.pi)
		self.C_pre *= (2*l+1) * spsf.factorial(l+m) / 2**m / spsf.factorial(m) / spsf.factorial(l-m)
		self.C_pre *= (0.5/F0)**(m+1.5)
		self.C_pow = m+1.5
		self.C_exp = -2*F0/3
		# Account for Coulomb correction
		self.C_pre *= (2*F0)**(2*nstar)
		self.C_pow -= 2*nstar
		if not averaging:
			self.ToStatic()

class HybridTheories(Ionization):
	'''Base class for theories that handle both multiphoton and tunneling limits'''
	def __init__(self,averaging,Uion,Z,lstar=0,l=0,m=0,w0=0,terms=1):
		if m!=0.0:
			raise ValueError('Hybrid theories all require m=0')
		if w0==0.0:
			raise ValueError('Hybrid theories require w>0')
		super().__init__(averaging,Uion,Z,lstar,l,m,w0,terms)
	def Params(self,E,w,Uion,Z):
		nstar = Z/np.sqrt(2*Uion)
		F = E*nstar**3/Z**3
		gam = Z*w/(nstar*E)
		g = (3/(2*gam))*((1+1/(2*gam**2))*np.arcsinh(gam)-np.sqrt(1+gam**2)/(2*gam))
		nu = (self.Uion/w) * (1 + 1/(2*gam**2))
		return nstar,F,gam,g,nu
	def FourierSum(self,gam,nu):
		A0 = 0.0*gam; # make this the same shape as gam, if gam is an array
		alpha = 2 * (np.arcsinh(gam)-gam/np.sqrt(1+gam**2))
		beta = 2*gam/np.sqrt(1+gam**2)
		dnu = np.ceil(nu) - nu
		for n in range(self.terms):
			A0 += np.exp(-alpha*(n+dnu))*spsf.dawsn(np.sqrt(beta*(n+dnu))) # only valid for m=0
		return A0;

class PPT(HybridTheories):
	'''PPT theory accounting for multiphoton and tunneling limits, but assuming m=0.
	This is PPT Eq. 54 multiplied by the standard Coulomb factor, and taking m=0.
	For the bound state coefficient we use ADK Eq. 19.
	If static rate is requested, the result is only valid in the limit gamma=0'''
	def Rate(self,Es):
		E = np.abs(Es) + self.cutoff_field
		w = self.w0
		lstar = self.lstar
		l = self.l
		nstar,F,gam,g,nu = self.Params(E,w,self.Uion,self.Z)
		ans = self.Uion
		ans *= 2**(2*nstar) / (nstar*spsf.gamma(nstar+lstar+1)*spsf.gamma(nstar-lstar)) # |C|^2
		ans *= np.sqrt(6/np.pi) * (2*l+1)
		ans *= (F*np.sqrt(1+gam**2)/2)**1.5
		ans *= (4/np.sqrt(3*np.pi)) * (gam**2/(1+gam**2)) * self.FourierSum(gam,nu)
		ans *= np.exp(-2*g/(3*F))
		# Following is the Coulomb correction factor
		ans *= (2/F)**(2*nstar)
		if self.avg:
			return ans
		else:
			# If explicit, simply unwind cycle averaging.
			# This does not account for changes in sub-cycle structure with adiabaticity, but gives right average behavior.
			return ans/CycleAveragingFactor(self.Uion,E)

class PMPB(HybridTheories):
	'''PMBP theory, similar to PPT but with improved Coulomb factor, assuming m=0.
	Take lstar=0,l=0,m=0 to strictly recover PMPB Eq. 6.'''
	def Rate(self,Es):
		# Note that |C_PPT|^2 = 4*|C_PMPB|^2
		# This can be seen in how the respective authors write the wavefunction.
		E = np.abs(Es) + self.cutoff_field
		w = self.w0
		lstar = self.lstar
		l = self.l
		nstar,F,gam,g,nu = self.Params(E,w,self.Uion,self.Z)
		ans = (2/np.pi)*self.Uion
		ans *= 2**(2*nstar-2) / (nstar*spsf.gamma(nstar+lstar+1)*spsf.gamma(nstar-lstar)) # |C_PMBP|^2 generalized
		ans *= 2*l + 1 # generalize also for l>0, but demand m=0
		ans *= (self.Uion/w)**-1.5
		ans *= np.sqrt(2*gam/np.sqrt(1+gam**2))
		ans *= self.FourierSum(gam,nu)
		ans *= np.exp(-2*g/(3*F))
		# Following is the improved Coulomb correction factor
		ans *= (2/F)**(2*nstar) * (1+gam*2/np.exp(1))**(-2*nstar)
		if self.avg:
			return ans
		else:
			# If explicit, simply unwind cycle averaging.
			# This does not account for changes in sub-cycle structure with adiabaticity, but gives right average behavior.
			return ans/CycleAveragingFactor(self.Uion,E)

class YI(HybridTheories):
	"""Yudin-Ivanov phase dependent ionization rate, Eq. 17., assuming m=0."""
	def Rate(self,amp,phase):
		"""Get the phase dependent rate.  The Phase can be extracted using the
		Ionization class method ExtractEikonalForm.

		:param numpy.array amp: Amplitude of the field, any shape
		:param numpy.array phase: Phase of the field, any shape"""
		amp = amp + self.cutoff_field
		theta = (phase - 0.5*np.pi)%np.pi - 0.5*np.pi
		w = self.w0
		lstar = self.lstar
		l = self.l
		nstar,F,gam,g,nu = self.Params(amp,w,self.Uion,self.Z)
		kappa = np.log(gam+np.sqrt(gam**2+1)) - gam/np.sqrt(1+gam**2)
		pre = self.Uion
		pre *= 2**(2*nstar) / (nstar*spsf.gamma(nstar+lstar+1)*spsf.gamma(nstar-lstar)) # A_nl = |C_PPT|^2
		pre *= (2*l + 1) * np.sqrt(3*kappa/gam**3) # here 2*l+1 = B_lm with m=0
		pre *= (1+gam**2)**0.75 * (4/np.sqrt(3*np.pi)) * (gam**2/(1+gam**2)) * self.FourierSum(gam,nu) # C
		pre *= (2/F)**(2*nstar-1)
		a = 1+gam*gam-np.sin(theta)**2
		b = np.sqrt(a*a+4*gam*gam*np.sin(theta)**2)
		c = np.sqrt((np.sqrt((b+a)/2)+gam)**2 + (np.sqrt((b-a)/2)+np.sin(np.abs(theta)))**2)
		Phi = (gam**2 + np.sin(theta)**2 + 0.5)*np.log(c)
		Phi -= 3*(np.sqrt(b-a)/(2*np.sqrt(2)))*np.sin(np.abs(theta))
		Phi -= (np.sqrt(b+a)/(2*np.sqrt(2)))*gam
		return pre * np.exp(-amp**2 * Phi / w**3)
