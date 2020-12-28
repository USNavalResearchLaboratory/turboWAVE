'''
Module :samp:`ionization`
-------------------------

Contains dictionary of ionization potentials ``ip`` and classes for computing
ionization rates.
'''


#[1] A.M. Perelemov, V.S. Popov, M.V. Terent'ev, Sov. Phys. JETP 23 (5), 924 (1966)
#[2] M.V. Ammosov, N.B. Delone, V.P. Krainov, Sov. Phys. JETP 64 (6), 1191 (1986)
#[3] L.V. Keldysh, Sov. Phys. JETP 20 (5), 1307 (1965)"""

import numpy as np
from scipy import special as spsf
from scipy import constants as C
import scipy.integrate


# data from photocopy of three pages (10-205--10-207) of unknown book
# refs listed in photocopy:
# Moore, C.E., Natl. Stand. Ref. Data Ser.--Natl. Bur. Stand. (U.S.) No. 34, 1970
# Martin, W.C. et al., Natl. Stand. Ref. Data Ser.---Natl. Bur. Stand. (U.S.) No. 60, 1978
# Sugar, J. and Corliss, C., J. Phys. Chem. Ref. Data Vol. 14, Suppl. 2, 1985
# Refs to papers in J. Phys. Chem. Ref. Data (1973-1991) may be found in cumulative index to that journal
# Cohen, E.R. and Taylor, B.N., J. Phys. Chem. Ref. Data 17, 1795 (1988)

ip = {'H' : [13.59844],
	'He' : [24.58741,54.41778],
	'Li' : [5.39172,75.64018,122.45429],
	'Be' : [9.32263,18.21116,153.89661,217.71865],
	'B' : [8.29803,25.15484,37.93064,259.37521,340.22580],
	'C' : [11.26030,24.38332,47.8878,64.4939,392.087,489.9934],
	'N' : [14.53414,29.6013,47.44924,77.4735,97.8902,552.0718,667.046],
	'O' : [13.61806,35.11730,54.9355,77.41353,113.8990,138.1197,739.29,871.4101],
	'Ne' : [21.56454,40.96328,63.45,97.12,126.21,157.93,207.2759,239.0989,1195.8286,1362.1995],
	'Ar' : [15.75962,27.62967,40.74,59.81,75.02,91.009,124.323,143.460,422.45,478.69,538.96,618.26,686.10,755.74,854.77,918.03,4120.8857,4426.2296],
	'Kr' : [13.99961,24.35985]}

def Shifted_Uion_H(E):
	"""Returns the Stark shifted ionization energy of Hydrogen"""
	return 0.5 + (9.0/4.0)*E**2 + 55.55*E**4 + 4908.0*E**6 + 794200*E**8

def Keldysh_Parameter(omega,Uion,E):
	'''Returns the Keldysh parameter

	:param double omega: radiation frequency in atomic units
	:param double Uion: ionization potential in atomic units
	:param double E: electric field in atomic units'''
	return omega*np.sqrt(2.0*Uion)/E

def wfunc(x):
	"""Exponential integral from PPT theory"""
	return np.real(0.5 * np.exp(-x**2) * np.sqrt(np.pi) * spsf.erf(x*1j) / 1j)

def CycleAveragingFactor(Uion,E):
	'''Returns the PPT cycle averaging factor, multiply static rate by this factor to get average rate.

	:param double Uion: ionization potential in atomic units
	:param double E: electric field in atomic units'''
	return np.sqrt(3.0/np.pi) * np.sqrt(E) / (2*Uion)**0.75

class AtomicUnits:
	def __init__(self,mks_length):
		# simulation units in mks units
		self.xs = mks_length
		self.ts = mks_length/C.c
		self.es = C.m_e*C.c/(self.ts*C.e)
		# atomic units in mks units
		self.xa = C.hbar / (C.m_e*C.c*C.alpha)
		self.ta = C.hbar / (C.m_e*C.c**2*C.alpha**2)
		self.ua = C.m_e*C.c**2*C.alpha**2
		self.ea = self.ua/self.xa/C.e
		# unit of density in mks
		w0 = C.c/mks_length
		self.ncrit = w0**2 * C.epsilon_0 * C.m_e / C.e**2
	def rate_sim_to_au(self,rate):
		return rate*self.ta/self.ts
	def rate_au_to_sim(self,rate):
		return rate*self.ts/self.ta
	def field_sim_to_au(self,field):
		return field*self.es/self.ea
	def field_au_to_sim(self,field):
		return field*self.ea/self.es

class Ionization(AtomicUnits):
	def __init__(self,mks_length,w0,Uion,Z,terms=1):
		'''Create an ionization object.  Base class should be treated as abstract.

		:param double mks_length: normalizing length in mks units
		:param double w0: radiation angular frequency in atomic units
		:param double Uion: ionization potential in atomic units
		:param double Z: residual charge in atomic units
		:param int terms: terms to keep in PPT expansion'''
		super().__init__(mks_length)
		self.w0 = w0
		self.Uion = Uion
		self.Z = Z
		self.terms = terms
		self.cutoff_field = 1e-3
	def ExtractEikonalForm(self,E,dt,w00=0.0,bandwidth=1.0):
		"""Extract amplitude, phase, and center frequency from a carrier resolved field E.
		The assumed form is E = Re(amp*exp(i*phase)).

		:param numpy.array E: Electric field, any shape, axis 0 is time.
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
		'''Precompute coefficients for ionization formulas'''
		return 0.0,1.0,-1.0
	def Rate(self,Es,averaging):
		'''Get the ionization rate over all space, returned value is in atomic units.

		:param np.array Es: Time domain electric field in atomic units, any shape.  Some subclasses demand that time be the first axis.
		:param bool averaging: Interpret Es a complex envelope and get cycle averaged rate.  In this case |Es| is expected to give the crest.'''
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
	"""Klaiber dressed coulomb corrected SFA (Eq. 97)"""
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
	def Coeff(self,averaging):
		'''Derive from ADK so that in the tunneling limit,
		we can simply use Popov's note to undo application of the Stirling formula'''
		nstar = self.Z/np.sqrt(2*self.Uion)
		NPPT = (2**(2*nstar-1)/spsf.gamma(nstar+1))**2
		NADK = (1/(8*np.pi*nstar))*(4*np.exp(1)/nstar)**(2*nstar)
		C1,C2,C3 = super().Coeff(averaging)
		C1 *= NPPT/NADK
		return C1,C2,C3

class PPT(PPT_Tunneling):
	'''Full PPT theory accounting for multiphoton and tunneling limits.
	If instantaneous rate is requested, the cycle averaging factor is divided out.'''
	def Rate(self,Es,averaging):
		F0 = np.sqrt(2*self.Uion)**3
		nstar = self.Z/np.sqrt(2*self.Uion)
		lstar = nstar - 1
		C2 = 2**(2*nstar) / (nstar*spsf.gamma(nstar+lstar+1)*spsf.gamma(nstar-lstar))
		E = np.abs(Es) + self.cutoff_field
		w = self.w0
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
		ans *= (2*F0/E)**(2*nstar) # coulomb correction
		ans *= self.Uion*C2*np.sqrt(6/np.pi) * np.exp(-2.0*F0*g/(3*E))
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
