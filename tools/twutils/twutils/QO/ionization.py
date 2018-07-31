"""Functions pertaining to photoionization
Contents:
ionization potentials dictionary
class atomicUnits
def Shifted_Uion_H
def Keldysh_Parameter
def ADK_Rate
def ADK_Rate_Avg
def Keldysh_Rate
def Hydrogen_Rate
def Hydrogen_Rate_Shifted
def PPT_Rate
def YI_Rate
def Klaiber_Rate"""

#[1] A.M. Perelemov, V.S. Popov, M.V. Terent'ev, Sov. Phys. JETP 23 (5), 924 (1966)
#[2] M.V. Ammosov, N.B. Delone, V.P. Krainov, Sov. Phys. JETP 64 (6), 1191 (1986)
#[3] L.V. Keldysh, Sov. Phys. JETP 20 (5), 1307 (1965)"""

import numpy as np
from scipy import special as s
from scipy import constants as C


# data from photocopy of three pages (10-205--10-207) of unknown book
# refs listed in photocopy:
# Moore, C.E., Natl. Stand. Ref. Data Ser.--Natl. Bur. Stand. (U.S.) No. 34, 1970
# Martin, W.C. et al., Natl. Stand. Ref. Data Ser.---Natl. Bur. Stand. (U.S.) No. 60, 1978
# Sugar, J. and Corliss, C., J. Phys. Chem. Ref. Data Vol. 14, Suppl. 2, 1985
# Refs to papers in J. Phys. Chem. Ref. Data (1973-1991) may be found in cumulative index to that journal
# Cohen, E.R. and Taylor, B.N., J. Phys. Chem. Ref. Data 17, 1795 (1988)

ip = {'H' : 13.59844,
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

class atomicUnits:
	"""Contains atomic units expressed in mks units"""
	length = C.hbar/(C.m_e*C.c*C.alpha)
	time = C.hbar/(C.m_e*C.c**2*C.alpha**2)
	energy = C.m_e*C.c**2*C.alpha**2
	field = energy/length/C.e

def Shifted_Uion_H(E):
	"""Returns the Stark shifted ionization energy of Hydrogen"""
	return 0.5 + (9.0/4.0)*E**2 + 55.55*E**4 + 4908.0*E**6 + 794200*E**8
	
def Keldysh_Parameter(omega,Uion,E):
	"""Returns the Keldysh parameter
	inputs (a.u.): omega , Uion , E"""
	return omega*np.sqrt(2.0*Uion)/E

def Cycle_Averaging_Factor(Uion,E):
	"""Returns the PPT cycle averaging factor
	Multiply static rate by this factor to get average rate
	inputs (a.u.):
	Uion = ionization potential
	E = electric field"""
	return np.sqrt(3.0/np.pi) * np.sqrt(E) / (2*Uion)**0.75
	
def ADK_Rate(Uion,Z,E):
	"""Returns ADK rate in a.u.
	inputs (a.u.):
	Uion = ionization potential
	Z = residual charge
	E = electric field
	This gives the instantaneous rate
	This can be obtained from Ref. 1, Eq. 4, l = m = 0, with C_kl coming from Ref. 2, Eq. 11."""
	nstar = Z/np.sqrt(2*Uion)
	D = ((4.0*np.exp(1.0)*Z**3)/(E*nstar**4))**nstar
	return (E*D*D/(8.0*np.pi*Z))*np.exp(-2*Z**3/(3.0*nstar**3*E))
	
def ADK_Rate_Avg(Uion,Z,E):
	"""Similar to ADK_Rate, except averages optical cycles"""
	return Cycle_Averaging_Factor(Uion,E)*ADK_Rate(Uion,Z,E)
	
def Keldysh_Rate(Uion,Z,E):
	"""Similar to ADK_Rate, but uses Keldysh formula
	The Keldysh rate has problems, and should only be used for comparison"""
	ans = np.sqrt(6.0*np.pi)/4.0
	ans *= Uion * np.sqrt(E/(Uion**1.5))
	ans *= np.exp(-(4.0/3.0)*np.sqrt(2.0)*(Uion**1.5)/E)
	return ans

def Hydrogen_Rate(E):
	"""Similar to ADK_Rate, but assumes hydrogen (Landau formula)"""
	return (4.0/E)*np.exp(-(2.0/3.0)/E)
	
def Hydrogen_Rate_Shifted(E):
	"""Ionization of Hydrogen, accounting for Stark shift"""
	Uion = Shifted_Uion_H(E)
	return (4.0/E) * ((2*Uion)**2.5) * np.exp(-(2.0/3.0)*((2*Uion)**1.5)/E)

def LL_Rate(Uion,E):
	"""Ionization of Hydrogen, Uion manually"""
	return (4.0/E) * ((2*Uion)**2.5) * np.exp(-(2.0/3.0)*((2*Uion)**1.5)/E)

def w0(x):
	"""Exponential integral from PPT theory"""
	return np.real(0.5 * np.exp(-x**2) * np.sqrt(np.pi) * s.erf(x*1j) / 1j)

def PPT_Rate(Uion,Z,E,w,terms):
	"""PPT cycle averaged ionization rate
	inputs (a.u.):
	Uion = ionization potential
	Z = residual charge
	E = electric field
	w = radiation frequency
	terms = terms to keep in PPT expansion"""
	F0 = np.sqrt(2*Uion)**3
	nstar = Z/np.sqrt(2*Uion)
	lstar = nstar - 1
	C2 = 2**(2*nstar) / (nstar*s.gamma(nstar+lstar+1)*s.gamma(nstar-lstar))
	gam = np.sqrt(2.0*Uion)*w/E
	alpha = 2 * (np.arcsinh(gam)-gam/np.sqrt(1+gam**2))
	beta = 2*gam/np.sqrt(1+gam**2)
	g = (3/(2*gam))*((1+1/(2*gam**2))*np.arcsinh(gam)-np.sqrt(1+gam**2)/(2*gam))
	nu = (Uion/w) * (1 + 1/(2*gam**2))
	A0 = 0
	for n in range(np.int(np.ceil(nu)),np.int(np.ceil(nu)+terms)):
		A0 += np.exp(-alpha*(n-nu))*w0(np.sqrt(beta*(n-nu)))
	A0 *= (4/np.sqrt(3*np.pi)) * (gam**2/(1+gam**2))
	ans = A0*(E*np.sqrt(1+gam**2)/(2*F0))**1.5
	ans *= (2*F0/E)**(2*nstar) # coulomb correction
	ans *= Uion*C2*np.sqrt(6/np.pi) * np.exp(-2.0*F0*g/(3*E))
	return ans
	
def YI_Rate(Uion,Z,E,w,phase,terms):
	"""Yudin-Ivanov instantaneous phase dependent ionization rate
	inputs (a.u.):
	Uion = ionization potential
	Z = residual charge
	E = electric field envelope (peak)
	w = radiation frequency
	phase = phase, defined such that 0 is at peak of field
	terms = terms to keep in PPT expansion"""
	nstar = Z/np.sqrt(2*Uion)
	lstar = nstar - 1
	Anl = 2**(2*nstar) / (nstar*s.gamma(nstar+lstar+1)*s.gamma(nstar-lstar))
	theta = (phase - 0.5*np.pi)%np.pi - 0.5*np.pi
	gam = np.sqrt(2.0*Uion)*w/E
	a = 1+gam*gam-np.sin(theta)**2
	b = np.sqrt(a*a+4*gam*gam*np.sin(theta)**2)
	c = np.sqrt((np.sqrt((b+a)/2)+gam)**2 + (np.sqrt((b-a)/2)+np.sin(np.abs(theta)))**2)
	Phi = (gam**2 + np.sin(theta)**2 + 0.5)*np.log(c)
	Phi -= 3*(np.sqrt(b-a)/(2*np.sqrt(2)))*np.sin(np.abs(theta))
	Phi -= (np.sqrt(b+a)/(2*np.sqrt(2)))*gam
	kappa = np.log(gam+np.sqrt(gam**2+1)) - gam/np.sqrt(1+gam**2)
	alpha = 2 * (np.arcsinh(gam)-gam/np.sqrt(1+gam**2))
	beta = 2*gam/np.sqrt(1+gam**2)
	nu = (Uion/w) * (1 + 1/(2*gam**2))
	A0 = 0
	for n in range(np.int(np.ceil(nu)),np.int(np.ceil(nu)+terms)):
		A0 += np.exp(-alpha*(n-nu))*w0(np.sqrt(beta*(n-nu)))
	A0 *= (4/np.sqrt(3*np.pi)) * (gam**2/(1+gam**2))
	pre = Anl * np.sqrt(3*kappa/gam**3)*(1+gam**2)**0.75 * A0 * Uion
	pre *= (2*(2*Uion)**1.5 / E)**(2*nstar-1)
	return pre * np.exp(-E**2 * Phi / w**3)

def Klaiber_Rate_Avg(Uion,Eau):
	"""Klaiber dressed coulomb corrected SFA (Eq. 97)
	inputs (a.u.)
	Uion = ionization potential
	Eau = static electric field"""
	
	Ua2 = Uion*C.alpha**2
	Ea = (2*Uion)**1.5 # definition of Ea applies to Eq. 97 but not Fig. 4 caption
	f1 = 2**(3-4*Ua2)*np.sqrt(3/np.pi)*(1-7*Ua2/72)*np.exp(4*Ua2)*(2*Uion)**(1.75-3*Ua2)
	f2 = s.gamma(3-2*Ua2)*Eau**(0.5-2*Ua2)
	f3 = np.exp(-0.667*Ea*(1-Ua2/12)/Eau)
	return (f1/f2)*f3

def Klaiber_Rate(Uion,Eau):
	"""Klaiber dressed coulomb corrected SFA (Eq. 97)
	Includes a factor from PPT to unwind cycle averaging
	inputs (a.u.)
	Uion = ionization potential
	Eau = static electric field"""

	return Klaiber_Rate_Avg(Uion,Eau)/Cycle_Averaging_Factor(Uion,Eau)	
