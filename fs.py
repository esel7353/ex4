import scipy.constants as sc
from math import *

def E_n(n):
	return - sc.m_e * sc.e**4/(8*sc.epsilon_0**2 * sc.h**2 * n**2)

def E_nj(n, j):
	return E_n(n) * (1- sc.alpha**2 / n * (1/(j+0.5) - 3/(4*n)))

def g_j(l, s, j):
	return 1 + (j*(j+1) + s*(s+1) - l*(l+1)) / (2 * j *(j+1))

def E_anomal(n, l, s, j):
	return (E_nj(n, j) + deltaE_anomal(l, s, j), E_nj(n, j) - deltaE_anomal(l, s, j))

def deltaE_anomal(l, s, j):
	return	 g_j(l, s, j) * bohrMag

def g_dE_Ea(n, l, s, j):
	return [g_j(l, s, j), deltaE_anomal(l, s, j) / sc.eV / sc.micro, E_nj(n, j) / sc.eV]

def square(z):
	return z**2


bohrMag = sc.physical_constants["Bohr magneton"][0]
eV = sc.eV
ueV = eV * sc.micro
meV = eV * sc.milli
neV = eV * sc.nano

keV = eV * sc.kilo
MeV = eV * sc.mega
GeV = eV * sc.giga

m_e = sc.m_e
alpha = sc.alpha
e = sc.e
h = sc.h
hbar = sc.hbar

