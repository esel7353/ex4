import scipy.constants as sp
import numpy as np

SI = 1
eV = sp.eV
keV = eV * sp.kilo
MeV = eV * sp.mega
GeV = eV * sp.giga
meV = eV * sp.milli
ueV = eV * sp.micro

class H_atom:
  #state of the atom
  j = 0.5
  n = 1
  F = 0
  I = 0.5

  #tools
  defaultUnitEnergy = SI

  #generally used values
  reducedMass = sp.e * sp.m_e / (sp.e + sp.m_e)
  Ry = reducedMass*sp.e**4/(8*sp.epsilon_0**2*sp.h**2)

  #syntax n 2s+1 S j F
  def parseState(self, state):
    tokens = state.split(" ")
    for i in range(0, len(tokens)):
      if tokens[i] in ["s", "p", "d", "f", "g", "h"]:
        lPos = i
        break
        
    for i in range(0, len(tokens)):
      if i-lPos = -2: self.n = tokens[i]; 
      if i-lPos = -1: self.n = (tokens[i]-1)/2; 
      if i-lPos = +1: self.j = tokens[i]; 
      if i-lPos = +2: self.F = tokens[i]; 
      if i-lPos = 0:
        if tokens[i] == "s":  self.l = 0;
        if tokens[i] == "p":  self.l = 1;
        if tokens[i] == "d":  self.l = 2;
        if tokens[i] == "f":  self.l = 3;
        if tokens[i] == "g":  self.l = 4;
        if tokens[i] == "h":  self.l = 5;
         


  def  EnergyBohr(self, n="this"):
    if isinstance(n, (list, tuple, range)): n = np.array(n);
    if n=="this": n=self.n; 
    
    #source Demtroeder 3, edit. 4, eq. 3.106 
    return -self.Ry/n**2 / self.defaultUnitEnergy 

  def EnergyDeltaFineStructure(self, n="this", j="this"):
    if n == "this": n=self.n;
    if j == "this": j=self.j;

    return self.EnergyBohr(n) * sp.alpha**2 / n * (1/(j+0.5) - 3/4.0/n)

  def EnergyFineStructure(self, n="this", j="this"):
    return self.EnergyBohr(n) + self.EnergyDeltaFineStructure(n, j)

  def EnergyDeltaHyperFine(self, F, j, I):
    #source Mayer-Kuckuk, edit. 4, eq. 5.51    
    return 0.5 * (F*(F+1) - j*(j+1) - I*(I+1))
