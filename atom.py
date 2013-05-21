import <cipy.constants as sp
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
  n = 1
  l = 0
  s = 0
  j = 0.5
  F = 0
  I = 0.5

  #tools
  defaultUnitEnergy = SI

  #generally used values
  reducedMass = sp.e * sp.m_e / (sp.e + sp.m_e)
  Ry = reducedMass*sp.e**4/(8*sp.epsilon_0**2*sp.h**2)

  #List for the Orbitals
  O=["S","P","D","F","G"]

  def SetState(self,n=-1,s=-1,X="S",j=-1 ):
    if n<=0 : self.n=n;
    if s<=0 : self.s=s;
    if j<=0 : self.j=j;
    l=O.index(X)
  
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

  def EnergyHyperFine(self, n, F, j ,I):
    #own idea
    return [self.EnergyFineStructure(n, j), self.EnergieDeltaHyperFine(F, j, I)]

