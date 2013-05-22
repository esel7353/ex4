import scipy.constants as sp

v = [0.236, 0.312, 0.391, 0.471, 0.551]
v = [x * 100 for x in v]

for F in range(0, 5):
  print(sp.h * sp.c * v[F] / (F+3) / sp.eV / sp.micro )
