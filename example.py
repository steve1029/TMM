import TMM
import numpy as np
import matplotlib.pyplot as plt

example = TMM.TMM()

nm = 1.e-9
um = 1.e-6
# freq = np.arange(400,800) * 1e9

wl = np.arange(550,750,0.2) * nm

example.set_wavelength(wl)

print('default incident angle : ', example.incangle)

brewsterangle = np.arctan(2)
example.set_incidentangle(angle=0., unit='radian')
print('modified incident angle : ', example.incangle)
# print '%.2e' %example.wavelength

print(example.set_mediumindex(1,2,1,2,1))
# print example.mediumtype('magnetic',[1,2,3,2,1])
print(example.set_mediumtype('nonmagnetic'))
print(example.set_mediumthick(102.4*nm, 153.6*nm, 102.4*nm))

size = (16,9)

matrixs = example.cal_spol_matrix()
ms = example.matrixs
Reflecs = example.Reflectance() 
Transms = example.Transmittance()
example.graph('wavelength',figsize=size)
example.graph('frequency',figuresize=size)
matrixp = example.cal_ppol_matrix()
mp = example.matrixp
Reflecp = example.Reflectance() 
Transmp = example.Transmittance()
example.graph('wavelength',figuresize=size)
example.graph('frequency',figsize=size)
# print np.array(mp)
# print Reflecs - Reflecp
# print Transms - Transmp
# print example.__doc__
# print example.cal_ppol_matrix.__doc__