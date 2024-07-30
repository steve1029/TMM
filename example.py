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
example.set_incidentangle(angle=np.pi/6, unit='radian')
print('modified incident angle : ', example.incangle)
# print '%.2e' %example.wavelength

example.set_mediumindex(1, 2.2, 1)
#example.set_mediumtype('magnetic',[1,2,3,2,1])
example.set_mediumtype('nonmagnetic')
example.set_mediumthick(500*nm)

size = (10,8)

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