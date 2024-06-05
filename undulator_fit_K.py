import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as codata
from scipy import optimize

#from orangecontrib.xoppy.util.xoppy_undulators import xoppy_calc_undulator_spectrum

# energy, flux, spectral_power, cumulated_power = xoppy_calc_undulator_spectrum(
#     ELECTRONENERGY=3.5,
#     ELECTRONENERGYSPREAD=0.0007,
#     ELECTRONCURRENT=0.3,
#     ELECTRONBEAMSIZEH=0.0001642,
#     ELECTRONBEAMSIZEV=1.122e-05,
#     ELECTRONBEAMDIVERGENCEH=3.347e-05,
#     ELECTRONBEAMDIVERGENCEV=4.02e-06,
#     PERIODID=0.024,
#     NPERIODS=65,
#     KV=1.955,
#     KH=0.0,
#     KPHASE=0.0,
#     DISTANCE=21.0,
#     GAPH=0.0028,
#     GAPV=0.0018,
#     GAPH_CENTER=0.0,
#     GAPV_CENTER=0.0,
#     PHOTONENERGYMIN=13000.0,
#     PHOTONENERGYMAX=16000.0,
#     PHOTONENERGYPOINTS=200,
#     METHOD=2,
#     USEEMITTANCES=1)
# K=1.955
# lamda=2.4#cm
# E=3.5#GeV
# base=0.95*(E*E)/( lamda*(1+K*K/2) )
# print(base)

# dat1=np.loadtxt('plot/1.955-e-30.txt')
# x1=dat1[:,0]
# y1=dat1[:,1]
# y1=y1/max(y1)
# dat2=np.loadtxt('plot/1.955-e-20.txt')
# x2=dat2[:,0]
# y2=dat2[:,1]
# y2=y2/max(y2)
#
# # dat3=np.loadtxt('plot/1.955-e-5.txt')
# # x3=dat3[:,0]
# # y3=dat3[:,1]
# # y3=y3/max(y3)
#
# dat4=np.loadtxt('plot/1.955-0.001.txt')
# x4=dat4[:,0]
# y4=dat4[:,1]
# y4=y4/max(y4)
#
#
# plt.plot(x1,y1)
# plt.plot(x2,y2)
# #plt.plot(x3,y3)
# plt.plot(x4,y4)
# plt.legend(['e-30','e-20','0.001'])
# plt.xlim([8100,8400])
# plt.show()

def K2eV(K,shift=-24,order=5):
    #order
    #shift from experiment slit width  -24eV
    ###
    EE=3.5#bl['ElectronEnergy']
    KV = K
    KH = 0.0
    PERIODID=0.024#unit m
    codata_mee = codata.m_e * codata.c ** 2 / codata.e  # electron mass in eV
    gamma = EE * 1e9 / codata_mee
    m2ev = codata.c * codata.h / codata.e  # lambda(m)  = m2eV / energy(eV)
    resonance_wavelength = (1 + (KV ** 2 + KH ** 2) / 2.0) / 2 / gamma ** 2 * PERIODID
    resonance_energy = m2ev / resonance_wavelength
    # print('resonance_energy',resonance_energy)
    res=resonance_energy*order+shift
    return res

def obj(K):
    res=abs(K2eV(K)- 8253.955046670162)
    return res
aa=optimize.fmin_cg(obj,[1.9])