import numpy as np
import scipy.constants as codata
import matplotlib.pyplot as plt
par_NM=np.array([ 2.8089956 , -4.6839549 ,  1.62960548])
def G2K(G,par=par_NM):
    lamda = 2.4  # 2.4cm 24mm
    G = G  # 6.332mm  0.6332cm
    C = G / lamda
    # B0=2.983*np.exp(   -C*(5.068-1.52*C)   )  ##B0=3.44*np.log(5.08*G/lamda-1.54*G/lamda*G/lamda)
    B0 = par[0] * np.exp(par[1]*C + par[2]*C*C)
    K = 0.9336 * B0 * lamda  # lamda convert to cm
    return K

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

def G2eV(G):
    K=G2K(G)
    e=K2eV(K)
    return e
x=np.linspace(6,10,40)
y=[]
for xi in x:
    xi=xi/10
    ei=G2eV(xi)
    y.append(ei)

y=np.array(y)
out=np.vstack([x,y])
np.savetxt('plot/gap2ev.txt',out.T)
plt.plot(x,y)
plt.show()
