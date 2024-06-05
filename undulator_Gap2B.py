import numpy as np
#import scipy.optimize
from scipy.optimize import minimize
from scipy.optimize import fmin_cg

# dat=np.loadtxt('1plot.txt')
# x=dat[:,0]
# y=dat[:,1]
# max=max(y)
# y1=y/max
# out=np.vstack([x,y1])
# np.savetxt('1plot_norm.txt',out.T)

lamda=2.4#2.4cm 24mm
G=0.6332#6.332mm  0.6332cm
#G=0.6
C=G/lamda
#B0=2.983*np.exp(   -C*(5.068-1.52*C)   )  ##B0=3.44*np.log(5.08*G/lamda-1.54*G/lamda*G/lamda)
B0=2.983*np.exp(   -5.068*C+1.52*C*C   )
K=0.9336*B0*lamda#lamda convert to cm
#print('B0',B0,' K ',K)

def G2K(G,par=np.array([2.983,-5.068,1.52])):
    lamda = 2.4  # 2.4cm 24mm
    G = G  # 6.332mm  0.6332cm
    C = G / lamda
    # B0=2.983*np.exp(   -C*(5.068-1.52*C)   )  ##B0=3.44*np.log(5.08*G/lamda-1.54*G/lamda*G/lamda)
    B0 = par[0] * np.exp(par[1]*C + par[2]*C*C)
    K = 0.9336 * B0 * lamda  # lamda convert to cm
    return K

def obj(par):
    #0.66 1.96352251
    #0.665 1.94768854
    #0.68  1.90272919
    K1=G2K(0.66,par)
    K2 = G2K(0.665, par)
    K3 = G2K(0.68, par)
    res=abs(K1-1.96352251)+abs(K2-1.94768854)+abs(K3-1.90272919)
    return res

def obj_print(par):
    #0.66 1.96352251
    #0.665 1.94768854
    #0.68  1.90272919
    K1=G2K(0.66,par)
    K2 = G2K(0.665, par)
    K3 = G2K(0.68, par)
    res=abs(K1-1.96352251)+abs(K2-1.94768854)+abs(K3-1.90272919)
    print(K1,'###1.96352251 ',K2,'###1.94768854',K3,'###1.90272919')
    return res
par_NM=np.array([ 2.8089956 , -4.6839549 ,  1.62960548])

x0=np.array([2.983,-5.068,1.52])
#aa=fmin_cg(obj,np.array([2.5,-5.068,1.52]))
algo='Nelder-Mead'#'Nelder-Mead'
aa=minimize(obj,x0,method=algo)
r1=obj_print(aa.x)
r2=obj_print(x0)
print(r1,'###initial ',r2)





# lamda=2.4#cm
# G=6.0#mm
# B0=0.96#T
# Br=B0*(np.log(100*np.pi*G/lamda))
