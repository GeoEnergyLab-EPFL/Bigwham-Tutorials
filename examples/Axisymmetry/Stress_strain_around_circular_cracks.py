
#%% imports etc.
import sys
import os
home = os.environ["HOME"]


from bigwham4py import BEMatrix
import numpy as np
from scipy.sparse.linalg import gmres

sys.path.append(os.path.join(os.getcwd(), '../..'))
from solutions.penny_shaped_crack_uniform import *
import matplotlib.pyplot as plt


# %%
def compliance_tensor(E,nu):
    """Returns the compliance tensor for isotropic linear elasticity in 3D.

    Args:
        E (float): Young's modulus
        nu (float): Poisson's ratio

    Returns:
        np.array: 6x6 compliance tensor
    """
    S = np.zeros((6, 6))
    factor = 1 / E
    S[0, 0] = S[1, 1] = S[2, 2] = factor
    S[0, 1] = S[0, 2] = S[1, 0] = S[1, 2] = S[2, 0] = S[2, 1] = -nu * factor
    S[3, 3] = S[4, 4] = S[5, 5] = 2 * (1 + nu) * factor
    return S


#%% material properties
G = 3.*1.e9  # shear modulus
nu = 0.25
E = (2 * G) * (1 + nu)

# compute S
S = compliance_tensor(E, nu)
# crack size
a = 0.04 

# line of observation points

rr =np.linspace(0.00001,3*a,1000)
z_obs = 0.015

#2d grid of points
x = np.linspace(-3*a,3*a,100)
y = np.linspace(-3*a,3*a,20)
xx,yy = np.meshgrid(x,y)
xy=np.vstack([xx.ravel(),yy.ravel()]).T

#%%%% SHear loading  GRID 
ct_shear = 1e6

sxx,syy,szz,sxy,sxz,syz  =stresses_shear_penny_shaped(xy[:,0],xy[:,1],z_obs,a=a,T_x=ct_shear,T_y=0.0,G=G,nu=nu)

#plt.plot(rr,sxx,label='sxx')
#plt.title(' Load over a finite patch )')

allStress= np.array([sxx,syy,szz,sxy,sxz,syz])

allStrain=allStress
for i in range(len(xy[:,0])):
    allStrain[:,i]=np.dot(S,allStress[:,i])


# plt.plot(rr,allStrain[0,:]*1e6,label='exx')
# plt.title('Eps_xx (micro-strain) - circular crack under shear ct loading')

Z=allStrain[0,:].reshape(xx.shape)*1e6


# This custom formatter removes trailing zeros, e.g. "1.0" becomes "1", and
# then adds a percent sign.
def fmt(x):
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return rf"{s} " if plt.rcParams["text.usetex"] else f"{s} "

# Basic contour plot
fig, ax = plt.subplots()
CS = ax.contour(xx,yy,Z)

ax.clabel(CS, CS.levels, fmt=fmt, fontsize=10)

#%% plot for different y 

plt.plot(xy[:,0],allStrain[0,:]*1e6,'b',label='exx')
plt.title('Circular crack under constant shear loading')
plt.xlabel('x (m)')
plt.ylabel('$\epsilon_{xx}$ (micro-strain)')

#%% plot for different x

plt.plot(xy[:,1],allStrain[1,:]*1e6,'b',label='eyy')
plt.title('$\epsilon_{xx}$ (micro-strain) - circular crack under tensile ct loading')


# %% just along one radius

sxx,syy,szz,sxy,sxz,syz  =stresses_shear_penny_shaped(rr*1,rr*0.,z_obs,a=a,T_x=ct_shear,T_y=0.0,G=G,nu=nu)

allStress= np.array([sxx,syy,szz,sxy,sxz,syz])

allStrain=allStress
for i in range(len(rr)):
    allStrain[:,i]=np.dot(S,allStress[:,i])


plt.plot(rr,allStrain[0,:]*1e6,label='exx')
plt.title('$\epsilon_{xx}$ (micro-strain) - circular crack under tensile ct loading')


#%%%% TENSILE loading 
ct_tens = ct_shear

sxx,syy,szz,sxy,sxz,syz  =stresses_tensile_penny_shaped(rr*1,rr*0.,z_obs,a=a,sig=ct_tens,G=G,nu=nu)

# # %%
# plt.plot(rr,sxx,label='sxx')
# plt.title(' Load over a finite patch )')

allStress= np.array([sxx,syy,szz,sxy,sxz,syz])

allStrain=allStress
for i in range(len(rr)):
    allStrain[:,i]=np.dot(S,allStress[:,i])

plt.plot(rr,allStrain[0,:]*1e6,label='exx')
plt.title('Eps_xx (micro-strain) - circular crack under tensile ct loading')


# %%
# estimate of compressuve disop
dl= (ct_shear/E)*0.25
# %%
