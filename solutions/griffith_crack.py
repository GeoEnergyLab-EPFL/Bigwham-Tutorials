#
# Solution for the elastic fields around a straight crack [-a,a] under uniaxial (tensile) uniform loading
#  Griffith crack solution 
# see e.g.  
#Tada, H., Paris P.C., and Irwin, G.R. The Stress Analysis of Cracks Handbook. 3rd edition, 2000.
#
#Pages 18-21
#
# and e.g.  Sneddon, I. N. The distribution of stress in the neighbourhood of a crack in an elastic solid. Proc. Roy. Soc. series A., 187(1009):229–260, 1946.
# section 2 (The 2D model)

import numpy as np


    
# solution for displacement and stresses  around a mode I Griffith crack ....
def Z_griffith(z,a=1,sig=1):
    return sig*(1/np.sqrt(1-(a/z)**2)-1)

def Zprime_griffith(z,a=1,sig=1):
    return -sig*(a**2)/(((1-(a**2/z**2))**(1.5))*(z**3))

def ZIb_griffith(z,a=1,sig=1):
    return sig*z*(-1+np.sqrt(1-a**2/z**2))

# displacements
def u1_griffith(z,a=1,sig=1,G=1,nu=0.25):
    return ((1-2*nu)*np.real(ZIb_griffith(z,a,sig)) -(np.imag(z))*np.imag(Z_griffith(z,a,sig)) )/(2*G)

def u2_griffith(z,a=1,sig=1,G=1,nu=0.25):
    return (2*(1-nu)*np.imag(ZIb_griffith(z,a,sig)) -(np.imag(z))*np.real(Z_griffith(z,a,sig)) )/(2*G)

def displacement_griffith(z,a=1,sig=1,G=1,nu=0.25):
    return u1_griffith(z,a,sig,G,nu),u2_griffith(z,a,sig,G,nu)

# displacement discontinuity - cod

def width_griffith(x,a=1,sig=1.,G=1,nu=0.25):
    return (sig*2*(1-nu)/G)*np.sqrt(a*a-x*x)

     
# stresses ...

def sig11_griffith(z,a=1,sig=1,nu=0.25):
    return np.real(Z_griffith(z,a,sig))-(np.imag(z))*np.imag(Zprime_griffith(z,a,sig))

def sig22_griffith(z,a=1,sig=1,nu=0.25):
    return np.real(Z_griffith(z,a,sig))+(np.imag(z))*np.imag(Zprime_griffith(z,a,sig))

def sig12_griffith(z,a=1,sig=1,nu=0.25):
    return -(np.imag(z))*np.real(Zprime_griffith(z,a,sig))

def stress_griffith(z,a=1,sig=1,nu=0.25):
    return sig11_griffith(z,a,sig,nu),sig22_griffith(z,a,sig,nu),sig12_griffith(z,a,sig,nu)

