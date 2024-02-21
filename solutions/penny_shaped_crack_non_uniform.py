
# A circular / penny-shaped crack in an infinite isotropic medium 
# 
#  MODE I solutions for different loading (and so-also some Dugdale CZM solutions....)
#  (or mode II/III for nu=0, uni-directional shear solution if nu=0)
import numpy as np

# the classical uniform loading
def width_uniform(r_corr,R=1,sig=1.,G=1.,nu=0.25):
    return (4*sig*(1-nu)/(np.pi*G))*np.sqrt(R*R-r_corr*r_corr)


# Point Load 
def width_ptload(r_corr,R=1,P=1.,G=1.,nu=0.25):
    fact = (2*P*(1-nu)/(np.pi*np.pi*G))
    return fact*np.arccos(r_corr/R)*(1/r_corr)
    
