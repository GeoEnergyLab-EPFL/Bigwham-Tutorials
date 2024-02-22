
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
    
import scipy 
from scipy.special import ellipeinc, ellipkinc
# finite size region loaded from 0 to R_w
def width_finite_region(r_corr,R=1,R_w=0.5,sig=1.,G=1.,nu=0.25):
    # r_corr must be a python array 
    fact = (4*sig*(1-nu)/(np.pi*G))
    # 
    rho_smaller_rw = r_corr[np.argwhere(r_corr<R_w).flatten()]/R
    rho_larger_rw = r_corr[np.argwhere(r_corr>=R_w).flatten()]/R
    m=R_w/R
    mu_1 =  np.sqrt((1-m*m)/(1-rho_smaller_rw**2))
    w_smaller_rw = np.sqrt(1-rho_smaller_rw*rho_smaller_rw)-mu_1+m*ellipeinc(np.arcsin(mu_1),(rho_smaller_rw/m)**2)
    mu_2 =  np.sqrt((1-rho_larger_rw**2)/(1-m*m))
    w_larger_rw = np.sqrt(1-rho_larger_rw**2) - mu_2+rho_larger_rw*ellipeinc(np.arcsin(mu_2),(m/rho_larger_rw)**2) + \
                ((m*m-rho_larger_rw**2)/rho_larger_rw)*ellipkinc(np.arcsin(mu_2),(m/rho_larger_rw)**2)
    w_all = r_corr*0
    w_all[np.argwhere(r_corr<R_w).flatten()]=fact*w_smaller_rw
    w_all[np.argwhere(r_corr>=R_w).flatten()]=fact*w_larger_rw
    return w_all
