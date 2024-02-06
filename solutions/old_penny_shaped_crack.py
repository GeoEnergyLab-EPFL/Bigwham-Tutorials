import numpy as np



def L1(r,z,a):
    return 0.5*(-np.sqrt ((a - r) ** 2 + z ** 2) + np.sqrt ((a + r) ** 2 + z ** 2))
def L2(r,z,a):
    return 0.5*(np.sqrt ((a - r) ** 2 + z ** 2) +  np.sqrt ((a + r) ** 2 + z ** 2))


### STRESSES and displacement Solutions from Kachanov's book !!! 
def stresses_tensile_penny_shaped(r,z,a=1.0,sig=1.0,nu=0.25):
    """Stresses (in cylindrical coordinate) due to a penny-shaped crack under uniform tensile loading
    Solution taken from M. Kachanov, B. Shafiro, and I. Tsukrov. Handbook of Elasticity Solutions. Kluwer, 2003. - pg 92-93
    Args:
        r (flota): radial coordinates 
        z (_type_): z coordinates
        a (float, optional): Crack radius . Defaults to 1.
        sig (float, optional): Applied internal Pressure. Defaults to 1.
        nu (float, optional): Material poisson's ratio. Defaults to 0.25.

    Returns:
       srr,stt,szz,srt,srz,stz (tuple of np.array): component of the stress tensor in  cylindrical coordinates
    """
    # sig_1 = s_xx+s_yy   
    # sig_2 = s_xx-xyy +2 i s_xy
    # sig_z = s_zz
    # tau_z = s_xz +i syz
    # we return for \theta = 0 (phi in Kachanov notation)
    fact = 2.*sig/np.pi
    l1 = L1(r,z,a)
    l2 = L2(r,z,a)
    s1= (1+2*nu)*(a*np.sqrt(l2**2-a**2)/(l2*l2-l1*l1) - np.arcsin(a/l2))+a*(z*z)*(l1**4+a*a*(2*a*a+2*z*z-3*r*r))/(((l2*l2-l1*l1)**3)*np.sqrt(l2*l2-a*a)) 
    s2 = a*l1*l1*np.sqrt(l2*l2-a*a)/(l2*l2*(l2*l2-l1*l1))*(1-2*nu+\
        z*(a*a*(6*l2*l2-2*l1*l1+r*r)-5*(l2**4))/(((l2*l2-l1*l1)**2)*(l2*l2-a*a) )    )  # at phi-=0
    sz = a*np.sqrt(l2**2-a**2)/(l2*l2-l1*l1)-np.arcsin(a/l2)-\
        a*(z*z)*(l1**4+a*a*(2*a*a+2*z*z-3*r*r))/(((l2*l2-l1*l1)**3)*np.sqrt(l2*l2-a*a))
    tz=-z*l1*np.sqrt(l2*l2-a*a)*(a*a*(4*l2*l2-5*r*r)+l1**4)/(l2*((l2*l2-l1*l1)**3) )
    # for theta = 0, we have
    # s_rr = 0.5*(s1+s2), s_tt = 0.5(s1-s2), srz_=tz
    return fact*(s1+s2)/2,fact*(s1-s2)/2,fact*sz,0.*tz,fact*tz,0.*tz


def displacement_tensile_penny_shaped(r,z,a=1.0,sig=1.0,nu=0.25,G=1.0):
    """Displacement (in cylindrical coordinates) due to a penny-shaped crack under uniform tensile loading
    Solution taken from M. Kachanov, B. Shafiro, and I. Tsukrov. Handbook of Elasticity Solutions. Kluwer, 2003. - pg 92-93

    Args:
        r (flota): radial coordinates 
        z (_type_): z coordinates
        a (float, optional): Crack radius . Defaults to 1.
        sig (float, optional): Applied internal Pressure. Defaults to 1.
        nu (float, optional): Material poisson's ratio. Defaults to 0.25.
        G (float, optional): Material  shear modulus. Defaults to 1.

    Returns:
        u_r, u_t, u_z (tuple of np.array of floats): component of displacements in  cylindrical coordinates
    """
    # kachanov u = ux + iuy = u_rr for theta =0
    
    fact = sig/(2.*np.pi*G)
    l1 = L1(r,z,a)
    l2 = L2(r,z,a)
    
    u_r=fact*r*((1-2*nu)*(a*np.sqrt(l2*l2-a*a)/(l2*l2)-np.arcsin(a/l2))+2*a*a*np.abs(z)*np.sqrt(a*a-l1*l1)/(l2*l2*np.sqrt(l2*l2-l1*l1)))
    u_t=u_r*0.
    
    u_z = 2*fact*(2*(1-nu)*((z/np.abs(z))*np.sqrt(a*a-l1*l1)-z*np.arcsin(a/l2))+z*(np.arcsin(a/l2)-a*np.sqrt(l2*l2-a*a)/(l2*l2-l1*l1)) )

    return u_r,u_t,u_z


# Solution taken from Nikolskiy's thesis (2016) appendix D 
# see also  A. Green. On boussinesq’s problem and penny-shaped cracks. In Mathematical Proceedings of the Cambridge Philosophical Society, volume 45, pages 251–257. Cambridge University Press, 1949.
#  Sneddon, I. N. The distribution of stress in the neighbourhood of a crack in an elastic solid. Proc. Roy. Soc. series A., 187(1009):229–260, 1946.
# there is seem to have a 'bug' in the u_r and s_rr and s_tt components ! -> Do not use
def psi(r,z,a):
    return  (((2*a**2 - (3*(-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2)/4.)*np.sqrt(-a**2 + (np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2/4.))/a +
             (-2*a**2 + r**2 - 2*z**2)*np.arcsin((-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))/(2.*r)))/4.

# first derivatives in cylindrical coordinates system 

def d_psi_dr(r,z,a):
    return  -0.5*((-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))*np.sqrt(-a**2 + (np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2/4.))/(np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2)) \
        + (r*np.arcsin((-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))/(2.*r)))/2.

def d_psi_dz(r,z,a):
    return    (a*z*(a**2 - 2*r**2 + z**2 + np.sqrt(a**2 - 2*a*r + r**2 + z**2)*np.sqrt((a + r)**2 + z**2)))/ \
        (2.*np.sqrt(2)*np.sqrt(a**2 - 2*a*r + r**2 + z**2)*np.sqrt((a + r)**2 + z**2)*np.sqrt(-a**2 + r**2 + z**2 + np.sqrt(a**2 - 2*a*r + r**2 + z**2)*np.sqrt((a + r)**2 + z**2))) + \
       ((-2*a**2 + r**2 - 2*z**2)*(-(z/np.sqrt((a - r)**2 + z**2)) + z/np.sqrt((a + r)**2 + z**2)))/(8.*r*np.sqrt(1 - (np.sqrt((a - r)**2 + z**2) - np.sqrt((a + r)**2 + z**2))**2/(4.*r**2))) -  \
       (3*(-(z/np.sqrt((a - r)**2 + z**2)) + z/np.sqrt((a + r)**2 + z**2))*(-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))*np.sqrt(-a**2 + (np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2/4.))/(8.*a) - \
       z*np.arcsin((-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))/(2.*r))

# second derivatives
def d_psi_drdr(r,z,a):
    return     -0.5*(a*(1 + (-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2/(np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2)*np.sqrt(-a**2 + (np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2/4.))/ \
        (-0.25*(-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2 + (np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2/4.) + np.arcsin((-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))/(2.*r))/2.

def d_psi_drdz(r,z,a):
    return      -((a*(-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))*np.sqrt(a**2 - (-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2/4.))/ \
        ((np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))*(-0.25*(-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2 + (np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2/4.)))

def d_psi_dzdz(r,z,a):
    return (a*np.sqrt(-a**2 + (np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2/4.))/(-0.25*(-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2 + (np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2/4.) - \
        np.arcsin((-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))/(2.*r))

# third derivatices
def d_psi_drdrdr(r,z,a):
    return  -0.5*(np.sqrt(-a**2 + (np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2/4.)* \
        ((-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))/(np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2)) + \
            (8*a**2*(a**2 - r**2 + z**2)*(-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2)))/ \
                ((np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**3*(-0.25*(-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2 + (np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2/4.)) + \
                    (a*r*(1 + (-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2/(np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2)* \
                        (1 - (2*(-a**2 + r**2 + z**2))/(-0.25*(-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2 + (np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2/4.)))/ \
                            (-0.25*(-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2 + (np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2/4.)))/ \
                                (-0.25*(-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2 + (np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2/4.)**3

def d_psi_drdrdz(r,z,a):
    return (4*a**2*z*((4*r**2*(-a**2 + (np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2/4.))/ \
        (-0.25*(-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2 + (np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2/4.)**3 - \
            (a**2 - z**2)/(-0.25*(-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2 + (np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2/4.)**2))/\
                ((np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2*np.sqrt(-a**2 + (np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2/4.))
    
def d_psi_drdzdz(r,z,a):
    return (2*a*np.sqrt(r**2 - (-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2/4.)* \
        ((-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**4/16. + a**2*(-5*r**2 + (np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2)))/ \
      ((np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))*(-0.25*(-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2 + (np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2/4.)**3)
   
def d_psi_dzdzdz(r,z,a):
    return  (np.sqrt(a**2 - (-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2/4.)*(a**2*(2*a**2 - 3*r**2 + 2*z**2) + (-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**4/16.))/ \
      (-0.25*(-np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2 + (np.sqrt((a - r)**2 + z**2) + np.sqrt((a + r)**2 + z**2))**2/4.)**3 

### derivatives in cartesian system 

### stresses in Cylindrical coordinates 
def Green_stresses_tensile_penny_shaped(r,z,a=1.0,sig=1.0,nu=0.25):
    """Stresses (in cylindrical coordinate) due to a penny-shaped crack under uniform tensile loading
        Green solution (srr and stt are BUGGY!!)
    Args:
        r (flota): radial coordinates 
        z (_type_): z coordinates
        a (float, optional): Crack radius . Defaults to 1.
        sig (float, optional): Applied internal Pressure. Defaults to 1.
        nu (float, optional): Material poisson's ratio. Defaults to 0.25.

    Returns:
       srr,stt,szz,srt,srz,stz (tuple of np.array): component of the stress tensor in  cylindrical coordinates
    """
    #
    fact = -1.*sig/np.pi
    kappa = 3-4*nu
    dpsi_rr=d_psi_drdr(r,z,a)
    dpsi_zz = d_psi_dzdz(r,z,a)
    dpsi_rrz = d_psi_drdrdz(r,z,a)
    dpsi_zzz = d_psi_dzdzdz(r,z,a)
    dpsi_rzz = d_psi_drdzdz(r,z,a)
    
    srr = fact*((kappa-1)*dpsi_rr+(kappa-3)*dpsi_zz+2*z*dpsi_rrz)
    srt =  0.*srr
    srz = fact*2*z*dpsi_rzz
    stt = fact*(-2*z*(dpsi_rrz+dpsi_zzz)-(kappa-1)*dpsi_rr-2*dpsi_zz) #fact*(  (kappa-5.)*dpsi_zz-2*z*dpsi_zzz)-srr 
    stz = 0.*srr
    szz = fact*(-2.*dpsi_zz+2*z*dpsi_zzz)

    return  srr,stt,szz,srt,srz,stz   

### displacements in Cylindrical coordinates 
def Green_displacement_tensile_penny_shaped(r,z,a=1.0,sig=1.0,nu=0.25,G=1.0):
    """Displacement (in cylindrical coordinates) due to a penny-shaped crack under uniform tensile loading
        u_r do not match's Kachanov's book solution
    Args:
        r (flota): radial coordinates 
        z (_type_): z coordinates
        a (float, optional): Crack radius . Defaults to 1.
        sig (float, optional): Applied internal Pressure. Defaults to 1.
        nu (float, optional): Material poisson's ratio. Defaults to 0.25.
        G (float, optional): Material  shear modulus. Defaults to 1.

    Returns:
        u_r, u_t, u_z (tuple of np.array of floats): component of displacements in  cylindrical coordinates
    """
    fact = -1.*sig/(2*np.pi*G)
    kappa = 3-4*nu
    dpsi_r=d_psi_dr(r,z,a)
    dpsi_rz=d_psi_drdz(r,z,a)
    dpsi_z= d_psi_dz(r,z,a)
    dpsi_zz = d_psi_dzdz(r,z,a)
    
    u_r = fact*((kappa-1)*dpsi_r+2*z*dpsi_rz)
    u_t = 0.*u_r
    u_z = fact*(-(kappa+1)*dpsi_z+2*z*dpsi_zz)
    return u_r,u_t,u_z


    

