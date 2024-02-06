
# Displacements and stress around a circular crack in an infinite medium 
#
# # translated from C++ using Google BARD by B. Lecampion
#  original C++ code & mma code - C. Peruzzo 
#
# references:
#   
# [1] Nikolskiy D., Three-dimensional boundary element analysis of fractured rock,PhD thesis; 2016
# [2] Barber J.R., Solid mechanics and its applications : Elasticity.London, New York : Springer Science; 2010
# [3] Fabrikant VI., Penny-shaped crack revisited: closed-form solutions. Philosophical Magazine A. 1987; 56 (2): 191-207
# [4] Green A.E., Zerna W., Theoretical Elasticity. London: Oxford University Press; 1968
# [5] Green A.E., On Boussinesq's problem and penny-shaped cracks. Mathematical Proceedings of the Cambridge Philosophical Society. 1949; 45 (2): 251-257
# [6] Segedin C.M., Note on a penny-shaped crack under shear. Mathematical Proceedings of the Cambridge Philosophical Society 1951; 47 (2): 396-400
#
# comments (C. Peruzzo):
#
# Ref [5], [4] express the solution of a penny shaped crack under constant internal pressure. 
# The expression of the potential \[Phi] is expressed in a compact way by [3]. 
# Ref [6] solved the problem considering a constant shear and notably pointed out that the problem solved in [5] reduces to the same differential equation 
# but with different boundary conditions. As a consequence \[Phi]  can be used to express the solution derived in [6]. 
# These concepts are briefly summarized by the appendix D in [1].
# Do note that the displacement field in [4] at $5.7 and in [2] $21.5.1 contain the same typo: a factor 2\[Mu] that multiply all the displacements. The error can easily be detected by comparison with [5] or by deriving the stresses from the displacements using the stress-displacement relations. 
#
# Additional comments (B. Lecampion)
# the solution can also be found in the book M. Kachanov, B. Shafiro, and I. Tsukrov. Handbook of Elasticity Solutions. Kluwer, 2003. - pg 92-93
# (possible typos nevertheless in that book)
#
# the expressions here have been tripled check.
#
#  
import numpy as np

def pi_():
    return np.pi

def r_(x, y):
    return np.sqrt(np.square(x) + np.square(y))

def theta_(x, y):
    return np.arctan2(y, x)

def f1_(x, y, z, a):
    r = r_(x, y)
    zz = np.square(z)
    f1 = np.sqrt((a + r) * (a + r) + zz)
    return f1

def f2_(x, y, z, a):
    r = r_(x, y)
    zz = np.square(z)
    f2 = np.sqrt((a - r) * (a - r) + zz)
    return f2

def l1_(x, y, z, a):
    f1 = f1_(x, y, z, a)
    f2 = f2_(x, y, z, a)
    return 0.5 * (f1 - f2)

def l2_(x, y, z, a):
    f1 = f1_(x, y, z, a)
    f2 = f2_(x, y, z, a)
    return 0.5 * (f1 + f2)

def l1_l1_(x, y, z, a):
    l1 = l1_(x, y, z, a)
    return np.square(l1)

def l2_l2_(x, y, z, a):
    l2 = l2_(x, y, z, a)
    return np.square(l2)

def one_m_2nu(nu):
    return 1 - 2 * nu

def f3_(x, y, z, a):
    l2l2 = l2_l2_(x, y, z, a)
    return np.sqrt(l2l2 - a * a)

def f4_(x, y, z, a):
    l1l1 = l1_l1_(x, y, z, a)
    return np.sqrt(a * a - l1l1)

def f5_(x, y, z, a):
    l2l2 = l2_l2_(x, y, z, a)
    f3 = f3_(x, y, z, a)
    return a * f3 / l2l2

def ux_ux_nl_base_function(x, y, z, a, G, nu, pz):
    l2 = l2_(x, y, z, a)
    l1l1 = l1_l1_(x, y, z, a)
    l2l2 = l2_l2_(x, y, z, a)
    pi = pi_()
    c1 = pz / (2 * pi * G)
    c2 = one_m_2nu(nu)
    c3 = f5_(x, y, z, a)
    c4 = 2 * a * a * np.abs(z) * f4_(x, y, z, a)
    c5 = l2l2 * (l2l2 - l1l1)
    return c1 * (c2 * (c3 - np.arcsin(a / l2)) + c4 / c5)

def sign_(x):
    if x >= 0: return 1
    else: return -1

def f6_(x, y, z, a, nu):
    l2 = l2_(x, y, z, a)
    c1 = 4 * nu - 5
    c2 = 4 * (1 - nu)
    return c1 * z * np.arcsin(a / l2) + c2 * np.sign(z) * f4_(x, y, z, a)

def f7_(x, y, z, a):
    l2l2 = l2_l2_(x, y, z, a)
    l1l1 = l1_l1_(x, y, z, a)
    return a * z * f3_(x, y, z, a) / (l2l2 - l1l1)
#    /*
#      *
#      * DISPLACEMENTS DUE TO AN UNIFORM NORMAL LOADING pz OVER A CIRCULAR CRACK OF RADIUS a
#      * expressions valid for z >=0
#      *
#      */

def ux_nl_(x, y, z, a, G, nu, pz):
    return x * ux_ux_nl_base_function(x, y, z, a, G, nu, pz)

def uy_nl_(x, y, z, a, G, nu, pz):
    return y * ux_ux_nl_base_function(x, y, z, a, G, nu, pz)

def uz_nl_(x, y, z, a, G, nu, pz):
    l2 = l2_(x, y, z, a)
    l1l1 = l1_l1_(x, y, z, a)
    l2l2 = l2_l2_(x, y, z, a)
    pi = pi_()
    c1 = pz / (pi * G)
    c2 = 2. * (1 - nu)
    f3 = f3_(x, y, z, a)  # sqrt(l2l2 - a*a)
    f4 = f4_(x, y, z, a)  # sqrt(a*a - l1l1)
    c3 = np.sign(z) * f4
    c4 = np.arcsin(a / l2)
    c5 = (a * f3) / (l2l2 - l1l1)
    return c1 * (c2 * (c3 - z * c4) + z * (c4 - c5))

#  /*
#      *
#      * DISPLACEMENTS DUE TO AN UNIFORM SHEAR LOADING px and py OVER A CIRCULAR CRACK OF RADIUS a
#      * expressions valid for z >=0
#      *
#      */

def ux_sl_(x, y, z, a, G, nu, px, py):
    pi = pi_()
    c1 = 1. / (pi * G * (2. - nu))
    f6 = f6_(x, y, z, a, nu)
    f7 = f7_(x, y, z, a)
    theta = theta_(x, y)
    l1l1 = l1_l1_(x, y, z, a)
    l2l2 = l2_l2_(x, y, z, a)
    c2 = (px + (px * np.cos(2 * theta) + py * np.sin(2 * theta)) * l1l1 / l2l2)

    return c1 * (f6 * px + f7 * c2)

def uy_sl_(x, y, z, a, G, nu, px, py):
    pi = pi_()
    c1 = 1. / (pi * G * (2. - nu))
    f6 = f6_(x, y, z, a, nu)
    f7 = f7_(x, y, z, a)
    theta = theta_(x, y)
    l1l1 = l1_l1_(x, y, z, a)
    l2l2 = l2_l2_(x, y, z, a)
    c2 = (py + (px * np.sin(2 * theta) - py * np.cos(2 * theta)) * l1l1 / l2l2)

    return c1 * (f6 * py + f7 * c2)

def uz_sl_(x, y, z, a, G, nu, px, py):
    pi = pi_()
    c1 = (2. * x * px + 2. * y * py) / (pi * G * (2. - nu))
    c2 = (1 - 2. * nu) / 2.
    l2 = l2_(x, y, z, a)
    l2l2 = l2 * l2
    f3 = f3_(x, y, z, a)
    c3 = (np.arcsin(a / l2) - a * f3 / l2l2)
    f4 = f4_(x, y, z, a)
    l1l1 = l1_l1_(x, y, z, a)
    c4 = z * a * a * f4 / (l2l2 * (l2l2 - l1l1))

    return c1 * (c2 * c3 + c4)


    # /*
    #  *
    #  * STRESSES DUE TO AN UNIFORM NORMAL LOADING pz OVER A CIRCULAR CRACK OF RADIUS a
    #  * expressions valid for z >=0
    #  *
    #  */


def f8_(x, y, z, a):
    f3 = f3_(x, y, z, a)  # sqrt(l2l2 - a*a)
    l2l2 = l2_l2_(x, y, z, a)
    l2 = l2_(x, y, z, a)
    l1l1 = l1_l1_(x, y, z, a)
    return a * f3 / (l2l2 - l1l1) - np.arcsin(a / l2)

def common_nl_sigxy_sig2R(x, y, z, a, G, nu, pz):
    pi = pi_()
    c1 = pz / pi

    l1l1 = l1_l1_(x, y, z, a)
    f3 = f3_(x, y, z, a)  # sqrt(l2l2 - a*a)
    c2 = a * l1l1 * f3

    l2l2 = l2_l2_(x, y, z, a)
    c3 = l2l2 * (l2l2 - l1l1)
    r = r_(x, y)
    c4 = z * z * (a * a * (6. * l2l2 - 2. * l1l1 + r * r) - 5. * l2l2 * l2l2)
    c5 = (l2l2 - l1l1) * (l2l2 - l1l1) * (l2l2 - a * a)

    return c1 * c2 / c3 * (1. - 2. * nu + c4 / c5)

def common_nl_sigzx_sigzy(x, y, z, a, G, nu, pz):
    pi = pi_()
    c1 = - 2. * pz / pi
    l1l1 = l1_l1_(x, y, z, a)
    l1 = l1_(x, y, z, a)
    l2 = l2_(x, y, z, a)
    f3 = f3_(x, y, z, a)  # sqrt(l2l2 - a*a)
    l2l2 = l2_l2_(x, y, z, a)
    r = r_(x, y)
    c2 = z * l1 * f3 * (a * a * (4. * l2l2 - 5. * r * r) + l1l1 * l1l1)

    c3 = l2 * (l2l2 - l1l1) * (l2l2 - l1l1) * (l2l2 - l1l1)
    return c1 * c2 / c3

def sig_1_nl_(x, y, z, a, G, nu, pz):
    pi = pi_()
    c1 = 2. * pz / pi
    c2 = 1 + 2 * nu
    c3 = f8_(x, y, z, a)
    l1l1 = l1_l1_(x, y, z, a)
    l2l2 = l2_l2_(x, y, z, a)
    f3 = f3_(x, y, z, a)  # sqrt(l2l2 - a*a)
    r = r_(x, y)
    aa = a * a
    zz = z * z
    rr = r * r
    c4 = a * zz * (l1l1 * l1l1 + aa * (2. * aa + 2. * zz - 3. * rr))
    c5 = (l2l2 - l1l1) * (l2l2 - l1l1) * (l2l2 - l1l1) * f3
    return c1 * (c2 * c3 + c4 / c5)


def sig_2R_nl_(x, y, z, a, G, nu, pz):
    theta = theta_(x, y)
    return 2 * np.cos(2 * theta) * common_nl_sigxy_sig2R(x, y, z, a, G, nu, pz)

def sig_xx_nl_(x, y, z, a, G, nu, pz):
    return 0.5 * (sig_1_nl_(x, y, z, a, G, nu, pz) + sig_2R_nl_(x, y, z, a, G, nu, pz))

def sig_yy_nl_(x, y, z, a, G, nu, pz):
    return 0.5 * (sig_1_nl_(x, y, z, a, G, nu, pz) - sig_2R_nl_(x, y, z, a, G, nu, pz))

def sig_xy_nl_(x, y, z, a, G, nu, pz):
    theta = theta_(x, y)
    return np.sin(2 * theta) * common_nl_sigxy_sig2R(x, y, z, a, G, nu, pz)

def sig_zy_nl_(x, y, z, a, G, nu, pz):
    theta = theta_(x, y)
    return np.sin(theta) * common_nl_sigzx_sigzy(x, y, z, a, G, nu, pz)

def sig_zx_nl_(x, y, z, a, G, nu, pz):
    theta = theta_(x, y)
    return np.cos(theta) * common_nl_sigzx_sigzy(x, y, z, a, G, nu, pz)

def sig_zz_nl_(x, y, z, a, G, nu, pz):
    pi = pi_()
    c1 = 2. * pz / pi
    c3 = f8_(x, y, z, a)
    l1l1 = l1_l1_(x, y, z, a)
    l2l2 = l2_l2_(x, y, z, a)
    f3 = f3_(x, y, z, a)  # sqrt(l2l2 - a*a)
    r = r_(x, y)
    aa = a * a
    zz = z * z
    rr = r * r
    c4 = a * zz * (l1l1 * l1l1 + aa * (2. * aa + 2. * zz - 3. * rr))
    c5 = (l2l2 - l1l1) * (l2l2 - l1l1) * (l2l2 - l1l1) * f3

    return c1 * (c3 - c4 / c5)



    # /*
    #  *
    #  * STRESSES DUE TO AN UNIFORM SHEAR LOADING px and py OVER A CIRCULAR CRACK OF RADIUS a
    #  * expressions valid for z >=0
    #  *
    #  */


def f9_(x, y, z, a):
    l1 = l1_(x, y, z, a)
    r = r_(x, y)
    f3 = f3_(x, y, z, a)  # sqrt(l2l2 - a*a)
    l2l2 = l2_l2_(x, y, z, a)
    l1l1 = l1_l1_(x, y, z, a)
    l2 = l2_(x, y, z, a)
    c1 = z * l1 * f3 * (a * a * (4. * l2l2 - 5. * r * r) + l1l1 * l1l1)
    c2 = l2 * (l2l2 - l1l1) * (l2l2 - l1l1) * (l2l2 - l1l1)
    return c1 / c2

def f10_(x, y, z, a):
    r = r_(x, y)
    f4 = f4_(x, y, z, a)  # sqrt(a*a - l1l1)
    l2l2 = l2_l2_(x, y, z, a)
    l1l1 = l1_l1_(x, y, z, a)
    aa = a * a
    c1 = z * f4 * (l1l1 * l1l1 + aa * (2. * aa + 2. * z * z - 3. * r * r))
    c2 = (l2l2 - l1l1) * (l2l2 - l1l1) * (l2l2 - l1l1)
    return c1 / c2

def f11_(x, y, z, a, nu):
    r = r_(x, y)
    aa = a * a
    rr = r * r
    l2l2 = l2_l2_(x, y, z, a)
    l1l1 = l1_l1_(x, y, z, a)
    f3 = f3_(x, y, z, a)  # sqrt(l2l2 - a*a)
    f4 = f4_(x, y, z, a)  # sqrt(a*a - l1l1)
    c1 = nu * a * f3
    c2 = z * f4 * (aa * (6. * l2l2 - 2 * l1l1 + rr) - 5. * l2l2 * l2l2)
    c3 = (l2l2 - l1l1) * (l2l2 - l1l1)
    return c1 + c2 / c3

def f12_(x, y, z, a):
    l2 = l2_(x, y, z, a)
    l1 = l1_(x, y, z, a)
    l2l2 = l2_l2_(x, y, z, a)
    l1l1 = l1_l1_(x, y, z, a)
    f4 = f4_(x, y, z, a)  # sqrt(a*a - l1l1)
    return  a * l1 * f4 / (l2 * (l2l2 - l1l1))
 

def f13_(x, y, z, a):
    r = r_(x, y)
    f3 = f3_(x, y, z, a)  # sqrt(l2l2 - a*a)
    l2l2 = l2_l2_(x, y, z, a)
    l1l1 = l1_l1_(x, y, z, a)
    l2 = l2_(x, y, z, a)
    l1 = l1_(x, y, z, a)
    c1 = a * a * (4. * l2l2 - 5. * r * r) + l1l1 * l1l1
    c2 = (l2l2 - l1l1) * (l2l2 - l1l1) * (l2l2 - l1l1)
    return z * l1 * f3 * c1 / (l2 * c2)

def f14_(x, y, z, a):
    l2l2 = l2_l2_(x, y, z, a)
    l1l1 = l1_l1_(x, y, z, a)
    l1 = l1_(x, y, z, a)
    l2 = l2_(x, y, z, a)
    f3 = f3_(x, y, z, a)  # sqrt(l2l2 - a*a)
    return 4. * z * a * a * l1 * f3 / (l2l2 * l2 * (l2l2 - l1l1))


def sig_1_sl_(x, y, z, a, G, nu, px, py):
    pi = pi_()
    theta = theta_(x, y)
    c1 = 4. * (px * np.cos(theta) + py * np.sin(theta)) / (pi * (2. - nu))
    c2 = f12_(x, y, z, a)
    c3 = f13_(x, y, z, a)
    return c1 * (-2. * (1. + nu) * c2 + c3)

def sig_2R_sl_(x, y, z, a, G, nu, px, py):
    pi = pi_()
    theta = theta_(x, y)
    c1 = -2./(pi * (2. - nu))
    c2 = 4. * (1. - nu) * f12_(x, y, z, a)
    c3 = f13_(x, y, z, a)
    c4 = px * np.cos(theta) - py * np.sin(theta)
    c5 = f14_(x, y, z, a)
    c6 = px * np.cos(3. * theta) + py * np.sin(3. * theta)

    return c1 * ((c2 - c3) * c4 + (c5 - c3) * c6)

def sig_xx_sl_(x, y, z, a, G, nu, px, py):
    return 0.5 * (sig_1_sl_(x, y, z, a, G, nu, px, py) + sig_2R_sl_(x, y, z, a, G, nu, px, py))

def sig_yy_sl_(x, y, z, a, G, nu, px, py):
    return 0.5 * (sig_1_sl_(x, y, z, a, G, nu, px, py) - sig_2R_sl_(x, y, z, a, G, nu, px, py))

def sig_xy_sl_(x, y, z, a, G, nu, px, py):
    pi = pi_()
    theta = theta_(x, y)
    c1 = -1./(pi * (2. - nu))
    c2 = 4. * (1. - nu) * f12_(x, y, z, a)
    c3 = f13_(x, y, z, a)
    c4 = px * np.sin(theta) + py * np.cos(theta)
    c5 = f14_(x, y, z, a)
    c6 = px * np.sin(3. * theta) - py * np.cos(3. * theta)

    return c1 * ((c2 - c3) * c4 + (c5 - c3) * c6)

 
def sig_zy_sl_(x, y, z, a, G, nu, px, py):
    c1 = 2. / (pi_() * (2. - nu))
    theta = theta_(x, y)
    c2 = f8_(x, y, z, a)
    c3 = f10_(x, y, z, a)
    c4 = f11_(x, y, z, a, nu)
    l2l2 = l2_l2_(x, y, z, a)
    l1l1 = l1_l1_(x, y, z, a)
    return c1 * ( ((2. - nu) * c2 + c3) * py + c4 * l1l1 * (px * np.sin(2. * theta) - py * np.cos(2. * theta)) / (l2l2 * (l2l2 - l1l1)) )

def sig_zx_sl_(x, y, z, a, G, nu, px, py):
    c1 = 2. / (pi_() * (2. - nu))
    theta = theta_(x, y)
    c2 = f8_(x, y, z, a)
    c3 = f10_(x, y, z, a)
    c4 = f11_(x, y, z, a, nu)
    l2l2 = l2_l2_(x, y, z, a)
    l1l1 = l1_l1_(x, y, z, a)
    return c1 * ( ((2. - nu) * c2 + c3) * px + c4 * l1l1 * (px * np.cos(2. * theta) + py * np.sin(2. * theta)) / (l2l2 * (l2l2 - l1l1)) )


def sig_zz_sl_(x, y, z, a, G, nu, px, py):
    theta = theta_(x, y)
    c1 = -4. * (px * np.cos(theta) + py * np.sin(theta)) / (pi_() * (2. - nu))
    return c1 * f9_(x, y, z, a)

    # // ------------------------------------------------------------------------------------
    # /*
    #  *
    #  * DISPLACEMENTS
    #  * - DUE TO AN UNIFORM SHEAR LOADING px and py OVER A CIRCULAR CRACK OF RADIUS a
    #  * - DUE TO AN UNIFORM SHEAR LOADING pz OVER A CIRCULAR CRACK OF RADIUS a
    #  * expressions valid for any z
    #  *
    #  */

    # // normal

def Ux_nl_(x, y, z, a, G, nu, pz):
    return ux_nl_(x, y, z, a, G, nu, pz)

def Uy_nl_(x, y, z, a, G, nu, pz):
    return uy_nl_(x, y, z, a, G, nu, pz)

def Uz_nl_(x, y, z, a, G, nu, pz):
    return uz_nl_(x, y, z, a, G, nu, pz)

def Ux_sl_(x, y, z, a, G, nu, px, py):
    value, abs_z = 0, abs(z)
    if abs_z > 0:
        value = ux_sl_(x, y, abs_z, a, G, nu, px, py)
    if z < 0:
        value = -value
    return value

def Uy_sl_(x, y, z, a, G, nu, px, py):
    value, abs_z = 0, abs(z)
    if abs_z > 0:
        value = uy_sl_(x, y, abs_z, a, G, nu, px, py)
    if z < 0:
        value = -value
    return value

def Uz_sl_(x, y, z, a, G, nu, px, py):
    return uz_sl_(x, y, z, a, G, nu, px, py)


    # /*
    #  *
    #  * STRESSES
    #  * - DUE TO AN UNIFORM SHEAR LOADING px and py OVER A CIRCULAR CRACK OF RADIUS a
    #  * - DUE TO AN UNIFORM SHEAR LOADING pz OVER A CIRCULAR CRACK OF RADIUS a
    #  * expressions valid for any z
    #  *
    #  */
    
def Sig_xx_nl_(x, y, z, a, G, nu, pz):
    return sig_xx_nl_(x, y, z, a, G, nu, pz)

def Sig_yy_nl_(x, y, z, a, G, nu, pz):
    return sig_yy_nl_(x, y, z, a, G, nu, pz)

def Sig_zz_nl_(x, y, z, a, G, nu, pz):
    return sig_zz_nl_(x, y, z, a, G, nu, pz)

def Sig_xy_nl_(x, y, z, a, G, nu, pz):
    return sig_xy_nl_(x, y, z, a, G, nu, pz)

def Sig_zy_nl_(x, y, z, a, G, nu, pz):
    value = 0
    value = sig_zy_nl_(x, y, z, a, G, nu, pz)
    if z < 0 and not ((x > 0 and y > 0) or (x > 0 and y < 0)):
        value = -value
    return value

def Sig_zx_nl_(x, y, z, a, G, nu, pz):
    return sig_zx_nl_(x, y, z, a, G, nu, pz)

# Shear 

def Sig_xx_sl_(x, y, z, a, G, nu, px, py):
    value, abs_z = 0, abs(z)
    if abs_z > 0:
        value = sig_xx_sl_(x, y, abs_z, a, G, nu, px, py)
    if z < 0 and ((x > 0 and y > 0) or (x > 0 and y < 0)):
        value = -value
    return value

def Sig_yy_sl_(x, y, z, a, G, nu, px, py):
    value, abs_z = 0, abs(z)
    if abs_z > 0:
        value = sig_yy_sl_(x, y, abs_z, a, G, nu, px, py)
    if z < 0 and ((x > 0 and y > 0) or (x > 0 and y < 0)):
        value = -value
    return value

def Sig_zz_sl_(x, y, z, a, G, nu, px, py):
    value, abs_z = 0, abs(z)
    if abs_z > 0:
        value = sig_zz_sl_(x, y, abs_z, a, G, nu, px, py)
    if z < 0 and ((x > 0 and y > 0) or (x > 0 and y < 0)):
        value = -value
    return value

def Sig_xy_sl_(x, y, z, a, G, nu, px, py):
    value, abs_z = 0, abs(z)
    if abs_z > 0:
        value = sig_xy_sl_(x, y, abs_z, a, G, nu, px, py)
    if z < 0 and ((x > 0 and y > 0) or (x > 0 and y < 0)):
        value = -value
    return value

def Sig_zy_sl_(x, y, z, a, G, nu, px, py):
    value = sig_zy_sl_(x, y, abs(z), a, G, nu, px, py)
    if z < 0:
        value = -value
    return value

def Sig_zx_sl_(x, y, z, a, G, nu, px, py):
    return sig_zx_sl_(x, y, abs(z), a, G, nu, px, py)




### STRESSES and displacement Solutions  
def stresses_tensile_penny_shaped(x,y,z,a=1.0,sig=1.0,G=1.,nu=0.25):
    """Stresses (in cartesian coordinate) due to a penny-shaped crack under uniform tensile loading
        Args:
        x (_type_): z coordinates
        y (_type_): z coordinates
        z (_type_): z coordinates
        a (float, optional): Crack radius . Defaults to 1.
        sig (float, optional): Applied internal Pressure. Defaults to 1.
        nu (float, optional): Material poisson's ratio. Defaults to 0.25.
        G (float, optional): Material  shear modulus. Defaults to 1. (note the stress solution does not depend on G)

    Returns:
       sxx,syy,szz,sxy,sxz,syz (tuple of np.array): component of the stress tensor in  Cartesian coordinates
    """
    s_xx = Sig_xx_nl_(x,y,z,a,G,nu,sig)
    s_yy = Sig_yy_nl_(x,y,z,a,G,nu,sig)
    s_zz = Sig_zz_nl_(x,y,z,a,G,nu,sig)
    s_xy = Sig_xy_nl_(x,y,z,a,G,nu,sig)
    s_xz = Sig_zx_nl_(x,y,z,a,G,nu,sig)
    s_yz = Sig_zy_nl_(x,y,z,a,G,nu,sig)
    return s_xx,s_yy,s_zz,s_xy,s_xz,s_yz

### STRESSES and displacement Solutions 
def stresses_shear_penny_shaped(x,y,z,a=1.0,T_x=1.0,T_y=1.0,G=1.,nu=0.25):
    """Stresses (in cartesian coordinate) due to a penny-shaped crack under uniform shear loading T_x, T_y
        Args:
        x (_type_): z coordinates
        y (_type_): z coordinates
        z (_type_): z coordinates
        a (float, optional): Crack radius . Defaults to 1.
        T_x (float, optional): Applied internal shear traction along x. Defaults to 1.
        T_y (float, optional): Applied internal shear traction along x. Defaults to 1.
        nu (float, optional): Material poisson's ratio. Defaults to 0.25.
        G (float, optional): Material  shear modulus. Defaults to 1. (note the stress solution does not depend on G)

    Returns:
       sxx,syy,szz,sxy,sxz,syz (tuple of np.array): component of the stress tensor in  Cartesian coordinates
    """
    s_xx = Sig_xx_sl_(x,y,z,a,G,nu,T_x,T_y)
    s_yy = Sig_yy_nl_(x,y,z,a,G,nu,T_x,T_y)
    s_zz = Sig_zz_nl_(x,y,z,a,G,nu,T_x,T_y)
    s_xy = Sig_xy_nl_(x,y,z,a,G,nu,T_x,T_y)
    s_xz = Sig_zx_nl_(x,y,z,a,G,nu,T_x,T_y)
    s_yz = Sig_zy_nl_(x,y,z,a,G,nu,T_x,T_y)
    return s_xx,s_yy,s_zz,s_xy,s_xz,s_yz


def displacement_tensile_penny_shaped(x,y,z,a=1.0,sig=1.0,nu=0.25,G=1.0):
    """Displacement (in Cartesian coordinates) due to a penny-shaped crack under uniform tensile loading

    Args:
        x (_type_): z coordinates
        y (_type_): z coordinates
        z (_type_): z coordinates
        a (float, optional): Crack radius . Defaults to 1.
        sig (float, optional): Applied internal Pressure. Defaults to 1.
        nu (float, optional): Material poisson's ratio. Defaults to 0.25.
        G (float, optional): Material  shear modulus. Defaults to 1.

    Returns:
        u_x, u_y, u_z (tuple of np.array of floats): component of displacements in Cartesian coordinates
    """ 
    u_x = Ux_nl_(x, y, z, a, G, nu, sig)
    u_y = Uy_nl_(x, y, z, a, G, nu, sig)
    u_z = Uz_nl_(x, y, z, a, G, nu, sig)
    return u_x,u_y,u_z


def displacement_shear_penny_shaped(x,y,z,a=1.0,T_x=1.0,T_y=1.0,nu=0.25,G=1.0):
    """Displacement (in Cartesian coordinates) due to a penny-shaped crack under uniform shear loading

    Args:
        x (_type_): z coordinates
        y (_type_): z coordinates
        z (_type_): z coordinates
        a (float, optional): Crack radius . Defaults to 1.
        T_x (float, optional): Applied internal shear traction along x. Defaults to 1.
        T_y (float, optional): Applied internal shear traction along x. Defaults to 1.
        nu (float, optional): Material poisson's ratio. Defaults to 0.25.
        G (float, optional): Material  shear modulus. Defaults to 1.

    Returns:
        u_x, u_y, u_z (tuple of np.array of floats): component of displacements in Cartesian coordinates
    """ 
    u_x = Ux_sl_(x, y, z, a, G, nu, T_x,T_y)
    u_y = Uy_sl_(x, y, z, a, G, nu, T_x,T_y)
    u_z = Uz_sl_(x, y, z, a, G, nu, T_x,T_y)
    return u_x,u_y,u_z





    