# %%
import sys
import os
home = os.environ["HOME"]
sys.path.append("/Users/bricelecampion/ClionProjects/BigWham_dev/build/interfaces/python")
from hmatrix import Hmatrix
import numpy as np
from scipy.sparse.linalg import gmres

sys.path.append(os.path.join(os.getcwd(), '../..'))
from solutions.griffith_crack import width_griffith


# %% Material properties
G = 1.0
nu = 0.25
E = (2 * G) * (1 + nu)

# %% function encapsulating the solution !

def solveGriffith(nelts,a=1,kernel = "2DP0"):
    """Wrapper function to solve the Griffith problem with a given number of elements + element type"""
    coor1D = np.linspace(-a, a, nelts + 1)
    coor = np.transpose(np.array([coor1D, coor1D * 0.0]))
    conn = np.fromfunction(lambda i, j: i + j, (nelts, 2), dtype=np.int_)
    # H-matrix parameters - choosen here to have an accurate solution
    max_leaf_size=100
    eta=2.
    eps_aca=1.e-5
    elas_prop = np.array([E, nu])
    h = Hmatrix(kernel, coor, conn, elas_prop, max_leaf_size, eta, eps_aca)
    t = np.ones(h.shape[0])
    t[0::2] = 0.
    jac_prec =h.H_jacobi_prec()
    jac_ilu= h.H_ILU_prec()
    d = gmres(h, t,M=jac_ilu,tol=1e-6)[0]
    dd = d.reshape((-1, 2))
    col_pts = h.getMeshCollocationPoints()
    x_coor_=col_pts[:,0]
    # Crack opening displacement discontinuity from analytical solution 
    # note that 2DP1 segment have their solution at nodes.... 
    # note that 2DP1 segment have displacement disconinuity at element vertex for each element
    if (kernel =="2DP1"):
        x_coor_ = np.ones(2*nelts)
        for e in range(nelts):
            x_coor_[e*2:e*2+2]=coor1D[conn[e]]
    
    w_true=width_griffith(x_coor_,a=1,sig=1,G=1,nu=0.25)
    rmse = np.sqrt((np.sum(dd[:,1]-w_true)**2)/(w_true.size))
    l2rel = (np.linalg.norm((dd[:,1]-w_true)))/(np.linalg.norm(w_true))
    return col_pts,dd,w_true,rmse,l2rel 


# %%
col_pts,dd_sol,w_true,rmse,l2rel = solveGriffith(30,a=1,kernel="2DP0")

import matplotlib.pyplot as plt
plt.plot(col_pts[:, 0], dd_sol[:,1], "-*")
plt.plot(col_pts[:, 0], w_true, "-r*")

rmse = np.sqrt((np.sum(dd_sol[:,1]-w_true)**2)/(w_true.size))


# %%
list_nelts=np.array([10,20,50,100,200,500,1000,2000,5000,10000])
rmse =[]
l2_rel = []
for n in list_nelts:
    col_pts,dd_sol,wtr,rmse_n,lerel_n = solveGriffith(n,a=1)
    w_true=width_griffith(col_pts[:, 0],a=1,sig=1,G=1,nu=0.25)
    rmse.append(rmse_n)
    l2_rel.append(lerel_n)
    
# computing the rate of convergence
beta_rmse = [1 * np.log(rmse[i]/rmse[i-1])/np.log(list_nelts[i-1]/list_nelts[i]) for i in range(2,list_nelts.size-1)]
beta_l2 = [1 * np.log(l2_rel[i]/l2_rel[i-1])/np.log(list_nelts[i-1]/list_nelts[i]) for i in range(2,list_nelts.size-1)]

# %%
plt.loglog(list_nelts, rmse, "-*")
plt.loglog(list_nelts, 1/list_nelts**0.5, "-")
plt.title('RMSE vs Nelts')
# rmese defined as  sqrt(1/Nelts Sum  (num_i-true_i)^2)
#%%
plt.loglog(list_nelts, l2_rel, "-*")
plt.loglog(list_nelts, 1/list_nelts, "-")
plt.title('L2 norm relative error vs Nelts')
# L2 norm relative defined as ||num_i - true_i||_2 / ||true_i||
