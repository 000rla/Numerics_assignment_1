import functions
import matplotlib.pyplot as plt
import numpy as np

import scipy.optimize as spo


def main():
    """runs the finite difference solvers for the linear advection equation
    """
    f=functions.solver()
    f.FTBS()
    f.FTCS()
    f.CTCS()

def main_animations():
    """runs like main, but outputs animations of the solutions instead of a plot.
    """
    f=functions.solver(animation=True)
    # f.FTBS()
    f.FTCS()
    # f.CTCS()

def linear_function(x,m,c):
    return m*x+c

def convergance_experiment():
    
    nx=np.array([50,100,200])
    nt=5*nx
    ftbs_error=np.zeros(len(nx))
    #ftcs_error=np.zeros(len(nx))
    ctcs_error=np.zeros(len(nx))
    
    for i in range(len(nx)):
        f=functions.solver(plotting=False, nx=nx[i], nt=nt[i])
        ftbs_error[i]=f.FTBS()
        #ftcs_error[i]=f.FTCS()
        ctcs_error[i]=f.CTCS()

    [m_ftbs_sp, c_ftbs_sp], pcov_ftbs_fit = spo.curve_fit(linear_function, np.log(1/nx), np.log(ftbs_error))
    [m_ctcs_sp, c_ctcs_sp], pcov_ctcs_fit = spo.curve_fit(linear_function, np.log(1/nx), np.log(ctcs_error))
    
    plt.loglog(1/nx, ftbs_error, '--bo', label='FTBS')
    plt.loglog(1/nx, 15*m_ftbs_sp/nx, 'b', label='Scipy fit: m ='+str(round(m_ftbs_sp,3)))
    plt.loglog(1/nx, 15*1/nx, '--c', label='A$\\Delta x^2$')
    
    plt.loglog(1/nx, ctcs_error, '--ro', label='CTCS')
    plt.loglog(1/nx, .1*m_ctcs_sp/nx, 'r', label='Scipy fit: m ='+str(round(m_ctcs_sp,3)))
    plt.loglog(1/nx, .1*2/nx,'--m', label='2A$\\Delta x^2$')

    plt.xlabel('$\\Delta x$')
    plt.ylabel('$l_2$ error')

    plt.grid()
    plt.legend()
    plt.savefig('convergance_experiment.pdf')
    plt.show()
    
convergance_experiment()