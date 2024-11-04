import functions
import matplotlib.pyplot as plt
import numpy as np

import scipy.optimize as spo


def main():
    """runs the finite difference solvers for the linear advection equation
    """
    f=functions.solver(conservation_check=True)
    f.FTBS()
    f.FTCS()
    f.CTCS()

def main_animations():
    """runs like main, but outputs animations of the solutions instead of a plot.
    """
    f=functions.solver(animation=True)
    f.FTBS()
    f.FTCS()
    f.CTCS()

def convergance_experiment():
    
    nx=np.array([50,100,200])
    nt=5*nx
    ftbs_error=np.zeros(len(nx))
    #ftcs_error=np.zeros(len(nx))
    ctcs_error=np.zeros(len(nx))
    
    for i in range(len(nx)):
        f=functions.solver(plotting=False, convergence_experiment=True, nx=nx[i], nt=nt[i])
        ftbs_error[i]=f.FTBS()
        #ftcs_error[i]=f.FTCS()
        ctcs_error[i]=f.CTCS()

    [m_ftbs_sp, c_ftbs_sp], pcov_ftbs_fit = spo.curve_fit(functions.linear_function, np.log(1/nx), np.log(ftbs_error))
    [m_ctcs_sp, c_ctcs_sp], pcov_ctcs_fit = spo.curve_fit(functions.linear_function, np.log(1/nx), np.log(ctcs_error))
    
    plt.loglog(1/nx, ftbs_error, '--bo', label='FTBS')
    plt.loglog(1/nx, 10/nx**m_ftbs_sp, 'b', label='Scipy fit: m ='+str(round(m_ftbs_sp,3)))
    plt.loglog(1/nx, 14/nx, '--c', label='$A\\Delta x$')

    plt.loglog(1/nx, ctcs_error, '--ro', label='CTCS')
    plt.loglog(1/nx, 45/nx**m_ctcs_sp, 'r', label='Scipy fit: m ='+str(round(m_ctcs_sp,3)))
    plt.loglog(1/nx, 45/nx**2,'--m', label='$A\\Delta x^2$')

    plt.xlabel('$\\Delta x$')
    plt.ylabel('$l_2$ error')

    plt.grid(which='both',alpha=0.5)
    plt.legend()
    plt.savefig('convergance_experiment.pdf')
    plt.show()


main()    
convergance_experiment()