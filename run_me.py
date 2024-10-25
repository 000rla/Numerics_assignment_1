import functions

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

def convergance_experiment():
    import matplotlib.pyplot as plt
    import numpy as np

    nx=[50,100,150]
    nt=[250,500,750]
    ftbs_error=[]
    ftcs_error=[]
    ctcs_error=[]

    f=functions.solver(plotting=True, nx=nx[0], nt=nt[0])
    err1=f.FTBS()
    ftbs_error.append(err1)
    err2=f.FTCS()
    ftcs_error.append(err2)
    err3=f.CTCS()
    ctcs_error.append(err3)

    g=functions.solver(plotting=False, nx=nx[1], nt=nt[1])
    err1=g.FTBS()
    ftbs_error.append(err1)
    err2=f.FTCS()
    ftcs_error.append(err2)
    err3=f.CTCS()
    ctcs_error.append(err3)

    h=functions.solver(plotting=False, nx=nx[2], nt=nt[2])
    err1=h.FTBS()
    ftbs_error.append(err1)
    err2=f.FTCS()
    ftcs_error.append(err2)
    err3=f.CTCS()
    ctcs_error.append(err3)

    fig = plt.figure()
    ax = plt.gca()
    
    plt.scatter(nx,ftbs_error,label='FTBS',c='b')
    plt.scatter(nx,ftcs_error,label='FTCS',c='r')
    plt.scatter(nx,ctcs_error,label='CTCS',c='g')
    plt.loglog([nx[0],nx[-1]],[ftbs_error[0],ftbs_error[-1]],c='c',linestyle='dashed',label='FTBS')
    plt.loglog([nx[0],nx[-1]],[ftcs_error[0],ftcs_error[-1]],c='m',linestyle='dashed',label='FTCS')
    plt.loglog([nx[0],nx[-1]],[ctcs_error[0],ctcs_error[-1]],c='olive',linestyle='dashed',label='CTCS')

    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.grid()
    plt.legend()
    plt.show()
    
convergance_experiment()