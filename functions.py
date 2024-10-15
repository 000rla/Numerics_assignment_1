import numpy as np
import matplotlib.pyplot as plt

N=100

def finite_differences():
    u = 1.0 # The wind speed
    nx = 40 # number of points in space
    x = np.linspace(0.0, 1.0, nx+1) # From zero to one inclusive
    nt = 80 # The number of time steps
    dx = 1./nx # The spacial resolution
    dt = 1./nt # The time step
    # The initial conditions and arrays for the old and new time steps
    phi = np.where(x%1. < 0.5, np.power(np.sin(2*x*np.pi), 2), 0.)
    phiOld = phi.copy()
    # Plot the initial conditions
    plt.plot(x, phi, 'k', label='initial conditions')
    plt.legend(loc='best')
    plt.ylabel('$\phi$')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([0,1])
    plt.pause(1)
    # Loop over all time steps
    for n in range(nt):
        for j in range(1,nx): # loop over space from 1 to nx−1
            # (avoiding boundary conditions)
            phi[j] = phiOld[j] - u*dt* ...
            # apply boundary conditions of your choice
            #
            # update phi for the next time−step
            phiOld = phi.copy()
            # Replot
            plt.cla()
            plt.plot(x, phi, 'b', label='Time '+str(n*dt))
            plt.legend(loc='best')
            plt.ylabel('$\phi$')
            plt.ylim([0,1])
            plt.pause(0.05)
    plt.show() # To keep the plot showing at the en