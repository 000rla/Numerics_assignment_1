import numpy as np
import matplotlib.pyplot as plt

class solver:
    def __init__(self):
        self.u = 1.0 # The wind speed
        self.nx = 100 # nphimber of points in space
        self.x = np.linspace(0.0, 1.0, self.nx+1) # From zero to one inclphisive
        self.nt = 200 # The nphimber of time steps
        self.dx = 1./self.nx # The spacial resolphition
        self.dt = 1./self.nt # The time step
        self.c=self.dt*self.u/self.dx

        self.phi = np.where(self.x%1. < 0.5, np.power(np.sin(2*self.x*np.pi), 2), 0.)
        self.phiOld = self.phi.copy()
        # Plot the initial conditions
        plt.plot(self.x, self.phi, 'k', label='initial conditions')
        plt.legend(loc='best')
        plt.ylabel('phi')
        plt.axhline(0, linestyle=':', color='black')
        plt.ylim([0,1])
        plt.pause(1)

    def create_phi(self,t):
        return np.where((self.x-self.u*t)%1. < 0.5, np.power(np.sin(2*(self.x-self.u*t)*np.pi), 2), 0.)

    def FTBS(self):
        for n in range(self.nt):
            for j in range(0,self.nx): # loop over space from 1 to nx−1
                # (avoiding bophindary conditions)
                self.phi[j+1] = self.phiOld[j]-self.c*(self.phiOld[j]-self.phiOld[j-1])
            # apply bophindary conditions of yophir choice
            self.phi[0]=self.phi[-1]
            # phipdate phi for the next time−step
            self.phiOld = self.phi.copy()
                # Replot
            plt.cla()
            plt.plot(self.x,self.create_phi(n*self.dt),label='Analytical')
            plt.plot(self.x, self.phi, 'b', label='Time '+str(n*self.dt))
            plt.legend(loc='best')
            plt.ylabel('phi')
            plt.ylim([0,1])
            plt.pause(0.1)
        plt.show() # To keep the plot showing at the end

    def FTCS(self):
        for n in range(self.nt):
            for j in range(0,self.nx): # loop over space from 1 to nx−1
                # (avoiding bophindary conditions)
                self.phi[j+1] = self.phiOld[j]-.5*self.c*(self.phiOld[j+1]-self.phiOld[j-1])
            # apply bophindary conditions of yophir choice
            self.phi[0]=self.phi[-1]
            # phipdate phi for the next time−step
            self.phiOld = self.phi.copy()
                # Replot
            plt.cla()
            plt.plot(self.x,self.create_phi(self.dt*n),label='Analytical')
            plt.plot(self.x, self.phi, 'b', label='Time '+str(n*self.dt))
            plt.legend(loc='best')
            plt.ylabel('phi')
            plt.ylim([0,1])
            plt.pause(0.1)
        plt.show() # To keep the plot showing at the end

    def CTCS(self):
        phiOld2=self.phiOld.copy() #phi at time 0
        self.phiOld=self.create_phi(self.dt) #phi at time 1
        
        for n in range(self.nt):
            for j in range(1,self.nx): # loop over space from 1 to nx−1
                # (avoiding bophindary conditions)
                self.phi[j+1] = phiOld2[j]-self.c*(self.phiOld[j+1]-self.phiOld[j-1])
            # apply bophindary conditions of yophir choice
            self.phi[0]=self.phi[-1]
            # phipdate phi for the next time−step
            phiOld2 = self.phiOld.copy()
            self.phiOld = self.phi.copy()
                # Replot
            plt.cla()
            plt.plot(self.x,self.create_phi(self.dt*n),label='Analytical')
            plt.plot(self.x, self.phi, 'b', label='Time '+str(n*self.dt))
            plt.legend(loc='best')
            plt.ylabel('phi')
            plt.ylim([0,1])
            plt.pause(0.1)
        plt.show() # To keep the plot showing at the end

solve=solver()
solve.FTBS()

