import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class solver:
    def __init__(self,animation=False):
        self.u = 1.0 # The wind speed
        self.nx = 50 # number of points in space
        self.x = np.linspace(0.0, 1.0, self.nx+1) # From zero to one inclusive
        self.nt = 250 # The number of time steps
        self.dx = 1./self.nx # The spacial resolution
        self.dt = 1./self.nt # The time step
        self.c = self.dt*self.u/self.dx #The Corant number

        self.phi = self.create_phi(0) #np.where(self.x%1. < 0.5, np.power(np.sin(2*self.x*np.pi), 2), 0.)
        self.phiOld = self.phi.copy()

        self.animation=animation

    def animation_plotting(self,n):
        plt.cla()
        plt.plot(self.x,self.create_phi(self.dt*n), linestyle='dashed', label='Analytical')
        plt.plot(self.x, self.phi, label='Time '+str(n*self.dt))
        plt.legend(loc='best')
        plt.title('c = '+str(self.c))
        plt.ylabel('phi')
        plt.ylim([0,1])
        plt.pause(0.1)

    def plotting(self,n,i):
        plt.plot(self.x,self.create_phi(self.dt*n), linestyle='dashed', label='Analytical',  c = sns.color_palette('tab10')[i])
        plt.plot(self.x, self.phi, label='Time '+str(n*self.dt), c = sns.color_palette('tab10')[i])
        plt.legend(loc='best')
        plt.title('c = '+str(self.c))
        plt.ylabel('phi')
        plt.ylim([0,1])

    def create_phi(self,t):
        return np.where((self.x-self.u*t)%1. < 0.5, np.power(np.sin(2*(self.x-self.u*t)*np.pi), 2), 0.)

    def FTBS(self):
        i=0
        for n in range(self.nt):
            for j in range(1,self.nx+1): # loop over space
                # (avoiding bophindary conditions)
                self.phi[j] = self.phiOld[j]-self.c*(self.phiOld[j]-self.phiOld[j-1])
            # apply bophindary conditions of yophir choice
            self.phi[0]=self.phi[-1]
            # phipdate phi for the next time−step
            self.phiOld = self.phi.copy()
            
            #plot
            if self.animation:
                self.animation_plotting(n)
            else:
                if n%(self.nt//5)==0:
                    self.plotting(n,i)
                    i+=1
        if self.animation:
            plt.show() 
        else:
            plt.savefig('FTBS_test.jpg')
            plt.show() 

    def FTCS(self):
        i=0
        for n in range(self.nt):
            for j in range(1,self.nx): # loop over space
                # (avoiding bophindary conditions)
                self.phi[j] = self.phiOld[j]-.5*self.c*(self.phiOld[j+1]-self.phiOld[j-1])
            # apply bophindary conditions of yophir choice
            self.phi[0] = self.phiOld[0]-.5*self.c*(self.phiOld[1]-self.phiOld[-2])
            self.phi[-1]=self.phi[0]
            # phipdate phi for the next time−step
            self.phiOld = self.phi.copy()
            
            #plot
            if self.animation:
                self.animation_plotting(n)
            else:
                if n%(self.nt//5)==0:
                    self.plotting(n,i)
                    i+=1
        if self.animation:
            plt.show() 
        else:
            plt.savefig('FTBS_test.jpg')
            plt.show() 


    def CTCS(self):
        i=0
        phiOlder=self.phiOld.copy() #phi at time 0
        self.phiOld=self.create_phi(self.dt) #phi at time 1
        
        for n in range(1,self.nt):
            for j in range(1,self.nx): # loop over space 
                # (avoiding bophindary conditions)
                self.phi[j] = phiOlder[j]-self.c*(self.phiOld[j+1]-self.phiOld[j-1])
            # apply bophindary conditions of yophir choice
            self.phi[0] = phiOlder[0]-self.c*(self.phiOld[1]-self.phiOld[-2])
            self.phi[-1] = self.phi[0]
            # phipdate phi for the next time−step
            phiOlder = self.phiOld.copy()
            self.phiOld = self.phi.copy()
            
            #plot
            if self.animation:
                self.animation_plotting(n)
            else:
                if n%(self.nt//5)==0:
                    self.plotting(n,i)
                    i+=1
        if self.animation:
            plt.show() 
        else:
            plt.savefig('FTBS_test.jpg')
            plt.show() 


solve=solver(animation=False)
solve.FTBS()
solve.FTCS()
solve.CTCS()

