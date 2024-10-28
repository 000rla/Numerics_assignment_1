import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sc

class solver:
    """Numerical solver for the linear advection equation using three different finite difference schemes:
        Forward time backwards space (FTBS)
        Forward time center space (FTCS)
        Center time center space (CTCS)

    Functions:
        plotting : outputs a plot of the numeric solution at five time points
        animation_plotting : outputs an animation of the numeric solution 
        create phi : find the analytical solution for the linear advection equation
        FTBS : solve the equation using the forward time backwards space scheme
        FTCS : solve the equation using the forward time center space scheme
        CTCS : solve the equation using the center time center space scheme

    Args:
        animation (bool, optional): Defaults to False. If True, animation_plotting will run, otherwise plotting will run. 
        u (float, optional): Defaults to 1.0. The wind speed
        nx (int, optional): Defaults to 50. The number of points in space. Keep in mind the Courant number when choosing a value.
        nt (int, optional): Defaults to 250. The number of time steps. Keep in mind the Courant number when choosing a value.
    """

    def __init__(self, plotting=True, animation=False, conservation_check=False, u=1, nx=50, nt=250, x=None):
        self.u = u # The wind speed
        self.nx = int(nx) # number of points in space
        self.x = np.linspace(0.0, 1.0, self.nx+1) # From zero to one inclusive
        self.nt = int(nt) # The number of time steps
        self.dx = 1./self.nx # The spacial resolution
        self.dt = 1./self.nt # The time step
        self.c = self.dt*self.u/self.dx #The Courant number

        print("The Courant number for this experiment is: ",self.c)

        self.phi = self.create_phi(0) #np.where(self.x%1. < 0.5, np.power(np.sin(2*self.x*np.pi), 2), 0.)
        self.phiOld = self.phi.copy()

        self.plot=plotting
        self.animation=animation
        self.conservation_check=conservation_check

        # #plotting analytical solution over time
        # fig, ax = plt.subplots(1,1,figsize=(14,6))
        # for n in [0.0,0.2,0.4,0.6,0.8]:
        #     plt.plot(self.x,self.create_phi(n), label='Time = '+str(n))

        # #making the legend
        # box = ax.get_position()
        # ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9]) # Shrink current axis by 20%
        # ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True, ncol=5) # Put a legend below current axis

        # plt.title('Analytical solution of $\\phi$')
        # plt.ylabel('$\\phi$')
        # plt.xlabel('x')
        # plt.grid(alpha=.5)
        # plt.ylim([0,1])
        # plt.xlim(0,1)

        # plt.savefig('analytical_solution.pdf')
        # plt.show()
        # plt.cla()

    def animation_plotting(self,n,T='title'):
        """Create an animation of the finite difference solution and the analytical solution evolving over time.

        Args:
            n (int): time step the plot is created at.
            T (str, optional): The title given to the plot and to the legend keys. Defaults to 'title'.
        """
        plt.cla()

        plt.plot(self.x,self.create_phi(self.dt*n), linestyle='dashed', label='Analytical')
        plt.plot(self.x, self.phi, label='Time '+str(n*self.dt))

        plt.legend(loc='best',title='c = '+str(self.c))
        plt.title(T)
        plt.ylabel('phi')
        plt.ylim([0,1])
        plt.pause(0.1)

    def plotting(self,n,i,T='title'):
        """Plots the finite difference solution and the analytical soltion together at five different time steps.

        Args:
            n (int): time step the plot is created at.
            i (int): Chooses the colour to be used for the plotting from seaborn package.
            T (str, optional): The title given to the plot and to the legend keys. Defaults to 'title'.
        """
        
        plt.plot(self.x,self.create_phi(self.dt*n), linestyle='dashed', label='Analytical at time '+str(n*self.dt),  c = sns.color_palette('tab10')[i])
        plt.plot(self.x, self.phi, label=T+' at time '+str(n*self.dt), c = sns.color_palette('tab10')[i])

        #making the legend
        box = self.ax.get_position()
        self.ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9]) # Shrink current axis by 20%
        self.ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True, ncol=5) # Put a legend below current axis

        plt.title(T)
        plt.ylabel('$\\phi$')
        plt.xlabel('x')
        plt.grid(alpha=.5)
        plt.ylim([0,1])
        plt.xlim(0,1)
    
    def mass_check(self,n):
        #intergrate using simpsons rule
        mass_anal=sc.integrate.simps(self.create_phi(self.dt*n),self.x)
        mass_fd = sc.integrate.simps(self.phi,self.x)

        #printing results
        print("at time "+str(n*self.dt)+" the volume of the analytical solution was "+str(mass_anal))
        print("the volume of the finite difference solution was "+str(mass_fd))
        print("that's a difference of "+str(abs(mass_anal-mass_fd)))

    def create_phi(self,t):
        """Finds the analytical soltion to the linear advection equation at given time t as given in the Kick-Off Camp notes.

        Args:
            t (int): time

        Returns:
            phi (numpy array): analytical solution to the linear advection equation at time t.
        """
        return np.where((self.x-self.u*t)%1. < 0.5, np.power(np.sin(2*(self.x-self.u*t)*np.pi), 2), 0.)
    
    def convergance_checker(self):
        plt.loglog((np.arange(1,self.nt)),(self.error_FTBS[1:]),label='FTBS')
        plt.loglog((np.arange(1,self.nt)),(self.error_FTCS[1:]),label='FTCS')
        plt.loglog((np.arange(3,self.nt)),(self.error_CTCS[1:]),label='CTCS')
        plt.grid()
        plt.legend()
        plt.show()

    def RMSE(self,aim,pred):
        return np.sqrt(np.mean((aim-pred)**2))
    

    def FTBS(self):
        """Solves the linear advection equation using the forward in time backwards in space scheme and plots the results.
        """
        # print("now solving using FTBS")

        if self.plot:
            # set up plotting
            fig, self.ax = plt.subplots(1,1,figsize=(14,6))
            i=0

        # solve for phi
        for n in range(self.nt): # loop over time
            for j in range(1,self.nx+1): # loop over space
                self.phi[j] = self.phiOld[j]-self.c*(self.phiOld[j]-self.phiOld[j-1])
        # setting periodic boundry
            self.phi[0]=self.phi[-1]

        # update phi for the next time−step
            self.phiOld = self.phi.copy()
            
        # plot solutions
            if self.animation:
                self.animation_plotting(n,'FTBS')
            if self.plot:
                if n%(self.nt//5)==0:
                    if self.conservation_check:
                        self.mass_check(n)
                    self.plotting(n,i,'FTBS')
                    i+=1
        if self.animation:
            plt.show() 
        if self.plot:
            plt.savefig('FTBS_test.pdf')
            plt.show() 

        return self.RMSE(self.create_phi(self.nt*self.dt),self.phi)

    def FTCS(self):
        """Solves the linear advection equation using the forward in time centered in space scheme and plots the results.
        """
        # print("now solving using FTCS")
        
        if self.plot:
            # set up plotting
            fig, self.ax = plt.subplots(1,1,figsize=(14,6))
            i=0

        # solve for phi
        for n in range(self.nt):
            for j in range(1,self.nx): # loop over space
                self.phi[j] = self.phiOld[j]-.5*self.c*(self.phiOld[j+1]-self.phiOld[j-1])

        # setting periodic boundry
            self.phi[0] = self.phiOld[0]-.5*self.c*(self.phiOld[1]-self.phiOld[-2])
            self.phi[-1]=self.phi[0]

        # update phi for the next time−step
            self.phiOld = self.phi.copy()
            
        #plot solutions
            if self.animation:
                self.animation_plotting(n,'FTCS')
            if self.plot:
                if n%(self.nt//5)==0:
                    if self.conservation_check:
                        self.mass_check(n)
                    self.plotting(n,i,'FTCS')
                    i+=1
        if self.animation:
            plt.show() 
        if self.plot:
            plt.savefig('FTCS_test.pdf')
            plt.show() 

        return self.RMSE(self.create_phi(self.nt*self.dt),self.phi)


    def CTCS(self):
        """Solves the linear advection equation using the centered in time centered in space scheme and plots the results.
        """
        # print('Now solving using CTCS')
        if self.plot:
            #set up plotting
            fig, self.ax = plt.subplots(1,1,figsize=(14,6))
            i=0

        # define phi at time 0 and time 1
        phiOlder = self.create_phi(0) #phi at time 0
        self.phiOld = self.create_phi(self.dt) #phi at time 1
        
        # solve for phi
        for n in range(2,self.nt): #loop over time
            for j in range(1,self.nx): # loop over space
                self.phi[j] = phiOlder[j]-self.c*(self.phiOld[j+1]-self.phiOld[j-1])

        # setting periodic boundry
            self.phi[0] = phiOlder[0]-self.c*(self.phiOld[1]-self.phiOld[-2])
            self.phi[-1] = self.phi[0]

        # update phi for the next time−step
            phiOlder = self.phiOld.copy()
            self.phiOld = self.phi.copy()
            
        # plot solutions
            if self.animation:
                self.animation_plotting(n,'CTCS')
            if self.plot:
                if n%(self.nt//5) == 0 or n == 2:
                    if self.conservation_check:
                        self.mass_check(n)
                    self.plotting(n, i, 'CTCS')
                    i+=1
        if self.animation:
            plt.show() 
        if self.plot:
            plt.savefig('CTCS_test.pdf')
            plt.show() 

        return self.RMSE(self.create_phi(self.nt*self.dt),self.phi)