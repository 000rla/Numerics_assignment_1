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
    f.FTBS()
    f.FTCS()
    f.CTCS()
    
main()