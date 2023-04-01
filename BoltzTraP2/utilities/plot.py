from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

def plot():
    %matplotlib inline
    fig = plt.figure()
    ax = fig.gca(projection='3d')
 
    data = np.genfromtxt('thermo.out')
    #zero coupon maturity dates
    y = data[:,0]
    #tenor
    x = data[:,1]
    #rates
    z = data[:,2]

    # plotting 
    ax.plot3D(x, y, z, 'green')
    #ax.plot_surface(x, y, z, rstride = 2, cstride = 1, cmap = plt.cm.Blues_r) 
    ax.set_title('Chemical potential')
    ax.set_ylabel('phase fraction')
    ax.set_xlabel('T (K)')
    ax.set_zlabel('$\mu$ (eV)')
 
    plt.show() 