import numpy as np
import matplotlib
matplotlib.use('tkagg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


def simplePlot(x,y,xlabel,ylabel,title,output):
    """
    this is a simple plot`
    """
    plt.plot(x,y)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(output, format='png', dpi=300)
    plt.close()
#    plt.show()


def plotLengths(xlabel,ylabel,title,output,args):
    """
    plots multiple data into one single plot
    """
    plt.plot(*args)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(output, format='png', dpi=300)
    plt.close()
