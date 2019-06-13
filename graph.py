import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
import matplotlib
matplotlib.use('tkagg')

"""
Simple plots to be called in other scripts
"""


def simplePlot(x, y, xlabel, ylabel, title, output):
    """
    This is a simple plot`
    """
    plt.plot(x, y)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(output, format='png', dpi=300)
    plt.close()
#    plt.show()


def plotLengths(xlabel, ylabel, title, output, args):
    """
    Plots multiple data into one single plot
    """
    plt.plot(*args)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(output, format='png', dpi=300)
    plt.close()
