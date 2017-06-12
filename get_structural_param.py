import argparse
import graph as gg
from itertools import chain
from parse_pv import parse_file
import subprocess
import os

def main(file_name, trimmed='trimmed.out'):
    # trimmed files
    cmd = 'grep -A16 "STEP NUMBER" {} > {}'.format(file_name, trimmed)
    subprocess.run(cmd, shell=True)    
    rs = parse_file(trimmed)
    
    # remove temporary file   
    os.remove(trimmed)
    xlabel = 'Step'
#    thingsToPlot = ['Pressure', 'Pressure[avg.]','Volume', 'Volume[Avg.]','a', 'b', 'c', 'alpha', 'beta', 'gamma']

    ylabel = 'Pressure'
    title  = 'Pressure [Avg.]'
    output = 'pressur.png'
    gg.simplePlot(rs[:, 0], rs[:, 2], xlabel, ylabel, title,output)
 
    ylabel = 'Volume (Bohr^3)'
    title = 'Volume [Avg.]'
    output = 'volume.png'
    gg.simplePlot(rs[:, 0], rs[:, 4], xlabel, ylabel, title,output)

    formats = ['g-', 'b-', 'r-']
    args1 = createArgs(rs[:, 0], rs[:, 5:8].transpose(), formats)

    ylabel = 'Cell lenghts (Bohr)'
    title = 'Cell lengths [Avg.]'
    output = 'lengths.png'
    gg.plotLengths(xlabel, ylabel, title, output, args1) 

    args2 = createArgs(rs[:, 0], rs[:, 8:11].transpose(), formats)
    ylabel = 'Cell angles (deg)'
    title = 'Cell angles [Avg.]'
    output = 'angles.png'
    gg.plotLengths(xlabel, ylabel, title, output, args2) 

def createArgs( x, ys, formats ):
    """
    generates a list for matplotlib
    [ x, y1, 'format', x, y2, 'format' ]
    
    """
    args = [(x, y, f) for y, f in zip(ys, formats)]
    return list(chain(*args))    


if __name__ == "__main__":
    msg = " get_structural_param -p <path/to/output>"
    parser = argparse.ArgumentParser(description=msg)
    parser.add_argument('-p', required=True, help='path to the cp2k output')
    args = parser.parse_args()
    main(args.p)


