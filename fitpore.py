# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 22:30:04 2011

@author: Pavlo Shchelokovskyy
"""
import sys
import argparse

from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import linregress

parser = argparse.ArgumentParser(description='Get pore tension in 2 steps.',
        epilog="""First run with only a filename as input and remember the boundaries of the linear stage.
        Then run again supplying all arguments to get the linear region fitted and pore tension displayed.""")
parser.add_argument('filename', help='Name of the input text file with pore radii and frame numbers')
parser.add_argument('-s', type=int, default=0, help='Start of the linear region (frame number)')
parser.add_argument('-e', type=int, default=0, help='End of the linear region (frame number)')
parser.add_argument('-r', type=float, help='Initial radius of the vesicle in microns.')
parser.add_argument('-f', type=float, help='Speed of image acquisition in frames per second')
parser.add_argument('-v', type=float, default=1.0e-3, 
                    help='Viscosity of the bulk solution in Pa*s (defaults 1e-3 Pa*s for water)')

if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
args = parser.parse_args()

radius, frame = np.loadtxt(args.filename, unpack=1, usecols=(0,5))
lnr = np.log(radius)

if not (args.r and args.f):
    plt.plot(frame, lnr, 'ro', label='experiment')
    plt.title("Choose the linear part and supply its range in the next run")
    plt.legend()
    plt.xlabel('frame number')
    plt.ylabel('$\\ln (r_p)$')
    plt.show()
else:
    
    lim1, lim2 = args.s, args.e
    if not lim1:
        ind1 = 0
    else:
        ind1, = np.where(frame == lim1)
        if not ind1:
            print "No frame number %i.\nExiting..."%lim1
            sys.exit(1)

    if not lim2:
        x = frame[ind1:]
        y = lnr[ind1:]
    else:
        ind2, = np.where(frame == lim2)
        if not ind2:
            print "No frame number %i.\nExiting..."%lim2
            sys.exit(1)
        if ind1 > ind2:
            ind1, ind2 = ind2, ind1
        x = frame[ind1:ind2+1]
        y = lnr[ind1:ind2+1]
    
    if len(x) < 2:
        print "Too small linear region, can not fit.\nExiting..."
        sys.exit(1)
        
    a, b, r, p, std = linregress(x, y)

    plt.plot(frame, lnr, 'ro', label='experiment')
    plt.plot(frame, b + a*frame, 'g-', label='linear fit')
    plt.axvline(lim1, c='b', ls='--')
    plt.axvline(lim2, c='b', ls='--')
    plt.legend()
    plt.xlabel('frame number')
    plt.ylabel('$\\ln (r_p)$')
    
    modelgamma = lambda x: -1.5 * np.pi * args.v * args.r*args.r * args.f *x
    gamma =  modelgamma(a)
    gammastd = np.abs(modelgamma(std))
    
    #since R is in micrometers and there is R**2, gamma is in picoNewtons
    plt.title('$\\gamma$ = %f $\\pm$ %f pN, from %i to %i'%(gamma, gammastd, lim1, lim2))
    plt.show()
    

