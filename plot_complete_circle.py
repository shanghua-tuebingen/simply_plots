
import shutil
import datetime
import os, errno
import csv
import numpy as np
import random
import subprocess
#import pandas
import math
import sys
import re
import pylab as plt
import scipy.stats as st
import datetime
import pandas as pd
from matplotlib import gridspec
import matplotlib.pylab as py
import matplotlib
# matplotlib.use("Agg")
from matplotlib import colors
# from matplotlib import colors
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.animation as manimation
# from matplotlib import pyplot as plt
# import cv2
from matplotlib.ticker import FormatStrFormatter
from scipy.stats.stats import pearsonr
from scipy.stats import linregress
import re
from scipy.optimize import curve_fit
import sympy as sym
from decimal import Decimal
from matplotlib.patches import Circle
import matplotlib.cm as cm
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition
import matplotlib.pyplot as plt; plt.rcdefaults()
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
from scipy.stats import sem, t

path_out=r'D:\Github_code'


def main():
    test()
    plot_test()
    
def test():

    circle_radius=np.zeros((5,3))
    circle_radius[:,0]=np.random.randint(0,5,5)
    circle_radius[:, 1] = np.random.randint(5,10,5)
    circle_radius[:, 2] = np.random.uniform(1,3,5)

    np.savetxt(os.path.join(path_out, 'circle_radius.txt'),circle_radius)

    print(circle_radius)
    # sys.exit()
   # find the circles which x+r oy+r larger than the 2000 cells
    radius=np.loadtxt(os.path.join(path_out,'circle_radius.txt'))

    line_radius= sum(1 for line in open(os.path.join(path_out,'circle_radius.txt')))

    data_r_extra=np.zeros((line_radius,3))
    for i in range (0,line_radius):#left, right, top, bottom [:,0]=5,[:,1]=5
        if (radius[i,0]+radius[i,2]>5)==True or (radius[i,0]-radius[i,2]>5)==True or (radius[i,1]+radius[i,2]>10)==True or (radius[i,1]-radius[i,2]>10)==True\
            or (radius[i, 0] + radius[i, 2] < 0) == True or (radius[i, 0] - radius[i, 2] < 0) == True or (radius[i, 1] + radius[i, 2] < 0) == True or (radius[i, 1] - radius[i, 2] < 0) == True:
            data_r_extra[i,0:3]=radius[i,0:3]
        i+=1
    radius_extra_nonzeros=data_r_extra[np.where(data_r_extra[:,2]!=0)]
    np.savetxt(os.path.join(path_out, 'radius_extra_nonzeros.txt'),
               radius_extra_nonzeros)

    radius_extra_nonzeros=np.loadtxt(os.path.join(path_out, 'radius_extra_nonzeros.txt'))
    line_radius_extra_nonzeros = sum(1 for line in open(os.path.join(path_out, 'radius_extra_nonzeros.txt')))
    radius_extra_nonzeros_left=np.zeros((line_radius_extra_nonzeros,3))
    radius_extra_nonzeros_right = np.zeros((line_radius_extra_nonzeros, 3))
    radius_extra_nonzeros_top = np.zeros((line_radius_extra_nonzeros, 3))
    radius_extra_nonzeros_bottom = np.zeros((line_radius_extra_nonzeros, 3))
    # find where the part of the circle is missing, for example, if the x mius the r smaller than zero then the circle must be at left side of the domain and missing one part
    for i in range (0,line_radius_extra_nonzeros):
        if radius_extra_nonzeros[i,1]-radius_extra_nonzeros[i,2]<0:# if x
            radius_extra_nonzeros_left[i,0:3]=radius_extra_nonzeros[i,0:3]
        elif radius_extra_nonzeros[i,1]+radius_extra_nonzeros[i,2]>10:
            radius_extra_nonzeros_right[i, 0:3] = radius_extra_nonzeros[i, 0:3]
        elif radius_extra_nonzeros[i, 0] + radius_extra_nonzeros[i, 2] >5:#if y
            radius_extra_nonzeros_top[i, 0:3] = radius_extra_nonzeros[i, 0:3]
        elif radius_extra_nonzeros[i, 0] - radius_extra_nonzeros[i, 2] <0:
            radius_extra_nonzeros_bottom[i, 0:3] = radius_extra_nonzeros[i, 0:3]
            i+=1

    radius_extra_nonzeros_left=radius_extra_nonzeros_left[np.where(radius_extra_nonzeros_left[:,2]!=0)]
    radius_extra_nonzeros_right=radius_extra_nonzeros_right[np.where(radius_extra_nonzeros_right[:,2]!=0)]
    radius_extra_nonzeros_top=radius_extra_nonzeros_top[np.where(radius_extra_nonzeros_top[:,2]!=0)]
    radius_extra_nonzeros_bottom=radius_extra_nonzeros_bottom[np.where(radius_extra_nonzeros_bottom[:,2]!=0)]

    np.savetxt(os.path.join(path_out, 'radius_extra_nonzeros_left.txt'), radius_extra_nonzeros_left)
    np.savetxt(os.path.join(path_out, 'radius_extra_nonzeros_right.txt'), radius_extra_nonzeros_right)
    np.savetxt(os.path.join(path_out, 'radius_extra_nonzeros_top.txt'), radius_extra_nonzeros_top)
    np.savetxt(os.path.join(path_out, 'radius_extra_nonzeros_bottom.txt'), radius_extra_nonzeros_bottom)

    line_radius_extra_left = sum(
        1 for line in open(os.path.join(path_out, 'radius_extra_nonzeros_left.txt')))
    line_radius_extra_right = sum(
        1 for line in open(os.path.join(path_out, 'radius_extra_nonzeros_right.txt')))
    line_radius_extra_top = sum(
        1 for line in open(os.path.join(path_out, 'radius_extra_nonzeros_top.txt')))
    line_radius_extra_bottom = sum(
        1 for line in open(os.path.join(path_out, 'radius_extra_nonzeros_bottom.txt')))

    #checking whether the total sum of the four directions are the same numbers lines in file of _extra.txt

    if line_radius_extra_left+line_radius_extra_right+line_radius_extra_top+line_radius_extra_bottom!=line_radius_extra_nonzeros:
        try:
            print("Something is wrong")
        except ImportError:
            raise SystemExit
    elif line_radius_extra_left+line_radius_extra_right+line_radius_extra_top+line_radius_extra_bottom==line_radius_extra_nonzeros:
        print("Great!")

def plot_test():
    radius=np.loadtxt(os.path.join(path_out,'circle_radius.txt'))
    radius_extra_nonzeros_left=np.loadtxt(os.path.join(path_out, 'radius_extra_nonzeros_left.txt'))
    radius_extra_nonzeros_right=np.loadtxt(os.path.join(path_out, 'radius_extra_nonzeros_right.txt'))
    radius_extra_nonzeros_top=np.loadtxt(os.path.join(path_out, 'radius_extra_nonzeros_top.txt'))
    radius_extra_nonzeros_bottom=np.loadtxt(os.path.join(path_out, 'radius_extra_nonzeros_bottom.txt'))

    plt.close('all')
    matplotlib.rc('xtick', labelsize=9)
    matplotlib.rc('ytick', labelsize=9)
    fig=plt.figure(figsize=(6,6))


    xyear1 = radius[:,0]
    yyear1 = radius[:,1]
    ryear1 = radius[:,2]

    line_radius_extra_left = sum(
        1 for line in open(os.path.join(path_out, 'radius_extra_nonzeros_left.txt')))
    line_radius_extra_right = sum(
        1 for line in open(os.path.join(path_out, 'radius_extra_nonzeros_right.txt')))
    line_radius_extra_top = sum(
        1 for line in open(os.path.join(path_out, 'radius_extra_nonzeros_top.txt')))
    line_radius_extra_bottom = sum(
        1 for line in open(os.path.join(path_out, 'radius_extra_nonzeros_bottom.txt')))
    if line_radius_extra_left>1:
    # so if the yzoi (which is in x direction when plotting)it is on the left side
        xextra_nonzeros_left=radius_extra_nonzeros_left[:,0]
        yextra_nonzeros_left =radius_extra_nonzeros_left[:, 1]+10
        rextra_nonzeros_left=radius_extra_nonzeros_left[:,2]
    elif line_radius_extra_left==1:
        xextra_nonzeros_left=radius_extra_nonzeros_left[0]
        yextra_nonzeros_left = radius_extra_nonzeros_left[1]+10
        rextra_nonzeros_left=radius_extra_nonzeros_left[2]
    elif line_radius_extra_left==0:
        xextra_nonzeros_left=0
        yextra_nonzeros_left =0
        rextra_nonzeros_left=0

    if line_radius_extra_right>1:
        # if the yzoi (which is in x direction when plotting)it is on the right side
        xextra_nonzeros_right=radius_extra_nonzeros_right[:,0]
        yextra_nonzeros_right=radius_extra_nonzeros_right[:,1]-10
        zextra_nonzeros_right=radius_extra_nonzeros_right[:,2]
    elif line_radius_extra_right==1:
        # if the yzoi (which is in x direction when plotting)it is on the right side
        xextra_nonzeros_right=radius_extra_nonzeros_right[0]
        yextra_nonzeros_right=radius_extra_nonzeros_right[1]-10
        zextra_nonzeros_right=radius_extra_nonzeros_right[2]
    elif line_radius_extra_right==0:
        xextra_nonzeros_right=0
        yextra_nonzeros_right =0
        zextra_nonzeros_right=0

    if line_radius_extra_top>1:
        # if the x (which is in y direction when plotting)it is on the top side
        xextra_nonzeros_top=radius_extra_nonzeros_top[:,0]-5
        yextra_nonzeros_top = radius_extra_nonzeros_top[:,1]
        zextra_nonzeros_top = radius_extra_nonzeros_top[:,2]
    elif line_radius_extra_top==1:
        xextra_nonzeros_top=radius_extra_nonzeros_top[0]-5
        yextra_nonzeros_top = radius_extra_nonzeros_top[1]
        zextra_nonzeros_top = radius_extra_nonzeros_top[2]
    elif line_radius_extra_top==0:
        xextra_nonzeros_top=0
        yextra_nonzeros_top =0
        zextra_nonzeros_top=0

    if line_radius_extra_bottom > 1:
    # if the x (which is in y direction when plotting)it is on the bottom side
        xextra_nonzeros_bottom=radius_extra_nonzeros_bottom[:,0]+5
        yextra_nonzeros_bottom = radius_extra_nonzeros_bottom[:,1]
        zextra_nonzeros_bottom = radius_extra_nonzeros_bottom[:,2]
    elif line_radius_extra_bottom==1:
        xextra_nonzeros_bottom=radius_extra_nonzeros_bottom[0]+5
        yextra_nonzeros_bottom = radius_extra_nonzeros_bottom[1]
        zextra_nonzeros_bottom = radius_extra_nonzeros_bottom[2]
    elif line_radius_extra_bottom==0:
        xextra_nonzeros_bottom=0
        yextra_nonzeros_bottom =0
        zextra_nonzeros_bottom=0

    ax1=fig.add_subplot(111)
    ax1.plot(yyear1, xyear1, 'o', color='g')
    circles(yyear1, xyear1, ryear1,c='none',fc='none',ec='black',ls='-')
    circles(yextra_nonzeros_left, xextra_nonzeros_left, rextra_nonzeros_left, c='none', fc='none', ec='black', ls='-')
    circles(yextra_nonzeros_right, xextra_nonzeros_right, zextra_nonzeros_right, c='none', fc='none', ec='black', ls='-')
    circles(yextra_nonzeros_top, xextra_nonzeros_top, zextra_nonzeros_top, c='none', fc='none', ec='black', ls='-')
    circles(yextra_nonzeros_bottom, xextra_nonzeros_bottom, zextra_nonzeros_bottom, c='none', fc='none', ec='black', ls='-')


    ax1.set_title("test data",size=11,y=1.02,x=0.5)
    plt.gca().set_aspect('equal', adjustable='box')

    #
    for sub_plot_axis in plt.gcf().get_axes():
        sub_plot_axis.set_xlim(0,10)
        sub_plot_axis.set_ylim(0,5)
    # position=fig.add_axes([0.1, 0.05, 0.8, 0.02])#xmin, ymin, dx, and dy for the subplot
    fig.text(0.5, 0.95, 'Periodic boundary circle plot', ha='center',size=13)
    fig.subplots_adjust(wspace=0.05)

    fig.tight_layout(rect=(0,0,1,0.95))
    plt.savefig(os.path.join(path_out,"Periodic_boundary_circle_complete.png"))
    plt.show()


