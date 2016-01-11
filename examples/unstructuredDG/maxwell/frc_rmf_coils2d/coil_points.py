# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 10:31:34 2016

@author: sousae
"""

from math import *

angle1=(pi/2.+2.*pi*(.0127/(2.*pi*.03))/2.)
angle2=(pi/2.-2.*pi*(.0127/(2.*pi*.03))/2.)

radius = 0.03

for k in range(8):
    x1 = radius*cos(angle1+k*pi/4.)
    y1 = radius*sin(angle1+k*pi/4.)
    
    x2 = radius*cos(angle2+k*pi/4.)
    y2 = radius*sin(angle2+k*pi/4.)
    print "========\n"
    print "== coil %1d ==\n" %(k+1)
    print "x1 =%5.5f\n " %x1
    print "y1 =%5.5f\n " %y1
    print "x2 =%5.5f\n " %x2
    print "y2 =%5.5f\n " %y2
    print "========\n"