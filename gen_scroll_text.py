#!/usr/bin/env python3
# ##########################################################
# This is a txt version of the 'scrollgenerator' program 
# Rolling up of a sheet placed at xy, xz or yz plane
# the rolling is made around x, y or z axis following an
# archimedean spiral 
# Version 1.01 - Nov/09/2023
# Ricardo Paupitz - ricardo.paupitz@unesp.br
# ---------------------------------------------------------- 
# ##########################################################


import sys 
import readline
import numpy as np
import pprint
import os

#########################################################
#    +++ Numerical functions  +++
#########################################################
# calculates the length of the spiral line
# for the given angle 'angle'
# Archimedean spiral: r=b*angle
def spirallength(angle,b):
    r = np.sqrt((angle**2.0)+1.0)
    lnn = np.log(angle + r)
    length = b*0.5*(angle*r + lnn)
    return(length)


###############################################
# when given the length 'length' for the line 
# on the spiral, this function returns the
# corresponding angle
# Archimedean spiral: r=b*angle
def findangle(length,b): 
    x=0
    n=0
    l=0
    dx=0.01
    while (l <= length):
        x = x+dx
        l = spirallength(x,b)
    ######################################
    tol = 0.001
    while (np.abs(l-length) > tol):
        if (l>length):
            dx = 0.5*dx
            x = x-dx
            l = spirallength(x,b)
            dl = np.abs(l-length)
        elif (l<length):
            dx = 0.5*dx
            x = x+dx
            l = spirallength(x,b)
            dl = np.abs(l-length)
    return(x)
    ######################################       
###############################################    

#######################################
# counts the number of lines of a file
def linecount(filename):
    liness = 0
    for lines in open(filename):
        liness += 1
    return liness
    filename.close()


###########################################################
# verify if there are negative values on x,y,z coordinates
# or if the sheet is too close to the origin
# and dislocates it as necessary
def verify(pl,f):
    disl = 5.0 # minimal distance from x=0,y=0 and z=0
    linn = open(f,'r')
    nl = linecount(f)
    listaux = linn.readlines()
    auxname = "%s_aux.xyz" % f 
    out = open(auxname,'w')
    out.write(listaux[0][:])
    out.write(listaux[1][:])
    linn.close()
    out.close()

    elem = []
    posx = []
    posy = []
    posz = []
    numx = 0
    numy = 0
    numz = 0
    coordinates = open(f,'r')
    for line in range(2,nl):
        d = listaux[line][:].split()
        elem.append(d[0])
        posx.append(float(d[1]))
        posy.append(float(d[2]))
        posz.append(float(d[3]))
    minx = min(posx)
    miny = min(posy)
    minz = min(posz)
    maxx = max(posx)
    maxy = max(posy)
    maxz = max(posz)
    ywidth = maxy -miny
    nx = len(posx)
    ny = len(posy)
    nz = len(posz)
    if ((minx < 0) and (pl == 'xy')):
        corrx = np.abs(minx)+disl
        while (numx < len(posx)):
            posx[numx] += corrx
            numx += 1
    else:
        pass
    
    if ((miny < 0) and (pl == 'xy')):
        corry = np.abs(miny)+disl
        while (numy < len(posy)):
            posy[numy] += corry
            numy += 1
    else:
        pass
    
    if ((minx < 0) and (pl == 'xz')):
        corrx = np.abs(minx)+disl
        while (numx < len(posx)):
            posx[numx] += corrx
            numx += 1
    else:
        pass
    
    if ((minz < 0) and (pl == 'xz')):
        corrz = np.abs(minz)+disl
        while (numz < len(posz)):
            posz[numz] += corrz
            numz += 1
    else:
        pass

    if ((miny < 0) and (pl == 'yz')):
        corry = np.abs(miny)+disl
        while (numy < len(posy)):
            posy[numy] += corry
            numy += 1
    else:
        pass
    
    if ((minz < 0) and (pl == 'yz')):
        corrz = np.abs(minz)+disl
        while (numz < len(posz)):
            posz[numz] += corrz
            numz += 1
    else:
        pass
    
# copy xyz file to an auxiliary temp file
# this file will be used by the 'map()' function
    if ((minx < 0) or (miny < 0) or (minz < 0)):
        auxname = "%s_aux.xyz" % f 
        aux = open(auxname,'a')
        ll=0
        while (ll < len(posx)):
            c = "%s  %s  %s  %s  \n" % (elem[ll],posx[ll],posy[ll],posz[ll])
            aux.write(c)
            ll += 1
    else:
        nlines = linecount(f)
        with open(f) as f1:
            auxname = "%s_aux.xyz" % f 
            with open(auxname, 'w') as f2:
                Lines = f1.readlines()
                for iline in range(nlines):
                    f2.write(Lines[iline])
            f2.close()
        f1.close()
#    clean = "rm -f %s_aux.xyz" % f
#    os.system(clean)
###########################################################
# map nanoribbon positions to a spiral configuration
def map(f,plane,axis,b,positions):
    a = []
    aa = []
    clean = "rm -f temp.txt"
    os.system(clean)
    linn = open(f,'r')
    lnn = linecount(f)
    listaux = linn.readlines()
    a.append(listaux[0][:])
    a.append(listaux[1][:])
    linn.close()
    i = 0
    lnumber = linecount(f)-2
    lin = open(f,'r')
    if (plane == 'xy'):
        if (axis == 'x'):
            for line in range(2,lnn):
                d = listaux[line][:].split()
                i += 1
                theta = findangle(float(d[2]),b)
                ###################################################
                # calculating the cartesian components of the unit
                # vector (versor) tangent to the curve pointing
                # to the increasing theta direction
                mod = float(b*np.sqrt(1.0 + (theta*theta)))
                yt = float((b*np.cos(theta) - b*theta*np.sin(theta))/mod)
                zt = float((b*np.sin(theta) + b*theta*np.cos(theta))/mod)
                # now, the perpendicular versor at the same point
                yp =  1.0*zt*float(d[3])  
                zp = -1.0*yt*float(d[3])  
                # now the new position of the atom, corrected  
                y = b*theta*np.cos(theta) + yp
                z = b*theta*np.sin(theta) + zp
                ###############################################################
                c = "%s  %s  %s  %s  \n" % (d[0],d[1],y,z)
                a.append(c)
        elif (axis == 'y'):
            for line in range(2,lnn):
                d = listaux[line][:].split()
                i += 1
                theta = findangle(float(d[1]),b)
                ###################################################
                # calculating the cartesian components of the unit
                # vector (versor) tangent to the curve pointing
                # to the increasing theta direction
                mod = float(b*np.sqrt(1.0 + (theta*theta)))
                xt = float((b*np.cos(theta) - b*theta*np.sin(theta))/mod)
                zt = float((b*np.sin(theta) + b*theta*np.cos(theta))/mod)
                # now, the perpendicular versor at the same point
                xp =  1.0*zt*float(d[3])  
                zp = -1.0*xt*float(d[3])  
                # now the new position of the atom, corrected  
                x = b*theta*np.cos(theta) + xp
                z = b*theta*np.sin(theta) + zp
                ###############################################################
                c = "%s  %s  %s  %s  \n" % (d[0],x,d[2],z)
                a.append(c)
    elif (plane == 'xz'):
        if (axis == 'x'):
            for line in range(2,lnn):
                d = listaux[line][:].split()
                i += 1
                theta = findangle(float(d[3]),b)
                ###################################################
                # calculating the cartesian components of the unit
                # vector (versor) tangent to the curve pointing
                # to the increasing theta direction
                mod = float(b*np.sqrt(1.0 + (theta*theta)))
                yt = float((b*np.cos(theta) - b*theta*np.sin(theta))/mod)
                zt = float((b*np.sin(theta) + b*theta*np.cos(theta))/mod)
                # now, the perpendicular versor at the same point
                yp =  1.0*zt*float(d[2])  
                zp = -1.0*yt*float(d[2])  
                # now the new position of the atom, corrected  
                y = b*theta*np.cos(theta) + yp
                z = b*theta*np.sin(theta) + zp
                c = "%s  %s  %s  %s  \n" % (d[0],d[1],y,z)
                a.append(c)
        elif (axis == 'z'):
            for line in range(2,lnn):
                d = listaux[line][:].split()
                i += 1       
                theta = findangle(float(d[1]),b)
                ###################################################
                # calculating the cartesian components of the unit
                # vector (versor) tangent to the curve pointing
                # to the increasing theta direction
                mod = float(b*np.sqrt(1.0 + (theta*theta)))
                xt = float((b*np.cos(theta) - b*theta*np.sin(theta))/mod)
                yt = float((b*np.sin(theta) + b*theta*np.cos(theta))/mod)
                # now, the perpendicular versor at the same point
                xp =  1.0*yt*float(d[2])
                yp = -1.0*xt*float(d[2])
                #print('mod_perpend: ',float(np.sqrt(xp*xp+yp*yp)))
                # now the new position of the atom, corrected  
                x = b*theta*np.cos(theta) + xp
                y = b*theta*np.sin(theta) + yp
                c = "%s  %s  %s  %s  \n" % (d[0],x,y,d[3])
                a.append(c)
    elif (plane == 'yz'):
        if (axis == 'y'):
            for line in range(2,lnn):
                d = listaux[line][:].split()
                i += 1
                theta = findangle(float(d[3]),b)
                ###################################################
                # calculating the cartesian components of the unit
                # vector (versor) tangent to the curve pointing
                # to the increasing theta direction
                mod = float(b*np.sqrt(1.0 + (theta*theta)))
                xt = float((b*np.cos(theta) - b*theta*np.sin(theta))/mod)
                zt = float((b*np.sin(theta) + b*theta*np.cos(theta))/mod)
                # now, the perpendicular versor at the same point
                xp =  1.0*zt*float(d[1])
                zp = -1.0*xt*float(d[1])
                #print('mod_perpend: ',float(np.sqrt(xp*xp+yp*yp)))
                # now the new position of the atom, corrected  
                x = b*theta*np.cos(theta) + xp
                z = b*theta*np.sin(theta) + zp                
                c = "%s  %s  %s  %s  \n" % (d[0],x,d[2],z)
                a.append(c)
        elif (axis == 'z'):
            for line in range(2,lnn):
                d = listaux[line][:].split()
                i += 1
                theta = findangle(float(d[2]),b)
                ###################################################
                # calculating the cartesian components of the unit
                # vector (versor) tangent to the curve pointing
                # to the increasing theta direction
                mod = float(b*np.sqrt(1.0 + (theta*theta)))
                xt = float((b*np.cos(theta) - b*theta*np.sin(theta))/mod)
                yt = float((b*np.sin(theta) + b*theta*np.cos(theta))/mod)
                # now, the perpendicular versor at the same point
                xp =  1.0*yt*float(d[1])
                yp = -1.0*xt*float(d[1])
                #print('mod_perpend: ',float(np.sqrt(xp*xp+yp*yp)))
                # now the new position of the atom, corrected  
                x = b*theta*np.cos(theta) + xp
                y = b*theta*np.sin(theta) + yp                
                c = "%s  %s  %s  %s  \n" % (d[0],x,y,d[3])
                a.append(c)
    else:
        print("unknown option")
    lin.close()
    ll = open(positions,'w')
    for text in a:
        ll.write(text)
    ll.close()
    clean = "rm -f %s" % f
    os.system(clean)

# +++ End of the numerical functions' definitions  +++
#########################################################
print('######################################################')
print('This is a txt version of the "scrollgenerator" program')
print('It can map a flat sheet around x, y or z axis')
print('note that:')
print('xy-plane -> x or y rolling up')
print('xz-plane -> x or z rolling up')
print('yz-plane -> y or z rolling up')
print('------------------------------------------------------')
print('    ')
sheet = input('Please, inform the .xyz input file:')
plane = input('plane with which the sheet is aligned (xy, xz or yz):')
axis  = input('axis of rotation:')
output_file = input('Please, inform the output .xyz file:')
distance = float(5.5)
separation = float(3.5)
sep = float(separation)/(2*np.pi)
name_aux = verify(plane,sheet)
input_aux = "%s_aux.xyz" % sheet
map(input_aux,plane,axis,sep,output_file)
print('Scroll is ready!')

 
