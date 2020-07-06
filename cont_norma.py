from numpy import pi, sin
from datetime import date
import glob
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
#from matplotlib import rc
#rc('text', usetex=True)


def mask(a,b,excl):
    nn = len(a)-excl
    xval = np.zeros(nn)
    yval = np.zeros(nn)

    pos = np.argsort(b)
    a1 = a[pos]
    b1 = b[pos]

    xval = a1[excl:]
    yval = b1[excl:]

    pos2 = np.argsort(xval)

    return (xval[pos2],yval[pos2])

def bines(x,y,n,m):
    global xbin,ybin, bin_x, bin_y
    deltax = (max(x) - min(x))/n
    bin_x = np.zeros(n)
    bin_y = np.zeros(n)

    for i in range(n):
        aa = np.where((x >= (min(x) + (deltax * i))) & ( x <= (min(x) + (deltax * (i+1)))))[0]
        bin_x[i] = min(x[aa]) + (max(x[aa]) - min(x[aa]))/2.
        bin_y[i] = clipping(y[aa])
    xbin,ybin = mask(bin_x,bin_y,m)

    
def clipping(yc,sig=0.5):
    l0 = 1 
    l1 = 0
    yold = yc

    
    med = np.nanmedian(yold)
    std = np.nanstd(yold)
    pos = np.where((yold >= (med - (sig * std))) & ( yold <= (med + (sig *std))))[0]
    ynew = yold[pos]
    return np.median(ynew)
        
def fit(nord,nbins,exbin):
    global x,y,xbin,ybin
    bines(x,y,nbins,exbin)
    a = np.polyfit(xbin,ybin, nord)
    p = np.poly1d(a)
    
    return(p(x))

#------- order functions
def order_plus(mouse_event):
    global x,y,nord,nbins, bin_x, bin_y
    nord = nord + 1
    if nbins <= nord:
        nbins = nord + 1

    update()
    t1.set_text('%d' % nord)
    

def order_init(mouse_event):
    global x,y,nord,nbins,nord_0
    
    nord = nord_0

    update()
    t1.set_text('%d' % nord)
   

def order_minus(mouse_event):
    global x,y,nord,nbins
   
    if nord > 1:
        nord = nord - 1
    else:
        nord = 1

    update()
    t1.set_text('%d' % nord)
    
#------- bin functions
def bin_plus(mouse_event):
    global x,y,nord, nbins
    nbins = nbins + 1

    update()
    t2.set_text('%d' % nbins)

def bin_init(mouse_event):
    global x,y,nord,nbins, nbins_0
    nbins = nbins_0

    update()
    t2.set_text('%d' % nbins)

def bin_minus(mouse_event):
    global x,y,nord, nbins
   
    if nbins > 2:
        nbins = nbins - 1
    else:
        nbins = 2

    update()
    t2.set_text('%d' % nbins)

#------- excluded bin functions
def exbin_plus(mouse_event):
    global x,y,nord, nbins, exbin
    if exbin < nbins -2:
        exbin = exbin + 1
    else:
        exbin = nbins - 2

    update()
    t3.set_text('%d' % exbin)

def exbin_init(mouse_event):
    global x,y,nord,nbins, exbin
    exbin = 0
    
    update()
    t3.set_text('%d' % exbin)

def exbin_minus(mouse_event):
    global x,y,nord, nbins, exbin
   
    if exbin > 0:
        exbin = exbin - 1
    else:
        exbin = 0

    update()
    t3.set_text('%d' % exbin)
    

#------- y-offset functions
def yoff_plus(mouse_event):
    global x,y,nord, nbins, exbin,yof
    if yof < 0.5:
        yof = yof + valoff
    else:
        yof = 0.5

    update()
    t4.set_text('%.3f' % yof)

def yoff_init(mouse_event):
    global x,y,nord,nbins, exbin, yof, yof_0
    yof = yof_0

    update()
    t4.set_text('%.2f' % yof)

def yoff_minus(mouse_event):
    global x,y,nord, nbins, exbin, yof
   
    if yof > -0.5:
        yof = yof - valoff
    else:
        yof = -0.5
    
    update()
    t4.set_text('%.3f' % yof)
    
#---- Update funtions
def update():
    global x,y,nord,nbins,exbin,yof,xbin,ybin, bin_x, bin_y, ax2
    
    ajus.set_ydata(fit(nord,nbins,exbin))
    norm.set_ydata(y/fit(nord,nbins,exbin) + yof)
    bins.set_data(xbin,ybin)
    bins2.set_data(bin_x,bin_y)
    ax2.set_ylim(min(y/fit(nord,nbins,exbin) + yof),1.1)
    fig.canvas.draw_idle()

def update2():
    global x,y,nord,nbins,exbin,yof,xbin,ybin, bin_x, bin_y, fig, nord_0,nbins_0,ajus,norm,bins2,bins

    fig = plt.figure(figsize=(13,6.5))
    ax = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    
    fig.subplots_adjust(hspace=0., bottom=0.25)
   

#-------select region function
def sel_reg(select):
    global xx,yy,x,y,zz,z, regions

    lim = regions[select-1]
    aa = np.where((xx >= lim[0]) & (xx <= lim[1]))[0]

    x = xx[aa]
    y = yy[aa]
    z = zz[aa]
    
    
    nreg = len(regions)
    plt.close()
    initial()
   
   
    
#---- Spectral region navigation functions   
def r_prev(mouse_event):
    global regions,reg,x,y, sal
    if sal == 0:
        save(True)
        
    if reg > 1:
        reg = reg - 1
        
    else:
        reg = 1

    sal = 0
    sel_reg(reg)
    
def r_next(mouse_event):
    global regions,reg,x,y, sal
    if sal == 0:
        save(True)

    nreg = len(regions)
    if reg < nreg:
        reg = reg + 1
        sal = 0
        sel_reg(reg)
    else:
        reg = len(regions)
        sal = 0
        sel_reg(reg)
       

#--- Initialization funtion
def initial():
    global nord,nord_0,nbins_0,nbins, yof, yof_0, sal, exbin

    nbins = nbins_0
    nord = nord_0
    yof = yof_0
    sal = 0
    exbin = 0

#-- Save region function
def save(mouse_event):
    global x,y,z, sal,nord, nbins, exbin, yof, name, reg, sal
    print('Region saved')
    print('-----')
    flux_norm = y/fit(nord,nbins,exbin) + yof
    ff = open(name+'-'+str(reg)+'.txt', 'w')
    for i in range(len(x)):
        if ~np.isnan(flux_norm[i]): 
            ff.write('%f,%f,%d\n' % (x[i], flux_norm[i], z[i]))
    ff.close()
    log(nord,nbins,exbin,yof,reg)
    sal = 1

#-- exit function
def salir(mouse_event):
    exit()

#-- next object function
def nobj(mouse_event):
    global reg, sal

    if sal == 0:
        save(True)
        
    val = glob.glob(name+'*.txt')
    files = ''
    for item in val:
        files = files + ' ' + item
    
    comand = 'cat ' + files + ' > ' + name+'_reg.dat'
    os.system(comand)

    comand = 'mv ' + name+'_reg.dat final/'
    os.system(comand)

    comand = 'mv ' + files + ' regions/'
    os.system(comand)

    plt.close()
    reg = 99

#--- normalization log function
def log(xo,bo,ebo,yo,pos):
    global name
    today = date.today()
    gg = open('log.txt', 'a')

    gg.write('%s,%d,%d,%d,%.3f,%d,%s\n' % (name,xo,bo,ebo,yo,pos,today))
    gg.close()
    
##---- main routine to normalize spectra
def norma(lista, direc, orden = 1, bins_ini = 20, ebins = 0, yof_ini = 0, valoff_ini = 0.005, reg_ini = 1, regs=[[15600,15650], [16720,16800], [22010,22130], [22170,22280],[22530,22730], [22920,23040]]):
    global nord, x, y, z, ax, nbins, xbin,ybin,xx, yy,zz, name,ajus,norm,bins,bins2,nbins_0,nord_0,yof_0,exbin,fig,t2,t3,t4,t1,reg,valoff, regions, ax2

    
    data = pd.read_csv(lista)
    for jj in range(len(data.obj)):

    ##initial values
        nord_0 = orden
        nbins_0 = bins_ini
        exbin = ebins
        yof_0 = yof_ini
        valoff = valoff_ini
        reg = reg_ini
        regions= regs
        name = data.obj[jj] 
        sal = 0
#--------------------

        xx,yy,zz = np.loadtxt(direc+ data.file[jj], unpack=True)   #Loading data



        while reg != 99:
            aa = np.where((xx >= regions[reg-1][0]) & (xx <= regions[reg-1][1]))[0]

            # constraining the data to the specific spectral region
            x = xx[aa]
            y = yy[aa]
            z = zz[aa]

            fig = plt.figure(figsize=(13,6.5))
            ax = fig.add_subplot(211)
            ax2 = fig.add_subplot(212)


            
            #Creating figure and drwing all the buttons
            fig.subplots_adjust(hspace=0., bottom=0.25)

            [ajus] = ax.plot(x, fit(nord_0,nbins_0,exbin), linewidth=3, color='k', ls='--')
            [norm] = ax2.step(x,y/fit(nord_0,nbins_0,exbin) + yof_0, color='g', where='mid', lw=2.)

            [bins2] = ax.plot(bin_x,bin_y,'rx',ms=9, mew=2.)
            [bins] = ax.plot(xbin,ybin,'ks',ms=12 ,mfc='w', mew=2.)
            ax.step(x, y, color='k', where='mid', lw=2., alpha=0.5)

            initial()

            # Order buttons
            ordplus_ax = fig.add_axes([0.28, 0.13, 0.05, 0.04])
            ordplus_button = Button(ordplus_ax, '+', color='#A2D9CE', hovercolor='0.975')
            ordplus_button.on_clicked(order_plus)

            ordinit_ax = fig.add_axes([0.28, 0.08, 0.05, 0.04])
            ordinit_button = Button(ordinit_ax, 'Reset', color='#A2D9CE', hovercolor='0.975')
            ordinit_button.on_clicked(order_init)

            ordmin_ax = fig.add_axes([0.28, 0.03, 0.05, 0.04])
            ordmin_button = Button(ordmin_ax, '-', color='#A2D9CE', hovercolor='0.975')
            ordmin_button.on_clicked(order_minus)

            fig.text(0.23,0.09, 'Order', size=15)
            t1 = fig.text(0.25, 0.055, '%d' % nord_0, size=15, color='r', fontweight='bold')
            
            # bins buttons
            binplus_ax = fig.add_axes([0.40, 0.13, 0.05, 0.04])
            binplus_button = Button(binplus_ax, '+', color='#A2D9CE', hovercolor='0.975')
            binplus_button.on_clicked(bin_plus)

            bininit_ax = fig.add_axes([0.40, 0.08, 0.05, 0.04])
            bininit_button = Button(bininit_ax, 'Reset', color='#A2D9CE', hovercolor='0.975')
            bininit_button.on_clicked(bin_init)

            binmin_ax = fig.add_axes([0.40, 0.03, 0.05, 0.04])
            binmin_button = Button(binmin_ax, '-', color='#A2D9CE', hovercolor='0.975')
            binmin_button.on_clicked(bin_minus)

            fig.text(0.36,0.09, 'Bins', size=15)
            t2 = fig.text(0.36, 0.055, '%d' % nbins_0, size=15, color='r', fontweight='bold')
    
            # excluded bins buttons
            exbinplus_ax = fig.add_axes([0.53, 0.13, 0.05, 0.04])
            exbinplus_button = Button(exbinplus_ax, '+', color='#A2D9CE', hovercolor='0.975')
            exbinplus_button.on_clicked(exbin_plus)

            exbininit_ax = fig.add_axes([0.53, 0.08, 0.05, 0.04])
            exbininit_button = Button(exbininit_ax, 'Reset', color='#A2D9CE', hovercolor='0.975')
            exbininit_button.on_clicked(exbin_init)

            exbinmin_ax = fig.add_axes([0.53, 0.03, 0.05, 0.04])
            exbinmin_button = Button(exbinmin_ax, '-', color='#A2D9CE', hovercolor='0.975')
            exbinmin_button.on_clicked(exbin_minus)

            fig.text(0.455,0.09, 'excl. bins', size=14)
            t3 = fig.text(0.50, 0.055,  '%d' % exbin, size=15, color='r', fontweight='bold')

            # yoffset buttons
            yplus_ax = fig.add_axes([0.65, 0.13, 0.05, 0.04])
            yplus_button = Button(yplus_ax, '+', color='#A2D9CE', hovercolor='0.975')
            yplus_button.on_clicked(yoff_plus)

            yinit_ax = fig.add_axes([0.65, 0.08, 0.05, 0.04])
            yinit_button = Button(yinit_ax, 'Reset', color='#A2D9CE', hovercolor='0.975')
            yinit_button.on_clicked(yoff_init)

            ymin_ax = fig.add_axes([0.65, 0.03, 0.05, 0.04])
            ymin_button = Button(ymin_ax, '-', color='#A2D9CE', hovercolor='0.975')
            ymin_button.on_clicked(yoff_minus)

            fig.text(0.585,0.09, 'y-offset', size=15)
            t4 = fig.text(0.59, 0.055, '%.3f' % yof_0, size=15, color='r', fontweight='bold')


            # previous region button
            prev_ax = fig.add_axes([0.75, 0.13, 0.05, 0.04])
            prev_button = Button(prev_ax, 'Prev. reg', color='#A2D9CE', hovercolor='0.975')
            prev_button.on_clicked(r_prev)

            # next region button
            next_ax = fig.add_axes([0.85, 0.13, 0.05, 0.04])
            next_button = Button(next_ax, 'Next reg', color='#A2D9CE', hovercolor='0.975')
            next_button.on_clicked(r_next)

            # save normalized region button
            save_ax = fig.add_axes([0.86, 0.03, 0.04, 0.04])
            save_button = Button(save_ax, 'Save', color='#A2D9CE', hovercolor='0.975')
            save_button.on_clicked(save)

            # exit button
            sal_ax = fig.add_axes([0.75, 0.03, 0.04, 0.04])
            sal_button = Button(sal_ax, 'Exit', color='#A2D9CE', hovercolor='0.975')
            sal_button.on_clicked(salir)

            # next object button
            nobj_ax = fig.add_axes([0.80, 0.03, 0.05, 0.04])
            nobj_button = Button(nobj_ax, 'Next obj', color='#A2D9CE', hovercolor='0.975')
            nobj_button.on_clicked(nobj)
        



            ax2.set_xlabel('wavelength (angstrom)', size=13)
            ax.set_title(name, fontsize=30) 
            ax.grid(True)
            ax2.grid(True)
            ax2.set_ylim(min(y/fit(nord,nbins,exbin) + yof),1.1)
            plt.show()
    exit()
    
