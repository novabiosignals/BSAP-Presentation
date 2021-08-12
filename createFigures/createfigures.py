#!/usr/bin/python
# -*- coding: latin-1 -*-


import matplotlib.font_manager
from pylab import *

from tools import *
from pandas import read_csv
from scipy.io import loadmat





def savefigs(name):
    import os
    savefig('../figures/new_figures/'+name+'.pdf', dpi=150)
    # os.system('epspdf figures/'+name+'.eps')


##########
#Initial configurations

def pylabconfig():

    rc('lines',linewidth=2,color='k')
    
    # rc('font',**{'family':'serif','serif':['Palatino']})
    rc('font',style='italic', size=10)
    
    
    rc('text', color='grey')
    
    #rc('text', usetex=True)
    
    rc('figure', figsize= (8, 5), dpi=80) 
    rc('axes', grid= False, edgecolor= '#ef4e23ff', labelsize=10,)
    # rc('grid', color= 'grey')
    rc('xtick',color= '#ef4e23ff', labelsize=10)
    rc('ytick',color= '#ef4e23ff', labelsize=10)
    
    close('all')



def plotwithhist(t,s,bins=50):

    from matplotlib.ticker import NullFormatter
    
    nullfmt=NullFormatter()
    # figure()
    ax1=axes([0.125,0.1,0.5,0.8])
    ax2=axes([0.125+0.5,0.1,0.2,0.8])

    tft_plot(t,s, [min(t), max(t), min(s), max(s)], "hist", "sighist", ax=ax1, save=False)
    ax1.set_xticks(ax1.get_xticks()[:-1])

    ax2.hist(s,bins,density=True, facecolor="#ef4e23ff", alpha=0.5, orientation='horizontal',lw=2)
    ax2.axis([0,1,ax1.axis()[2],ax1.axis()[3]])
    ax2.yaxis.set_major_formatter(nullfmt)
    ax2.set_xticks([0,0.5,1])
    savefigs("emgwithhists")


    ###########

def tft_plot(t, s, axis_i, title_i, name, linestyle="-", ax=None, save=True):
    # figure()
    xlabel('time (s)', color="#ef4e23ff")
    ylabel('amplitude', color="#ef4e23ff")
    title(title_i, fontsize=12, color="#ef4e23ff")

    if(ax==None):
        ax = subplot(111)

    ax.plot(t, s, c="k", linestyle=linestyle)
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.axis(axis_i)
    if(save):
        savefigs(name)
        close()

def plotwithstats(t,s):

    from matplotlib.ticker import NullFormatter

    nullfmt=NullFormatter()
    figure()
    ax1=axes([0.125,0.1,0.5,0.8])
    ax2=axes([0.125+0.5,0.1,0.2,0.8])


    tft_plot(t,s, [min(t), max(t), min(s), max(s)], "EMG with Stats", "emg_stats", ax=ax1, save=False)
    ax1.set_xticks(ax1.get_xticks()[:-1])

    meanv=s.mean()
    mi=s.min()
    mx=s.max()
    sd=s.std()

    ax2.bar(-0.5,mx-mi,1,mi,lw=2,color="#ef4e23ff", alpha=0.2)
    ax2.bar(-0.5,sd*2,1,meanv-sd,lw=2,color="#ef4e23ff", alpha=0.5)
    ax2.bar(-0.5,0.2,1,meanv-0.1,lw=2,color="#ef4e23ff")
    ax2.axis([-1,1,ax1.axis()[2],ax1.axis()[3]])

    ax2.yaxis.set_major_formatter(nullfmt)
    ax2.set_xticks([])
    savefigs("emgwithstats")





#########
def syntheticdata():
    figure()
    t = arange(0.0, 8, 0.05)
    l = len (t)
    x = sin(t*2) 
    tft_plot(t, x, axis_i=[min(t), max(t), min(x), max(x)], title_i="Periodic", name="sinewave", save=False)
    axis('equal')

    xlabel('time (s)',color='grey')
    ylabel('amplitude',color='grey')
    title('Periodic', fontsize=12, color='grey')

    savefigs('sinewave')

    ###########

    #non periodic from emg

    # discrete

    figure()
    xlabel('time (s)',color='grey')
    ylabel('amplitude',color='grey')
    title('Discrete signal', fontsize=12, color='grey')

    tft_plot(t[::4],(x*8).astype('i')[::4], [min(t[::4]), max(t[::4]), min(x), max(x)], "Discrete signal", "discretesinewave", save=False)
    axis('equal')

    savefigs('discretesinewave')

def temp():
    #temperature
    s = loadtxt("../data/temp.txt")[:,2]
    s = smooth(s, 100)
    fs = 1000
    vcc = 3
    n = 16
    ntcV = (s*vcc)/(2**n)
    ntcR = (10000 * ntcV) / (vcc-ntcV)
    a0 = 0.00112764514
    a1 =0.000234282709
    a2 = 0.0000000877303013
    tempK = 1/(a0+a1*log(ntcR)+a2*(log(ntcR))**3)
    tempT = tempK - 273.15
    t = linspace(0, len(s)/fs, len(s))

    figure()
    tft_plot(t, tempT, [0, max(t), 23.6, 24.8], "Temperature", "temperature")
    close()

def lux():
    s = loadtxt("../data/lux.txt")[:,2]
    n = 16
    fs = 100
    s = (s/(2**n))*100
    t = linspace(0, len(s) / fs, len(s))
    s = smooth(s, 10)
    # s = s - mean(s)
    # s = smooth(s[:2500], 100)
    tft_plot(t, s, [0, max(t), 50, 100], "Lux", "lux")
    close()

def goniometer():
    s = loadtxt("../data/goniometer.txt")[:, 3]
    n = 16
    fs = 100
    vcc =3
    s = (((vcc/ ((2 ** n)-1)) * s)-vcc/2)/((vcc/2)*606*0.00001)
    s = smooth(s, 10)
    t = linspace(0, len(s) / fs, len(s))

    # s = s - mean(s)
    # s = smooth(s[:2500], 100)

    tft_plot(t, s, [0, max(t), 0, 100], "Goniometer", "goniometer")
    close()


def signalsdone():
################ BVP
    color_i = "#ef4e23ff"
    [t,s]=loadtxt('../data/bvp.txt.gz')

    t=t[:1024]
    s=s-mean(s)
    s=smooth(s[:1024],40)

    tft_plot(t, s, [0, max(t), min(s), max(s)], 'Blood Volume Pressure', "bvp")
    close()
    
    #ecg clean
    
    [t,s]= loadtxt('../data/cleanecg.txt.gz')
    tft_plot(t, s*100, [0, max(t), -25000, 30000], 'ECG', "cleanecg")
    close()

    #ecg noisy
    [t,s]= loadtxt('../data/noisyecg.txt.gz')
    tft_plot(t, s*100, [0,8,-3,3], 'ECG', "noisyecg")
    close()

    #eda noisy
    [t,s]= loadtxt('../data/eda.txt.gz')
    tft_plot(t, s, [0,40,0,1200], 'Actividade electrodermica', "eda")
    close()

    #xyz
    #xyz cal
    [t,x,y,z]= loadtxt('../data/xyzcal.txt.gz')


    ix = arange(len(t))
    ix = ix[::10]

    figure()
    ax = plt.subplot(111)
    ax.set_xlabel('time (s)',color=color_i)
    ax.set_ylabel('amplitude (mV)',color=color_i)
    #    title('Electrodermal activity', fontsize=12, color='grey')
    ax.set_title(u'Acelerómetro', fontsize=12, color=color_i)


    ax.plot(t[ix],smooth(x[ix]), color="k")
    ax.plot(t[ix],smooth(y[ix]),'--',color='k')
    ax.plot(t[ix],smooth(z[ix]),'-.',color='k')
    ax.axis([0,25,400,1200])
    legend(('x','y','z'))
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')

    savefigs('xyzcal')


    #xyz 5 steps and fall

    [t,x,y,z]= loadtxt('../data/xyz.txt.gz')
    figure()
    ax = plt.subplot(111)
    ax.set_xlabel('time (s)',color=color_i)
    ax.set_ylabel(u'aceleração ($N/s^2$)',color=color_i)
    #    title('Electrodermal activity', fontsize=12, color='grey')
    ax.set_title(u'Acelerómetro', fontsize=12, color=color_i)
    
    ix=arange(len(t))
    ix=ix[::1]
    ax.plot(t[ix],smooth(x[ix],3), color="k")
    ax.plot(t[ix],smooth(y[ix],3),'--',color='k')
    ax.plot(t[ix],smooth(z[ix],3),'-.',color='k')
    ax.axis([0,8,-4,4])
    legend(('x','y','z'))
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    savefigs('xyzsteps')


   #respirataion 
    [t,s]= loadtxt('../data/resp.txt.gz')
    tft_plot(t, s, [0,30,700,900], 'Respiração', "resp")
    close()

    #force
    [t,s]= loadtxt('../data/force.txt.gz')
    tft_plot(t[::20], s[::20], [0, max(t[::20]), 0, 120], 'Força', "força")
    close()

    #emg

    [t,s]= loadtxt('../data/emg.txt.gz')
    tft_plot(t[::2], s[::2], [0,max(t[::2]), min(s[::2])-(0.1*(max(s)-min(s))),max(s[::2])+(0.1*(max(s)-min(s)))], u'Electromiografia', "emg")
    close()

    
    #switch
    [t,s]= loadtxt('../data/switch.txt.gz')
    tft_plot(t, s,[0,1,0,5000],u'Interruptor', "switch")
    close()
    
def signalsdoing():
    pass    
    
    
def stats():

    #ecg + hist
    # figure()
    plotwithhist(arange(100),randn(100),bins=20)
    # figure()
    plotwithstats(arange(100),randn(100))


def noise():
    ###########
    figure()
    t = arange(0.0, 8, 0.05)
    l = len (t)
    x=randn(l)*.5+4
    tft_plot(t,x, [0,8,-10,10], "$\mu = 4, \sigma = \\frac{1}{2}$", "noise1")

    #######
    figure()
    t = arange(0.0, 8, 0.05)
    x=randn(l)*2
    tft_plot(t,x, [0,8,-10,10], "$\mu = 0, \sigma = 2$", "noise2")

def processing():
    #derivative
    pass

def applications():
    #Linear regression
    
    #Linear regression example
    # This is a very simple example of using two scipy tools 
    # for linear regression, polyfit and stats.linregress
    
    #Sample data creation
    #number of points 
    n=50
    t=linspace(-5,5,n)
    #parameters
    a=0.8; b=-4
    x=polyval([a,b],t)
    #add some noise
    xn=x+randn(n)
    
    #Linear regressison -polyfit - polyfit can be used other orders polys
    (ar,br)=polyfit(t,xn,1)
    xr=polyval([ar,br],t)
    #compute the mean square error
    err=sqrt(sum((xr-xn)**2)/n)
    
    print('Linear regression using polyfit')
    print('parameters: a=%.2f b=%.2f \nregression: a=%.2f b=%.2f, ms error= %.3f' % (a,b,ar,br,err))
    
    #matplotlib ploting
    title(u'Regressão linear')
    ax = plt.subplot(111)
    tft_plot(t,x,[min(t), max(t), -9, 2], "", "", linestyle="--", ax=ax, save=False)
    tft_plot(t,xn,[min(t), max(t), -9, 2], "", "",linestyle=':', ax=ax, save=False)
    tft_plot(t,xr,[min(t), max(t), -9, 2], "", "", linestyle='-.', ax=ax, save=False)
    legend(['original',u'com ruído', 'regression'])

    
    #interpolation up and down
    pass

if __name__ == '__main__':
    pylabconfig()
    # syntheticdata()

    # signalsdone()
    # signalsdoing()
    # lux()
    #
    # goniometer()
    # temp()
    # stats()
    # noise()
    # processing()
    # applications()

    show()
