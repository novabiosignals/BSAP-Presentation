# version 1.0

# Includes
# filtfilt
# smooth

# Hz
# filterhum50

import scipy as sp

from scipy import *
from scipy import signal

def lowpass(s,f,order=2,fs=1000.0):
    b,a=signal.butter(order,f/fs)
    return signal.lfilter(b,a,s)

def highpass(s,f,order=2,fs=1000.0):
    b,a=signal.butter(order,f/fs,btype='highpass')
    return signal.lfilter(b,a,s)

def bandpass(s,f1,f2,order=2,fs=1000.0):
    b,a=signal.butter(order,[f1/fs,f2/fs],btype='bandpass')
    return signal.lfilter(b,a,s)

def bandstop(s,f1,f2,order=2,fs=1000.0):
    b,a=signal.butter(order,[f1/fs,f2/fs],btype='bandstop')
    return signal.lfilter(b,a,s)


def point(x,y,color='r',ims=15):
    plot([x],[y],color=color,ms=ims,marker='.')


def plotfft(s,fmax,doplot=True):
    fs=abs(fft(s))
    f=linspace(0,fmax/2,len(s)/2)
    if doplot:
        plot(f[1:len(s)/2],fs[1:len(s)/2])
    return f,fs,



from numpy import vstack, hstack, eye, ones, zeros, linalg, \
newaxis, r_, flipud, convolve, matrix, array
from scipy.signal import lfilter

def lfilter_zi(b,a):
    #compute the zi state from the filter parameters. see [Gust96].

    #Based on:
    # [Gust96] Fredrik Gustafsson, Determining the initial states in forward-backward 
    # filtering, IEEE Transactions on Signal Processing, pp. 988--992, April 1996, 
    # Volume 44, Issue 4

    n=max(len(a),len(b))

    zin = (  eye(n-1) - hstack( (-a[1:n,newaxis],
                                 vstack((eye(n-2), zeros(n-2))))))

    zid=  b[1:n] - a[1:n]*b[0]

    zi_matrix=linalg.inv(zin)*(matrix(zid).transpose())
    zi_return=[]

    #convert the result into a regular array (not a matrix)
    for i in range(len(zi_matrix)):
      zi_return.append(float(zi_matrix[i][0]))

    return array(zi_return)




def filtfilt(b,a,x):
    #For now only accepting 1d arrays
    ntaps=max(len(a),len(b))
    edge=ntaps*3

    if x.ndim != 1:
        raise(ValueError, "Filiflit is only accepting 1 dimension arrays.")

    #x must be bigger than edge
    if x.size < edge:
        raise(ValueError, "Input vector needs to be bigger than 3 * max(len(a),len(b).")


    if len(a)!=len(b):
        b=r_[b,zeros(len(a)-len(b))]


    zi=lfilter_zi(b,a)

    #Grow the signal to have edges for stabilizing 
    #the filter with inverted replicas of the signal
    s=r_[2*x[0]-x[edge:1:-1],x,2*x[-1]-x[-1:-edge:-1]]
    #in the case of one go we only need one of the extrems 
    # both are needed for filtfilt

    (y,zf)=lfilter(b,a,s,-1,zi*s[0])

    (y,zf)=lfilter(b,a,flipud(y),-1,zi*y[-1])

    return flipud(y[edge-1:-edge+1])


def Hz(f,FS=1000.0):
	return f/FS*2.0


def filterhum50(s, filter_order=2):


    hum_notch50=array([0.95,1.05])*Hz(50)
    hum_notch100=hum_notch50*2



    [b,a]=signal.butter(filter_order, hum_notch50, btype="bandstop")
    fdata=signal.lfilter(b, a, s)

    [b,a]=signal.butter(filter_order, hum_notch100, btype="bandstop")
    fdata=signal.lfilter(b, a, fdata)
    return fdata


import numpy

def smooth(x,window_len=10,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """ 
     
    if x.ndim != 1:
        raise(ValueError, "smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise(ValueError, "Input vector needs to be bigger than window size.")
        

    if window_len<3:
        return x
    
    
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise(ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
    

    s=numpy.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')
    
    y=numpy.convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]

    
    
    
from numpy import *
from pylab import *
    
def smooth_demo():    
    
    t=linspace(-4,4,100)
    x=sin(t)
    xn=x+randn(len(t))*0.1
    y=smooth(x)
    
    ws=31
    
    subplot(211)
    plot(ones(ws))
    
    windows=['flat', 'hanning', 'hamming', 'bartlett', 'blackman']
    
    hold(True)
    for w in windows[1:]:
        eval('plot('+w+'(ws) )')
    
    axis([0,30,0,1.1])
    
    legend(windows)
    title("The smoothing windows")
    subplot(212)
    plot(x)
    plot(xn)
    for w in windows:
        plot(smooth(xn,10,w))
    l=['original signal', 'signal with noise']
    l.extend(windows)
    
    legend(l)
    title("Smoothing a noisy signal")
    show()


    
    
# Todo compose a demo and test files

if __name__=='__main__':

    from scipy.signal import butter
    from scipy import sin, arange, pi, randn

    from pylab import plot, legend, show, hold

    t=arange(-1,1,.01)
    x=sin(2*pi*t*.5+2)
    #xn=x + sin(2*pi*t*10)*.1
    xn=x+randn(len(t))*0.05

    [b,a]=butter(3,0.05)

    z=lfilter(b,a,xn)
    y=filtfilt(b,a,xn)



    plot(x,'c')
    hold(True)
    plot(xn,'k')
    plot(z,'r')
    plot(y,'g')

    legend(('original','noisy signal','lfilter - butter 3 order','filtfilt - butter 3 order'))
    hold(False)
    show()




def peakDetectPos(x,tol):

    b=+x
    b[b<tol]=0
    return find((diff(numpy.sign(diff(b))))==-2)+1



def priorPeak(V,Vprior):
    return numpy.array([Vprior[find(Vprior<i)][-1] for i in V])
    
    
def postPeak(V,Vpost):
    return numpy.array([Vpost[find(Vpost>i)][-1] for i in V])
