import numpy as np
from matplotlib import pyplot as plt


def HandD(x):
    c = (1-np.cos(x))/2
    return (1/2) - (1/4)*c+ (3/2)*c*np.log(c)

def Dipole(x):
    return np.cos(x)

def Monopole(x):
    return 1 + 0*x

#--------------------------------------------------------

def plotBinnedCrossCor(xi,rho,sig,os_A,orf=HandD,bins=10):
    temp = np.arange(0,len(xi),len(xi)/bins,dtype=np.int16)
    ranges = np.zeros(bins+1)
    ranges[0:bins]=temp
    ranges[bins]=len(xi)
    
    xiAvg = np.zeros(bins)
    rhoAvg = np.zeros(bins)
    sigmaComb = np.zeros(bins)
    
    #Need to sort by pulsar separation
    sortMask = np.argsort(xi)
    
    for i in range(bins):
        #Mask and select range of values to average
        subXi = xi[sortMask]
        subXi = subXi[int(ranges[i]):int(ranges[i+1])]
        subRho = rho[sortMask]
        subRho = subRho[int(ranges[i]):int(ranges[i+1])]
        subSig = sig[sortMask]
        subSig = subSig[int(ranges[i]):int(ranges[i+1])]
        
        #Useful to not type this out every time
        subSigSquare = np.square(subSig)
        
        #Average the separations, no weighting
        xiAvg[i] = np.average(subXi)
        
        #Average the correlated amplitude, with weighting
        rhoAvg[i] = np.sum(subRho/subSigSquare)/np.sum(1/subSigSquare)
        
        #Averaging the uncertanties
        sigmaComb[i] = 1/np.sqrt(np.sum(1/subSigSquare))
    
    #Plot data
    plt.errorbar(xiAvg,rhoAvg,yerr=sigmaComb,fmt='ro',label='Data')
    
    #Model plotting
    xvals = np.arange(0.05,np.pi,.05)
    yvals = os_A*orf(xvals)
        
    #Plot sum
    plt.plot(xvals,yvals,'b-',label='ORF')
    plt.legend()
    plt.grid()
    plt.title(f'Binned correlation data with ORF')
    plt.xlabel('Pulsar separation')
    plt.ylabel('Correlated $A^2$')

