import numpy as np
import glob
import matplotlib.pyplot as plt
import scipy.linalg as sl
import pandas as pd

import enterprise
from enterprise.pulsar import Pulsar

import enterprise_extensions
from enterprise_extensions.frequentist import optimal_statistic
from enterprise_extensions import models, model_utils, hypermodel

import corner
from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc
from simFuncs import *
from enterprisePTAs import crnPTA, mcSample
from OS_and_plotting import HandD,Dipole,Monopole,plotBinnedCrossCor

import os, glob
import libstempo as T2
from libstempo import toasim as LT
from libstempo import plot as LP

from simFuncs import *
import argparse


# This script should take the paramaters needed to create single sources and orf paramaters and return the orf distribution

# Run this with 
# run combine_run.py -data /fred/oz002/users/mmiles/VIPER_SummerSchool/mdc2 -orf hd -results /fred/oz002/users/mmiles/VIPER_SummerSchool/results_mass5e9_freq2e-8_dist60_real2 -mass 5e9 -freq 2e-8 -distance 60 -realisations 2

#But with other paramaters

parser = argparse.ArgumentParser(description="Creates a bunch of different ORF results")
parser.add_argument("-data", dest="data", help="data directory to use", required = True)
parser.add_argument("-orf", dest="orf", help="overlap reduction function to use", required = True)
parser.add_argument("-results", dest="results", help="result directory to use", required = True)
parser.add_argument("-mass", dest="mass", help="(Optional) single source mass", required = False)
parser.add_argument("-freq", dest="freq", help="(Optional) single source orbital frequency", required = False)
parser.add_argument("-distance", dest="distance", help="(Optional) single source distance (kpc)", required = False)
parser.add_argument("-theta", dest="theta", help="(Optional) Theta position angle", required = False)
parser.add_argument("-phi", dest="phi", help="(Optional) Phi position angle", required = False)

parser.add_argument("-realisations", dest="realisations", help="(Optional) Adding this flag will let you change the number of realisations. If not it'll be 10.", required=False)

args = parser.parse_args()

orf = str(args.orf)
results_dir = str(args.results)
datadir = str(args.data)

ss_mass = np.float128(args.mass)
ss_freq = np.float128(args.freq)
ss_distance = np.float128(args.distance)

ss_theta = np.float128(args.theta)
ss_phi = np.float128(args.phi)

if args.theta is None:
    ss_theta = 1.75
if args.phi is None:
    ss_phi = 5.

if not os.path.exists(results_dir):
    try:
        os.makedirs(results_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

if args.realisations is not None:
    realisations =int(args.realisations)


parfiles = sorted(glob.glob(datadir + '/*.par'))

# define observation times and duration
Tobs = 20.0 # years
deltaT = 20.0 # points per year
obstimes = 53000.0 + np.arange(0.0, Tobs*365.25, 365.25/20.0)
tref = obstimes[0]*86400.0

# create pulsar objects from parfiles
# default TOA errors are set to 0.4 microseconds


def takeClosest(num,collection):
    snr_select =  min(collection,key=lambda x:abs(x-num))
    return np.where(collection == snr_select)[0]



ii = 0
snrbin = []
xi_bin = []
rho_bin = []
sig_bin = []
os_bin = []
os_sig_bin = []


if realisations is not None:
    while ii < realisations+1:
        
        psrs = create_psrs(parfiles, obstimes=obstimes)


        Epsrs = [Pulsar(p) for p in psrs]

        th = [p.theta for p in Epsrs]
        ph = [p.phi for p in Epsrs]

        # add GWB to our pulsars
        # have the option of changing amplitude and spectral index (gamma)
        add_gwb(psrs)

        if ss_mass is not None:

            # create parameter dictionary for CW signal
            # change params as needed
            pdict = {'gwtheta': ss_theta,
                     'gwphi': ss_phi,
                     'mc': ss_mass,
                     'dist': ss_distance, 
                     'fgw': ss_freq,
                     'phase0': 0.0,
                     'psi': np.pi/4.0,
                     'inc': 0.0}

            # add CW signal to our pulsars
            # prints out the name of each pulsar as it loops through
            # included option to change the number of iterations for the timing model fit
            add_cgw(psrs, pdict, tref, iters=2)

            name = "{}_{}_{}_{}_{}_{}".format(orf, pdict["gwtheta"], pdict["gwphi"], pdict["mc"], pdict["dist"], pdict["fgw"])

            save_sims(psrs, outdir=results_dir+'/')

        else:
            save_sims(psrs, outdir=results_dir+'/sims_no_source')
        
        psrs = lt2ent(psrs)
        
        pta = crnPTA(psrs,fixedGamma=True)
    
        OS = optimal_statistic.OptimalStatistic(psrs,bayesephem=False,gamma_common=4.33,
                                                orf=orf,pta=pta)

        xi,rho,sig,os,os_sig = OS.compute_os(params={'gw_log10_A':-14})
        
        xi_bin.append(xi)
        rho_bin.append(rho)
        sig_bin.append(sig)
        os_bin.append(os)
        os_sig_bin.append(os_sig)

        snr = os/os_sig
        snrbin.append(snr)

        ii = ii+1

    snrbin = np.array(snrbin)
    
    median_index = takeClosest(np.median(snr),snrbin)[0]
    snr = snrbin[median_index]
    
    data = [xi_bin[median_index], rho_bin[median_index], sig_bin[median_index]]
    
    df = pd.DataFrame(np.array(data).T,columns = ["xi","rho","sig"])

else:   
    while ii < 11:
        
        psrs = create_psrs(parfiles, obstimes=obstimes)


        Epsrs = [Pulsar(p) for p in psrs]

        th = [p.theta for p in Epsrs]
        ph = [p.phi for p in Epsrs]

        # add GWB to our pulsars
        # have the option of changing amplitude and spectral index (gamma)
        add_gwb(psrs)

        if ss_mass is not None:

            # create parameter dictionary for CW signal
            # change params as needed
            pdict = {'gwtheta': ss_theta,
                     'gwphi': ss_phi,
                     'mc': ss_mass,
                     'dist': ss_distance, 
                     'fgw': ss_freq,
                     'phase0': 0.0,
                     'psi': np.pi/4.0,
                     'inc': 0.0}

            # add CW signal to our pulsars
            # prints out the name of each pulsar as it loops through
            # included option to change the number of iterations for the timing model fit
            add_cgw(psrs, pdict, tref, iters=2)

            name = "{}_{}_{}_{}_{}_{}".format(orf, pdict["gwtheta"], pdict["gwphi"], pdict["mc"], pdict["dist"], pdict["fgw"])

            save_sims(psrs, outdir=results_dir+'/')

        else:
            save_sims(psrs, outdir=results_dir+'/sims_no_source')
        
        psrs = lt2ent(psrs)
        
        pta = crnPTA(psrs,fixedGamma=True)
    
        OS = optimal_statistic.OptimalStatistic(psrs,bayesephem=False,gamma_common=4.33,
                                                orf=orf,pta=pta)

        xi,rho,sig,os,os_sig = OS.compute_os(params={'gw_log10_A':-14})
        
        xi_bin.append(xi)
        rho_bin.append(rho)
        sig_bin.append(sig)
        os_bin.append(os)
        os_sig_bin.append(os_sig)

        snr = os/os_sig
        snrbin.append(snr)

        ii = ii+1

    snrbin = np.array(snrbin)
    
    median_index = takeClosest(np.median(snr),snrbin)[0]
    snr = snrbin[median_index]
    
    data = [xi_bin[median_index], rho_bin[median_index], sig_bin[median_index]]
    
    df = pd.DataFrame(np.array(data).T,columns = ["xi","rho","sig"])
    

df.to_pickle(results_dir+'/data_'+name)

os = os_bin[median_index]
o_sig = os_sig_bin[median_index]
    
with open(results_dir+'/params_'+name+'.txt', 'w') as f:
    
    f.write("OS_sig: {} \n".format(os_sig))
    f.write("OS: {} \n".format(os))
    f.write("SNR: {} \n".format(snr))

if ss_mass is not None:
    with open(results_dir+'/params_'+name+'.txt', 'w') as f:
        for key, value in pdict.items(): 

            f.write('%s:%s\n' % (key, value))

        f.close()
else:
    f.close()

plotBinnedCrossCor(xi,rho,sig,os)
plt.savefig(results_dir+"/"+name+"ORF.jpg")
plt.close()

