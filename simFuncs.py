#Script to create simulated data

import os, glob
import numpy as np
import libstempo as T2
from libstempo import toasim as LT
from libstempo import plot as LP
import enterprise
from enterprise.pulsar import Pulsar


def create_psrs(parfiles, obstimes, toaerr=0.4):
    psrs = []
    for file in parfiles:
        psrs.append(LT.fakepulsar(parfile=file, obstimes=obstimes, toaerr=toaerr))
    return psrs

def add_gwb(psrs, amp=1e-14, gam=13./3.):
    """ Add a GWB.
        Takes in pulsar objects, amplitude of RN, and spectral index gamma.
    """
    LT.createGWB(psrs, amp, gam) #modifies pulsars in place, no need to return anything
    
def add_cgw(psrs, pdict, tref, iters):
    """ Add a continuous wave signal.
        Takes in pulsar objects, a parameter dictionary for the single source, and
        a reference time for observations.
    """
    for psr in psrs:
        #also modifying pulsars in place
        LT.add_cgw(psr, gwtheta=pdict['gwtheta'], gwphi=pdict['gwphi'], mc=pdict['mc'], dist=pdict['dist'], fgw=pdict['fgw'], phase0=pdict['phase0'], psi=pdict['psi'], inc=pdict['inc'], pdist=1., pphase=None, psrTerm=False, evolve=False, phase_approx=False, tref=tref)
        
        #iterate the timing model fit a few times
        try:
            psr.fit(iters=iters)
            print(psr.name)
        except:
            print(psr.name, 'had timing model fit issue. Excluding from PTA.')

def lt2ent(psrs):
    """ Converts libstempo pulsar objects to enterprise pulsar objects.
    """
    Epsrs = [Pulsar(psr) for psr in psrs]
    return Epsrs
    
def save_sims(psrs, outdir):
    """ Save simulated timing files.
    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for p in psrs:
        p.savetim(outdir + p.name + '_simulated.tim')
