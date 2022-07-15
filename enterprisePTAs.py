import numpy as np

import enterprise

from enterprise.pulsar import Pulsar
import enterprise.signals.parameter as parameter
from enterprise.signals import utils
from enterprise.signals import signal_base
from enterprise.signals import selections
from enterprise.signals.selections import Selection
from enterprise.signals import white_signals
from enterprise.signals import gp_signals
from enterprise.signals import gp_priors

from enterprise_extensions import model_utils

from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc

def crnPTA(psrs,fixedGamma=True):
    Tspan = model_utils.get_tspan(psrs)

    efac = parameter.Constant(1.0)
    ef = white_signals.MeasurementNoise(efac=efac)

    log10_A_gw = parameter.Uniform(-18,-12)('gw_log10_A')
    if fixedGamma:
        gamma_gw = parameter.Constant(13.0/3.0)('gw_gamma')
    else:
        gamma_gw = parameter.Uniform(0,7)('gw_gamma')
        
    plgw = utils.powerlaw(log10_A=log10_A_gw, gamma=gamma_gw)
    crn = gp_signals.FourierBasisGP(plgw, components=14, Tspan=Tspan, name='gw')

    tm = gp_signals.TimingModel(use_svd=True)

    model = ef + tm  + crn

    # initialize PTA
    pta = signal_base.PTA([model(psr) for psr in psrs])
    return pta

def mcSample(pta,n=1000000,outDir='./PerFreq2A/',resume=False):
    # dimension of parameter space
    ndim = len(pta.param_names)

    # initial jump covariance matrix
    cov = np.diag(np.ones(ndim) * 0.01**2)

    # set up jump groups by red noise groups 
    groups  = [range(0, ndim)] #CHECK

    # intialize sampler
    sampler = ptmcmc(ndim, pta.get_lnlikelihood, pta.get_lnprior, cov, groups=groups,
                     outDir=outDir,resume=False)

    # sampler for N steps
    x0 = np.hstack(p.sample() for p in pta.params)
    
    sampler.sample(x0, n, SCAMweight=30, AMweight=15, DEweight=50)
