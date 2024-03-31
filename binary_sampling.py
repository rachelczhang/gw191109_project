import numpy as np
import kick as kc

import random
from datetime import datetime
random.seed(datetime.now())

# ---------- constants (units AU=1, Msun=1, G=1)

pi = 3.14159
pc = 2.0*10**5      # pc in AU
vsc = 30.0          # scale of velocity 30 km/s
c = 10**4           # velocity of light

dyear = 365.24	  # from yr to day

bper = 1.0 - 0.55	  # Sana+ 2012 parameters
sanamin = 0.15
sanamax = 5.5

#---------- calculations

# binary parameters
# input min max slope mass function, # mass ratio = 0 independent , = 1 uniform, # binary eccentricity = 1 thermal , = 2 circular, = 3 uniform, # min max sma, # binary sma = 1 log-uniform , = 2 Sana+ 2012 , = 3 uniform
# output mass1, mass2, sma (AU), period (day), eccentricity

def bin_sample(mmin, mmax, slope, iratio, ife, amin, amax, ifa):

    mb1 = kc.sample_mass(mmin, mmax, slope)

    mb2 = 0
    if iratio == 0: # mass ratio = 0 independent , = 1 uniform
        mb2 = kc.sample(mmin, mmax, slope)
    elif iratio == 1:
        while mb2 < mmin:
            mb2 = random.random()
            mb2 = mb2 * mb1
    if(mb2 > mb1):
        mbb = mb1
        mb1 = mb2
        mb2 = mbb
    mb12 = mb1+mb2

    if ife == 1:  # binary eccentricity = 1 thermal , = 2 circular, = 3 uniform
        e12 = random.random()
        e12 = np.sqrt(e12)
    elif ife == 2:
        e12 = 0
    elif ife == 3:
        e12 = random.random()

    if ifa == 1:  # binary sma = 1 log-uniform , = 2 Sana+ 2012 , = 3 uniform
        sem12 = random.random()
        sem12 = amin * np.exp(sem12 * np.log(amax / amin))
        per12 = dyear * (sem12 / mb12) ** 1.5
    elif ifa == 2:
        per12 = random.random()
        per12 = (per12 * (sanamax ** bper - sanamin ** bper) +
                 sanamin ** bper) ** (1. / bper)
        per12 = 10 ** per12
        sem12 = mb12 * (per12 / dyear) ** (2.0 / 3.0)
    elif ifa == 3:
        sem12 = random.random()
        sem12 = amin + sem12 * (amax - amin)
        per12 = dyear * (sem12 / mb12) ** 1.5

    return(
        mb1,
        mb2,
        sem12,
        per12,
        e12,
    )
