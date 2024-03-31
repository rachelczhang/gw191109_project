import kick as kc
import binary_sampling as binsam
import numpy as np
import os
import sys
from scipy.stats import maxwell

import random
from datetime import datetime
random.seed(datetime.now())

# Anomaly
# input: - orbital eccentricity; - tolerance for convergence
# output: - eccentric anomaly ; - true anomaly
def anomaly(e, tol):

    tol = 10 ** (-3)

    mean = 2.0 * np.pi * random.random()
    ecc1, ecc2 = 0., 2.0 * np.pi
    effe1, effe2 = ecc1 - e * np.sin(ecc1) - mean, ecc2 - e * np.sin(ecc2) - mean

    while abs(effe1 - effe2) > tol:
        ecc3 = (ecc1 + ecc2) / 2.0
        effe3 = ecc3 - e * np.sin(ecc3) - mean
        if (effe3 > 0 and effe1 > 0) or (effe3 < 0 and effe1 < 0):
            ecc1, effe1 = ecc3, effe3
        else:
            ecc2, effe2 = ecc3, effe3

    gammaa = ecc1 # eccentric anomaly
    theta = 2.0 * np.arctan(np.sqrt((1.0 + e) / (1.0 - e)) * np.tan(ecc1 / 2.0)) # true anomaly

    return(
	    gammaa,
	    theta,
	)

# Change in orbital elements after SN kick
# input: - initial mass [Msun]; - final mass [Msun]; - companion mass [Msun]; - semi-major axis [AU]; - eccentricity; - kick velocity [km/s]
# input optional : - mass third companion[Msun]; - outer semi-major axis [AU]; - outer eccentricty; - orbital inclination [rad]
# output : - new semi-major axis [AU], eccentricity; - tilt angle [rad]; - new outer semi-major axis [AU], eccentricity, inclination
def SNkick(m1old, m1new, m2, a, e, vkick, m3=0., a3=0., e3=0., i3=0.):

    # scale vkick so that G=1, Msun=1, AU=1

    vkick = vkick / 30.0

    # definitions important quantities

    mbold = m2 + m1old
    mbnew = m2 + m1new
    deltam = m1old - m1new
    muold = m2 / mbold
    munew = m2 / mbnew
    vorb = np.sqrt(mbold / a)

    phi = np.arcsin(2.0 * random.random() - 1.0)
    omega = 2.0 * np.pi * random.random()
    vk = np.array([vkick * np.cos(phi) * np.sin(omega), -vkick * np.cos(phi) * np.cos(omega), vkick * np.sin(phi)])
    vkmod = np.sqrt(kc.scalar(vk, vk))
    alpha = 2.0 * np.pi * random.random()

    # anomaly

    tol = 10 ** (-3)
    gammaa, theta = anomaly(e, tol)

    # r and v vectors

    rr = np.array([0., -a * (1.0 - e ** 2.0) / (1.0 + e * np.cos(theta)), 0.])
    rrmod = np.sqrt(kc.scalar(rr, rr))

    v0 = vorb * np.sqrt((1.0 + e * np.cos(gammaa)) / (1.0 - e * np.cos(gammaa)))
    x = a * np.sqrt(1.0 - e ** 2.0) * np.cos(gammaa) * np.cos(theta) + a * np.sin(gammaa) * np.sin(theta)
    y = -a * np.sqrt(1.0 - e ** 2.0) * np.cos(gammaa) * np.sin(theta) + a * np.sin(gammaa) * np.cos(theta)
    vx = v0 * x / np.sqrt(x ** 2.0 + y ** 2.0)
    vy = v0 * y / np.sqrt(x ** 2.0 + y ** 2.0)
    vv = np.array([vx, vy, 0.])
    vvn = np.add(vv, vk)
    vvnmod = np.sqrt(kc.scalar(vvn, vvn))

    # new sma, ecc, tilt

    anew = 2.0 / rrmod - vvnmod ** 2.0 / mbnew
    anew = 1.0 / anew

    rcv = kc.cross(rr, vvn)
    rcvmod = np.sqrt(kc.scalar(rcv, rcv))
    enew = rcvmod ** 2.0 / mbnew / anew
    enew = np.sqrt(1.0 - enew)

    rcvold = kc.cross(rr, vv)
    rcvoldmod = np.sqrt(kc.scalar(rcvold, rcvold))
    arg = kc.scalar(rcv, rcvold) / rcvmod / rcvoldmod
    tilt = float('%s' % float('%.5g' % arg))
    tilt = np.arccos(tilt)

	# compute shift in the outer orbit due to SN in the inner if hierarchy (m3>0)

    if m3 > 0.0:

        # useful quantities

        mtotold = mbold + m3
        mtotnew = mbnew + m3
        vorb3 = np.sqrt(mtotold / a3)

        # shift

        vsys = np.array([(1.0 - munew) * (vx * (muold - munew) / (1.0 - munew) + vk[0]), (1.0 - munew) * (vy * (muold - munew) / (1.0 - munew) + vk[1]), (1.0 - munew) * vk[2]])
        deltarr = np.array([0., (munew - muold) * a * (1.0 - e ** 2.0) /(1.0 + e * np.cos(theta)), 0.])

        # anomaly outer orbit

        tol = 10 ** (-3)
        gammaa3, theta3 = anomaly(e3, tol)

        # r3 and v3 vectors

        rr3 = np.array([a3 * (1.0 - e3 ** 2.0) / (1.0 + e3 * np.cos(theta3)) * np.cos(i3) * np.sin(alpha), -a3 * (1.0 - e3 ** 2.0) / (1.0 + e3 * np.cos(theta3)) * np.cos(i3) * np.cos(alpha), a3 * (1.0 - e3 ** 2.0) / (1.0 + e3 * np.cos(theta3)) * np.sin(i3)])
        rr3mod = np.sqrt(kc.scalar(rr3, rr3))

        v03 = vorb3 * np.sqrt((1.0 + e3 * np.cos(gammaa3)) / (1.0 - e3 * np.cos(gammaa3)))
        x3 = a3 * np.sqrt(1.0 - e3 ** 2.0) * np.cos(gammaa3) * np.cos(theta3) + a3 * np.sin(gammaa3) * np.sin(theta3)
        y3 = -a3 * np.sqrt(1.0 - e3 ** 2.0) * np.cos(gammaa3) * np.sin(theta3) + a3 * np.sin(gammaa3) * np.cos(theta3)
        xx3 = x3 * np.cos(alpha) - y3 * np.cos(i3) * np.sin(alpha)
        yy3 = x3 * np.sin(alpha) + y3 * np.cos(i3) * np.cos(alpha)
        zz3 = y3 * np.sin(i3)
        vv3 = np.array([v03 * xx3 / np.sqrt(xx3 ** 2.0 + yy3 ** 2.0 + zz3 ** 2.0), v03 * yy3 / np.sqrt(xx3 ** 2.0 + yy3 ** 2.0 + zz3 ** 2.0), v03 * zz3 / np.sqrt(xx3 ** 2.0 + yy3 ** 2.0 + zz3 ** 2.0)])
        vsum = np.add(vv3, vsys)
        vsummod = np.sqrt(kc.scalar(vsum, vsum))

        rr3n = np.add(rr3, deltarr)
        rr3nmod = np.sqrt(kc.scalar(rr3n, rr3n))

        # new sma, ecc, incl 

        a3new = 2.0 / rr3nmod - vsummod ** 2.0 / mtotnew
        a3new = 1.0 / a3new

        rcv3 = kc.cross(rr3n, vsum)
        rcv3mod = np.sqrt(kc.scalar(rcv3, rcv3))
        e3new = rcv3mod ** 2.0 / mtotnew / a3new
        e3new = np.sqrt(1.0 - e3new)

        i3new = np.arccos(kc.scalar(rcv, rcv3) / rcvmod / rcv3mod)

    else:

        a3new = 0.0
        e3new = 0.0
        i3new = 0.0

    return(
        anew,
        enew,
        tilt,
        a3new,
        e3new,
        i3new,
    )




