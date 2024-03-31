import numpy as np
import os
import sys

import random
from datetime import datetime
random.seed(datetime.now())

# Scalar product
# input: - vector 1, vector 2
# output: - scalar product


def scalar(v1, v2):

    scalar_product = np.dot(v1, v2)

    return scalar_product

# Cross product
# input: - vector 1, vector 2
# output: - vector cross product


def cross(v1, v2):

    cross_product = []
    cross_product.append(v1[1] * v2[2] - v1[2] * v2[1])
    cross_product.append(v1[2] * v2[0] - v1[0] * v2[2])
    cross_product.append(v1[0] * v2[1] - v1[1] * v2[0])
    cross_product = np.array(cross_product)

    return cross_product

# Rotation assuming isotropy on a sphere
### input: - vector
# output: - vector rotated


def rotation_iso(vec):

    phi = 2.0 * np.pi * random.random()
    theta = np.arccos(1.0 - 2.0 * random.random())

    result = []
    result.append(vec[0] * np.sin(theta) * np.cos(phi))
    result.append(vec[0] * np.sin(theta) * np.sin(phi))
    result.append(vec[0] * np.cos(theta))
    result = np.array(result)

    return result

# Rotation of a vector
### input: - vector
# output: - vector rotated


def rotation(vec):

    phi = 2.0 * np.pi * random.random()
    theta = np.pi * random.random()

    result = []
    result.append(vec[0] * np.sin(theta) * np.cos(phi))
    result.append(vec[0] * np.sin(theta) * np.sin(phi))
    result.append(vec[0] * np.cos(theta))
    result = np.array(result)

    return result

# Sample spin
# input: - spin sampling mode, maximum spin, mass
### output: - spin


def sample_spin(ispin, sspin, mb):

    spin = []
    if ispin == 0:
        spin = 0.0
    elif ispin == 1:
        spin = sspin * random.random()
    elif ispin == 2:
        spin = sspin
    elif ispin == 3:
        spinup = 0.33 * np.tanh(38.0-mb) + 0.59
        spindown = 0.40 * np.tanh(21.0-mb) + 0.40
        spin = random.random() * (spinup - spindown) + spindown

    return spin

# Sample mass assuming power law
# input: - min mass, max mass, slope
### output: - mass


def sample_mass(mmin, mmax, slope):

    gammaa = 1.0 - slope

    if slope == 1:
        mb = mmin * np.exp(random.random() * np.log(mmax / mmin))
    else:
        mb = (random.random() * (mmax ** gammaa - mmin **
                                 gammaa) + mmin ** gammaa) ** (1.0 / gammaa)

    return mb

# Kick after GW merger using Lousto+ 2010, 2011, 2012
# input: - mass 1, mass2, spin1, spin2, ecc
# output: - kick velocity, effective spin, chip, final mass, final spin


def kick(mb1, mb2, s1, s2, ecc):

    # remove the following when done
    #ffile = open('input.txt', 'w')

    # params for kick velocity

    A = 1.2*10**4  # km/s
    B = -0.93
    H = 6.9*10**3  # km/s
    xi = (145.0/180.0)*np.pi
    V11 = 3678.0  # km/s
    VA = 2481.0  # km/s
    VB = 1793.0  # km/s
    VC = 1507.0  # km/s

    ex = np.array([1.0, 0.0, 0.0])
    ey = np.array([0.0, 1.0, 0.0])
    ez = np.array([0.0, 0.0, 1.0])

    # remove the following when done
    #print(mseed, mb2, ispin, s1, sspin, file=ffile)

    mmb1, mmb2 = max(mb1, mb2), min(mb1, mb2)
    mb1, mb2 = mmb1, mmb2

    ratio12 = mb2 / mb1
    eta = ratio12 / (1.0+ratio12) ** 2.0

    # sampling spin BH1 if rep = 0 and BH2

    spin1 = np.array([s1, 0.0, 0.0])
    mod_spin1 = s1
    spin1 = rotation_iso(spin1)

    spin2 = np.array([s2, 0.0, 0.0])
    mod_spin2 = s2
    spin2 = rotation_iso(spin2)

    # phase BH1-BH2 and direction of merger

    iota0 = 2.0 * np.pi * random.random()

    dirmer = np.array([1.0, 0.0, 0.0])
    dirmer = rotation(dirmer)

    # calculating vectors

    a1p = scalar(spin1, ez)
    a2p = scalar(spin2, ez)

    if mod_spin1 == 0 and mod_spin2 == 0:
        iotaD = 0.0
        diffs = spin2 - ratio12 * spin1
        diffsp = cross(diffs, ez)
        mod_diffsp = np.sqrt(scalar(diffsp, diffsp))
    else:
        Delta = spin2 - ratio12 * spin1
        Deltap = cross(Delta, ez)
        mod_Deltap = np.sqrt(scalar(Deltap, Deltap))
        iotaD = np.arccos(scalar(Deltap, dirmer) / mod_Deltap)
        mod_diffsp = mod_Deltap

    stilde = 2.0 * (spin2 + spin1 * ratio12 ** 2.0) / (1.0 + ratio12) ** 2.0
    stildepar = scalar(stilde, ez)

    # kick velocity

    vm = A * (eta ** 2.0) * np.sqrt(1.0 - 4.0 * eta) * (1 + B * eta)
    vperp = H * (eta ** 2.0) / (1.0 + ratio12) * (a2p - ratio12 * a1p)
    vpar = (16.0 * (eta ** 2.0) / (1.0 + ratio12)) * (V11 + VA * stildepar + VB *
                                                      stildepar ** 2.0 + VC * stildepar ** 3.0) * mod_diffsp * np.cos(iotaD - iota0)

    vk1 = vm + vperp * np.cos(xi)
    vk2 = vperp * np.sin(xi)
    vk3 = vpar
    mod_vk = (1.0 + ecc) * np.sqrt(vk1 ** 2.0 + vk2 ** 2.0 + vk3 ** 2.0)

    chieff = (mb1 * a1p + mb2 * a2p) / (mb1 + mb2)

    if s1 > 0.0:
        sin1 = np.arccos(a1p / s1)
    else:
        sin1 = 0.0
    if s2 > 0.0:
        sin2 = np.arccos(a2p / s2)
    else:
        sin2 = 0.0
    kq = ratio12 * (4.0 * ratio12 + 3.0) / (4.0 + 3.0 * ratio12)
    chip = max(s1 * np.sin(sin1), kq * s2 * np.sin(sin2))

    #########

    if mod_spin1 == 0 or mod_spin2 == 0:
        aalpha = 0.0
    else:
        aalpha = np.arccos(scalar(spin1,spin2)/mod_spin1/mod_spin2)
    if mod_spin1 == 0:
        bbeta = 0.0
    else:
        bbeta = np.arccos(a1p/mod_spin1)
    if mod_spin2 == 0:
        ggamma = 0.0
    else:
        ggamma = np.arccos(a2p/mod_spin2)

    mfin = massrem(mb1, mb2, s1, s2)
    sfin = spinrem(mb1, mb2, s1, s2, aalpha, bbeta, ggamma)

    #########

    return(
        mod_vk,
        chieff,
        chip,
        mfin,
        sfin,
    )

# Remnant spin following Jimenez-Forteza+ PHYSICAL REVIEW D 95, 064024 (2017)
# input: - mass 1, mass 2 , spin 1, spin 2, angle spin1-spin2, angle Lorb-spin1, angle Lorb-spin2
# output: - remnant spin


def spinrem(mb1, mb2, s1, s2, aalpha, bbeta, ggamma):

    ratio12 = mb2 / mb1
    eta = ratio12 / (1.0+ratio12) ** 2.0

    shat = ((mb1 ** 2.0) * s1 + (mb2 ** 2.0) * s2) / (mb1 ** 2.0 + mb2 ** 2.0)
    stot = ((mb1 ** 2.0) * s1 + (mb2 ** 2.0) * s2) / (mb1 + mb2) ** 2.0

    a2, a3, a5 = 3.833, -9.49, 2.513
    lorb1 = (1.3 * a3 * eta ** 3.0 + 5.24 * a2 * eta ** 2.0 +
             2.0 * np.sqrt(3.0) * eta) / (2.88 * a5 * eta + 1.0)

    lorb2 = 0.68637

    b1, f10, f11, f12 = 1.00096, 0.0, 4.4092, 0.512
    f13 = 64.0 - 64.0 * f10 - 16.0 * f11 - 4.0 * f12
    b2, f20, f21, f22 = 0.788, 0.0, 8.774, -32.1
    f23 = 64.0 - 64.0 * f20 - 16.0 * f21 - 4.0 * f22
    b3, f30, f31, f32 = 0.654, 0.0, 22.83, -154.0
    f33 = 64.0 - 64.0 * f30 - 16.0 * f31 - 4.0 * f32
    b5, f50, f51, f52 = 0.840, 1.8805, -4.77, 0.0
    f53 = 64.0 - 64.0 * f50 - 16.0 * f51 - 4.0 * f52
    b1 = b1 * (f10 + f11 * eta + f12 * eta ** 2.0 + f13 * eta ** 3.0)
    b2 = b2 * (f20 + f21 * eta + f22 * eta ** 2.0 + f23 * eta ** 3.0)
    b3 = b3 * (f30 + f31 * eta + f32 * eta ** 2.0 + f33 * eta ** 3.0)
    b5 = b5 * (f50 + f51 * eta + f52 * eta ** 2.0 + f53 * eta ** 3.0)
    lorb3 = 0.68637 + (0.00954 * b3 * shat ** 3.0 + 0.0851 * b2 *
                       shat ** 2.0 - 0.194 * b1 * shat) / (1.0 - 0.579 * b5 * shat)

    lorb = lorb1 - lorb2 + lorb3

    d10, d11, d20, d30, d31 = 0.322, 9.33, -0.0598, 2.32, -3.26
    deltachi = s1 - s2

    a1 = d10 * (1.0 - 4.0 * eta) ** 0.5 * eta ** 2.0 * (d11 * eta + 1.0)
    a2 = d20 * eta ** 3.0
    a3 = d30 * (1.0 - 4.0 * eta) ** 0.5 * eta ** 3.0 * (d31 * eta + 1.0)
    deltal = a1 * deltachi + a2 * deltachi ** 2.0 + a3 * shat * deltachi

    lorb = lorb + deltal

    #srem = lorb + stot # parallel/anti-parallel spins
    lorb = lorb / eta
    srem = np.sqrt(s1 ** 2 + s2 ** 2 * (ratio12 ** 4) + 2.0 * s1 * s2 * (ratio12 ** 2) * np.cos(aalpha) + 2 * (s1 * np.cos(bbeta) + s2 * (ratio12 ** 2) * np.cos(ggamma)) * lorb * ratio12 + lorb ** 2 * (ratio12 ** 2)) / (1 + ratio12) ** 2

    return srem

# Remnant mass following Jimenez-Forteza+ PHYSICAL REVIEW D 95, 064024 (2017)
# input: - mass 1, mass 2
# output: - remnant mass


def massrem(mb1, mb2, s1, s2):

    ratio12 = mb2 / mb1
    eta = ratio12 / (1.0+ratio12) ** 2.0

    shat = ((mb1 ** 2.0) * s1 + (mb2 ** 2.0) * s2) / (mb1 ** 2.0 + mb2 ** 2.0)

    a2, a3, a4 = 0.5610, -0.847, 3.145
    erad1 = a4 * eta ** 4.0 + a3 * eta ** 3.0 + a2 * \
        eta ** 2.0 + (1.0 - 2.0 * np.sqrt(2.0) / 3.0) * eta

    erad3 = 0.0484161

    b1, b2, b3, b5 = -0.209, -0.197, -0.159, 2.985
    f10, f11 = 1.81, 15.7
    f12 = 16.0 - 16.0 * f10 - 4.0 * f11
    f20, f21 = 4.27, 0.0
    f22 = 16.0 - 16.0 * f20 - 4.0 * f21
    f30, f31 = 31.09, 243.6
    f32 = 16.0 - 16.0 * f30 - 4.0 * f31
    f50, f51 = 1.56735, -0.58
    f52 = 16.0 - 16.0 * f50 - 4.0 * f51
    b1 = b1 * (f10 + f11 * eta + f12 * eta ** 2.0)
    b2 = b2 * (f20 + f21 * eta + f22 * eta ** 2.0)
    b3 = b3 * (f30 + f31 * eta + f32 * eta ** 2.0)
    b5 = b5 * (f50 + f51 * eta + f52 * eta ** 2.0)
    erad2 = 0.0484161 / (1.0 - 0.212 * b5 * shat) * (0.128 * b3 *
                                                     shat ** 3.0 + 0.211 * b2 * shat ** 2.0 + 0.346 * b1 * shat + 1.0)

    erad = erad1 * erad2 / erad3

    d10, d11, d20, d30, d31 = -0.098, -3.23, 0.0112, -0.0198, -4.92
    deltachi = s1 - s2

    a1 = d10 * (1.0 - 4.0 * eta) ** 0.5 * eta ** 2.0 * (d11 * eta + 1.0)
    a2 = d20 * eta ** 3.0
    a3 = d30 * (1.0 - 4.0 * eta) ** 0.5 * eta * (d31 * eta + 1.0)
    deltae = a1 * deltachi + a2 * deltachi ** 2.0 + a3 * shat * deltachi

    erad = erad + deltae

    mrem = (mb1 + mb2) * (1.0 - erad)

    return mrem
