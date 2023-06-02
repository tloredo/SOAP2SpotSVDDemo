"""
Experiment with SVD using SOAP-2 data in the Ca H & K line region.

Created 2023-05-26 by Tom Loredo
Based on a script from May 2017, subsequently revised in 2019, 2020, 2021.
"""

import numpy as np
from numpy import *
import scipy
from scipy.sparse import linalg
import matplotlib as mpl
from matplotlib.pyplot import *

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
from matplotlib import colors as mcolors

import colorcet as cc

from soap2_data import prep_data, fetch_ca_spec, fetch_full_spec
from outerproduct import OuterProduct


# This imports Tom's plotting defaults
try:
    import myplot
    from myplot import ipy_ion, close_all, csavefig
    ipy_ion()
    #myplot.tex_on()
    csavefig.save = False
except ImportError:
    ion()

# If having trouble with plots, try switching the backend:
# matplotlib.use('TkAgg')
# mpl.style.use('classic')


# Fetch the data:
fetcher = prep_data('SOAP2-1Spot')

# Full spectrum, just 4 phases:
full_spec = fetch_full_spec(fetcher)

# Zoomed into the calcium H & K line region, 100 phases:
ca_spec = fetch_ca_spec(fetcher)


# We'll use the many-phase zoomed-in spectrum here,
# aliased so we can swap in other spectra later.
dynspec = ca_spec

# Compute average spectrum, and an "image" of the avg-subtracted 
# spectra (with time as the vertical dimension), rescaled.
avg = dynspec.active.sum(0) / dynspec.nphases
delta_image = dynspec.active - avg
l, u = delta_image.min(), delta_image.max()
ll, lu = dynspec.lambdas[0], dynspec.lambdas[-1]
delta_image = (delta_image - l)/(u - l)


# Lines in this region of the spectrum:
# For info on air/vacuum conversino, see:
# Spectra - SDSS-III
# http://www.sdss3.org/dr8/spectro/spectra.php
# Atomic Data for Resonance Absorption Lines... - ADS
# https://ui.adsabs.harvard.edu/abs/1991ApJS...77..119M/abstract
# Ca II K line, air & vacuum:
CaK = [3933.663, 3934.777]
# Ca II H line, air & vacuum:
CaH = [3968.468, 3969.591]
lines = CaK + CaH
lines = (CaK[0], CaH[0])  # air only; appears to be what SOAP uses


# Plot the quiet-star spectrum.
if True:
    fig = figure('quiet-full', figsize=(12,3.5))
    fig.subplots_adjust(left=.1, right=.9, bottom=0.21, top=.9)
    title('Quiet Sun integrated spectrum')
    plot(dynspec.lambdas, dynspec.quiet, label='Integrated spectrum no activity at full resolution')
    xlabel(r'$\lambda$ ($\AA$)')
    ylabel('Relative flux')
    for line in lines:
        axvline(line, c='r', lw=.75, alpha=.5)
    # Label vacuum lines:
    text(CaK[0]+.7, 6200, 'Ca II K', horizontalalignment='left',
         verticalalignment='top')
    text(CaH[0]+.7, 6200, 'Ca II H', horizontalalignment='left',
         verticalalignment='top')


# Plot spectra at a few phases.
if False:
    fig = figure('phases', figsize=(12,3.5))
    fig.subplots_adjust(left=.1, right=.9, bottom=0.21, top=.9)
    title('Selected phases')
    for i in range(0, dynspec.nphases, dynspec.nphases//4):
        plot(dynspec.lambdas, dynspec.active[i,:], label='$\phi = {:0.2f}$'.format(dynspec.phases[i]),
             lw=1, alpha=.5)
    xlabel(r'$\lambda$ ($\AA$)')
    ylabel('Relative flux')
    legend()
    for line in lines:
        axvline(line, c='r', lw=.75, alpha=.5)

    # Plot mean-subtracted spectra at a few phases.
    fig = figure('diff', figsize=(12,3.5))
    fig.subplots_adjust(left=.1, right=.9, bottom=0.21, top=.9)
    title('Mean-subtracted spectra')

    # If there are many phases, plot a subset.
    if dynspec.nphases >= 10:
        step = dynspec.nphases//10
    else:
        step = 1

    for i in range(0, dynspec.nphases, step):
        plot(dynspec.lambdas, dynspec.active[i,:] - avg, label='$\phi = {:0.2f}$'.format(dynspec.phases[i]),
             lw=1.5, alpha=.5)
    xlabel(r'$\lambda$ ($\AA$)')
    ylabel('Difference')
    legend()
    for line in lines:
        axvline(line, c='r', lw=.75, alpha=.5)


# Try a waterfall-style plot of spectra at a few phases.
if True:
    fig = figure(figsize=(8, 4))
    # This breaks after mpl-3.6; use add_subplot instead:
    # ax = fig.gca(projection='3d')
    ax = fig.add_subplot(projection='3d')

    verts = []

    # If there are many phases, plot a subset.
    if dynspec.nphases >= 10:
        step = dynspec.nphases//10
    else:
        step = 1

    ind = list(range(0, dynspec.nphases, step))
    zs = dynspec.phases[ind]
    y_l, y_u = 0., 0.
    for i, z in zip(ind, zs):
        # ys = dynspec.active[i,:] - avg
        ys = dynspec.active[i,:]
        y_l = min(y_l, ys.min())
        y_u = max(y_u, ys.max())
        # ys[0], ys[-1] = 0, 0
        verts.append(list(zip(dynspec.lambdas, ys)))

    def cc2(arg):
        return mcolors.to_rgba(arg, alpha=0.6)

    poly = PolyCollection(verts, facecolors=[cc2('r'), cc2('g'), cc2('b'),
                                             cc2('y')])
    poly.set_alpha(0.7)
    ax.add_collection3d(poly, zs=zs, zdir='y')

    ax.set_xlabel(r'$\lambda$ ($\AA$)')
    ax.set_xlim3d(dynspec.lambdas[0], dynspec.lambdas[-1])
    ax.set_ylabel('$\phi$')
    ax.set_ylim3d(-.5, .5)
    ax.set_zlabel(r'$\Delta$')
    ax.set_zlim3d(y_l, y_u)


# Plot an image of the difference of spectra from the time-averaged spectrum.
# *** Double-check the origin='lower' setting here (matches use in OuterProduct,
# which was verified w/ test case).
if False:
    if True:  # larger labels
        rc('axes', labelsize=18)
        rc('xtick.major', pad=8)
        rc('xtick', labelsize=14)
        rc('ytick.major', pad=8)
        rc('ytick', labelsize=14)
        rc('figure.subplot', bottom=.19, top=.925)

    # fig = figure('diff-img', figsize=(12,5))  # orig size
    fig = figure('diff-img', figsize=(12,3.5))  # shorter for proposal
    imshow(delta_image, cmap=cc.cm.bwy_r, interpolation='nearest', aspect='auto',
           origin='lower', extent=(ll, lu, dynspec.phases[0], dynspec.phases[-1]))
    xlim(ll, lu)
    for line in lines:
        # axvline(line, c='r', lw=.75, alpha=.5)
        axvline(line, c='w', lw=1, alpha=.5)
    xlabel(r'$\lambda$ ($\AA$)')
    ylabel("Spot rotational phase")
    title('Difference from time-averaged spectrum')

    # Zoomed for inset:
    # fig = figure('diff-img2', figsize=(12,5))
    fig = figure('diff-img2', figsize=(12,3.5))
    imshow(delta_image, cmap=cc.cm.bwy_r, interpolation='nearest', aspect='auto',
           origin='lower', extent=(ll, lu, dynspec.phases[0], dynspec.phases[-1]))
    xlim(ll, lu)
    for line in lines:
        # axvline(line, c='r', lw=.75, alpha=.5)
        axvline(line, c='w', lw=1, alpha=.5)
    xlim(3940., 3947.)
    ylim(-.25, .25)
    # xlabel(r'$\lambda$ ($\AA$)')
    #ylabel("Spot rotational phase")
    # title('Difference from time-averaged spectrum')


# Test of outer product grid alignment:
if False:
    spec_x = linspace(1., 5., 5)
    spec_y = array([0., 1., 0., -1., 0.])
    temp_x = linspace(0., 10., 6)
    temp_y = array([-2., -1., 0., 1., 2., 1.])
    outer = OuterProduct(spec_x, spec_y, temp_x, temp_y)
    f5 = outer.plot()


# Compute and plot evolving spectral components via SVD.
if True:

    # Compute top 3 singular values and vectors; the largest singular value
    # has index 2 (i.e., they are increasing in s[:]).
    # U[phi,k] contains the phase basis function for singular value # k
    # Vt[k,lambda] contains the wavelength basis function for singular value # k
    U, s, Vt = linalg.svds(delta_image, k=3)

    outer = OuterProduct(dynspec.lambdas, Vt[-1,:], dynspec.phases, U[:,-1])
    f1 = outer.plot(lines)
    f1.axes[1].set_title('1st singular value')

    outer = OuterProduct(dynspec.lambdas, Vt[-2,:], dynspec.phases, U[:,-2])
    # f2 = outer.plot(lines)
    f2 = outer.plot(lines, zoom=((3959.5, 3965.5), (-.24, .24)))
    f2.axes[1].set_title('2nd singular value')

    outer = OuterProduct(dynspec.lambdas, Vt[-3,:], dynspec.phases, U[:,-3])
    # f3 = outer.plot(lines)
    f3 = outer.plot(lines, zoom=((3946.5, 3950.), (-.24, .24)))
    f3.axes[1].set_title('3rd singular value')

    # *** DANGER: Changing xlim does not appear to affect the image as expected.
    # Just pass subsets of the data to zoom.
    if False:
        # Zoomed version of SVD3:
        outer = OuterProduct(dynspec.lambdas, Vt[-3,:], dynspec.phases, U[:,-3])
        f4 = outer.plot(lines)
        f4.axes[1].set_title('3rd singular value')
        xrng = (3946.5, 3950.)
        xrng = (3947.25, 3947.9)
        f4.axes[0].set_xlim(*xrng)
        f4.axes[1].set_xlim(*xrng)
        spec = dynspec.active[0,:]
        spec = spec/spec.max()
        l, u = Vt[-3,:].min(), Vt[-3,:].max()
        spec = spec*(u-l) + l
        f4.axes[1].plot(dynspec.lambdas, spec, 'g--', lw=1, alpha=.5)
        f4.axes[0].axvline(3948.05)
        f4.axes[1].axvline(3948.05)

    # Zoomed version of SVD3:
    i1, i2 = 7544, 7753  # very narrow region, for checking alignments
    i1, i2 = 7520, 7890
    lambdas = dynspec.lambdas[i1:i2]
    vvals = Vt[-3,i1:i2]
    outer = OuterProduct(lambdas, vvals, dynspec.phases, U[:,-3])
    f5 = outer.plot()
    spec = dynspec.active[0,:]
    spec = spec/spec.max()
    l, u = Vt[-3,:].min(), Vt[-3,:].max()
    # spec = spec*(u-l) + l  # v. center the spec curve
    spec = 1.2*u*spec  # keep spec curve > 0
    f5.axes[1].plot(dynspec.lambdas, spec, 'g--', lw=1, alpha=.5)
    f5.axes[1].axhline(0., ls='-', c='gray', lw=0.75)
    f5.axes[1].set_title('3rd singular value')
