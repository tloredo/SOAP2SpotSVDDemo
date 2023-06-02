"""
Provide a class for computing and displaying the outer product of two
tabulated functions.

2023-06-02:  Refactored from earlier `soap-process` scripts
"""
import numpy as np
from numpy import *

import matplotlib as mpl
from matplotlib.pyplot import *
import colorcet as cc


class OuterProduct:
    """
    Compute and display the outer product of two tabulated functions, as
    components of a reduced-rank approximation of a 2-D function.
    """

    def __init__(self, xvals, uvals, yvals, vvals):
        """
        Define a 2-D function f(x,y) as the outer product of separate
        functions of x and y:

            f(x,y) = u(x)*v(y)
        """
        self.xvals, self.uvals = xvals, uvals
        self.yvals, self.vvals = yvals, vvals
        self.outer = np.outer(uvals, vvals)

    def plot(self, lines=None, zoom=None):
        """
        Plot the outer product as an image, with the component functions
        plotted in the upper and right margins.
        """
        # Set up a figure with the 3 sets of axes.
        fig = figure(figsize=(10,7))
        gs = GridSpec(4,4)
        ax_xy = fig.add_subplot(gs[1:4,0:3])
        ax_x = fig.add_subplot(gs[0,0:3])
        ax_y = fig.add_subplot(gs[1:4,3])

        # Plot outer product as an image.
        x_l, x_u = self.xvals[0], self.xvals[-1]
        y_l, y_u = self.yvals[0], self.yvals[-1]
        # Since the functions are evaluated on grid points, but the image
        # uses pixels, shift the extent end points by 1/2 px to put the
        # function values at pixel centers.
        x_l -= 0.5*(self.xvals[1]-self.xvals[0])
        x_u += 0.5*(self.xvals[-1]-self.xvals[-2])
        y_l -= 0.5*(self.yvals[1]-self.yvals[0])
        y_u += 0.5*(self.yvals[-1]-self.yvals[-2])
        # Note the transposition, as the row index runs along the vertical
        # dimension of mpl images.
        ax_xy.imshow(self.outer.transpose(), interpolation='bilinear', aspect='auto',
                     cmap=cc.cm.bwy_r, origin='lower', extent=(x_l, x_u, y_l, y_u))
        ax_xy.set_xlim(x_l, x_u)
        ax_xy.set_ylim(y_l, y_u)
        ax_xy.set_xlabel(r'$\lambda$ ($\AA$)')
        ax_xy.set_ylabel("Spot rotational phase")

        ax_x.set_xlim(x_l, x_u)
        ax_y.set_ylim(y_l, y_u)

        # x and y components:
        ax_x.plot(self.xvals, self.uvals)
        ax_x.set_ylabel('Ampl.')
        ax_y.plot(self.vvals, self.yvals)
        ax_y.set_xlabel('Ampl.')

        # Hide tick labels on components in margins.
        setp(ax_x.get_xticklabels(), visible=False)
        setp(ax_y.get_yticklabels(), visible=False)

        if lines:
            for line in lines:
                ax_xy.axvline(line, c='r', lw=.75, alpha=.5)

        if zoom:
            #x_l, x_u = zoom[0]
            #y_l, y_u = zoom[1]
            zfig = figure(figsize=(7.5,7*.75))
            imshow(self.outer.transpose(), interpolation='bilinear', aspect='auto',
                   cmap=cc.cm.bwy_r, origin='lower', extent=(x_l, x_u, y_l, y_u))
            print(x_l, x_u, y_l, y_u)
            print(zoom[0], zoom[1])
            xlim(*zoom[0])
            ylim(*zoom[1])

        return fig
