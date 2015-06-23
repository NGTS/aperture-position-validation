#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import
import argparse
import logging
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import wcs
import matplotlib.animation as manimation

logging.basicConfig(level='INFO', format='%(levelname)7s %(message)s')
logger = logging.getLogger(__name__)

half_width = 20

def get_data(filename, ra, dec):
    header = fits.getheader(filename)
    w = wcs.WCS(header)
    x, y = w.all_world2pix(ra, dec, 0)
    x_offset, y_offset = x % 1, y % 1
    xlim = (x - half_width, x + half_width)
    ylim = (y - half_width, y + half_width)
    vignette = np.log10(fits.getdata(filename)[ylim[0]:ylim[1], xlim[0]:xlim[1]] + 1)
    return vignette, x_offset, y_offset

def fetch_ra_dec_from_catalogue(catalogue, index):
    with fits.open(catalogue) as infile:
        if 'catalogue' in infile:
            cat = infile['catalogue'].data
            return (cat['ra'][index], cat['dec'][index])
        else:
            cat = infile[1].data
            return (
                    np.degrees(cat['ra'][index]),
                    np.degrees(cat['dec'][index]))


def main(args):
    if args.verbose:
        logger.setLevel('DEBUG')
    logger.debug(args)

    if not args.ra and not args.dec:
        ra, dec = fetch_ra_dec_from_catalogue(args.catalogue, args.index)
    else:
        ra, dec = args.ra, args.dec

    files = sorted(args.filename)
    nfiles = len(files)

    fig = plt.figure()
    initial_data, initial_x_offset, initial_y_offset = get_data(files[0], ra, dec)

    im = plt.imshow(initial_data, interpolation='None')
    ax = plt.gca()
    circle = plt.Circle(
            (half_width, + initial_x_offset, half_width + initial_y_offset),
            3, edgecolor='g', facecolor='None')

    ax.add_patch(circle)


    def update(i):
        fname = files[i]
        vignette, offset_x, offset_y = get_data(fname, ra, dec)
        im.set_data(vignette)
        im.set_clim(vignette.min(), vignette.max())
        circle.center = (half_width + offset_x, half_width + offset_y)

        print('\r{}/{}'.format(i, nfiles), end='')
        return im, circle

    ani = manimation.FuncAnimation(fig, update, frames=nfiles, interval=15, blit=True)
    if args.output is not None:
        ani.save(args.output, writer='ffmpeg')
    else:
        plt.show()





if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', nargs='+')
    parser.add_argument('-r', '--ra', type=float)
    parser.add_argument('-d', '--dec', type=float)
    parser.add_argument('-i', '--index', required=False, type=int)
    parser.add_argument('-c', '--catalogue', required=False)
    parser.add_argument('-o', '--output', required=False)
    parser.add_argument('-v', '--verbose', action='store_true')
    main(parser.parse_args())
