#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import
import argparse
import logging
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from astropy.io import fits
from astropy import wcs

logging.basicConfig(level='INFO', format='%(levelname)7s %(message)s')
logger = logging.getLogger(__name__)


def get_data(filename, ra, dec, half_width):
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
            ra, dec = cat['ra'][index], cat['dec'][index]
        else:
            cat = infile[1].data
            ra, dec = (np.degrees(cat['ra'][index]),
                    np.degrees(cat['dec'][index]))
    return ra, dec


def main(args):
    if args.verbose:
        logger.setLevel('DEBUG')
    logger.debug(args)

    if not args.ra and not args.dec:
        ra, dec = fetch_ra_dec_from_catalogue(args.catalogue, args.index)
    else:
        ra, dec = args.ra, args.dec

    logger.info('RA: %s, DEC: %s', ra, dec)

    files = sorted(args.filename)
    nfiles = len(files)
    half_width = args.half_width

    fig = plt.figure()
    initial_data, initial_x_offset, initial_y_offset = get_data(files[0], ra, dec,
                                                                half_width)

    im = plt.imshow(initial_data, interpolation='None',
            origin='lower')
    ax = plt.gca()
    circle = plt.Circle(
        (half_width, +initial_x_offset, half_width + initial_y_offset), 3,
        edgecolor='g',
        facecolor='None')

    ax.add_patch(circle)

    def update(i):
        fname = files[i]
        vignette, offset_x, offset_y = get_data(fname, ra, dec, half_width)
        im.set_data(vignette)
        med_vignette = np.average(vignette)
        im.set_clim(0.99 * med_vignette, 1.2 * med_vignette)
        circle.center = (half_width + offset_x, half_width + offset_y)

        print('\r{}/{}'.format(i + 1, nfiles), end='')
        return im, circle

    interval = 1000. / args.fps

    ani = manimation.FuncAnimation(fig, update, frames=nfiles, interval=interval, blit=True)
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
    parser.add_argument('-w', '--half-width', required=False, type=int, default=20)
    parser.add_argument('--fps', required=False, default=24, type=int)
    main(parser.parse_args())
