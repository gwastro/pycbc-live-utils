#!/usr/bin/env python

"""Plot a GW skymap in HEALPix FITS format and show the corresponding timing
information alongside."""

import argparse
import itertools
import numpy as np
import matplotlib.pyplot as pp
from scipy.optimize import brentq
from pycbc.distributions import HealpixSky
from pycbc.detector import Detector


cli_parser = argparse.ArgumentParser(description=__doc__)
cli_parser.add_argument(
    '--input-file', required=True
)
cli_parser.add_argument(
    '--output-file', required=True
)
# FIXME get the time and detectors from the FITS header
cli_parser.add_argument(
    '--event-gps-time', required=True, type=float
)
cli_parser.add_argument(
    '--detectors', nargs='+', default=['H', 'L', 'V']
)
cli_args = cli_parser.parse_args()

print('Sampling skymap')

distr = HealpixSky(healpix_file=cli_args.input_file)
samples = distr.rvs(100000)

time = cli_args.event_gps_time
detector_names = cli_args.detectors

detectors = [Detector(name + '1') for name in detector_names]

detector_pairs = list(itertools.combinations(detector_names, 2))

print('Calculating delays')

delays = []
median_delays = []
for d1, d2 in itertools.combinations(detectors, 2):
    delays.append(
        d1.time_delay_from_detector(
            d2, samples['ra'], samples['dec'], time
        )
    )
    median_delays.append(np.median(delays[-1]))

print('Calculating rings')

ring_ras = np.random.uniform(0, 2 * np.pi, 10000)
ring_decs = np.random.uniform(-np.pi / 2, np.pi / 2, 10000)

rings = []
for i, (d1, d2) in enumerate(itertools.combinations(detectors, 2)):
    rings.append(
        d1.time_delay_from_detector(d2, ring_ras, ring_decs, time) - median_delays[i]
    )

print('Plotting')

pp.figure(figsize=(10,10))

pp.subplot(2, 2, 1)
pp.hexbin(samples['ra'], samples['dec'], cmap='plasma_r', mincnt=1, lw=0)
for names, ring in zip(detector_pairs, rings):
    contours = pp.tricontour(
        ring_ras,
        ring_decs,
        ring,
        levels=[0],
        colors='k',
        linewidths=0.25
    )
    pp.clabel(
        contours,
        inline=True,
        fontsize=8,
        fmt=''.join(names)
    )
pp.xlabel('Right ascension [rad]')
pp.ylabel('Declination [rad]')

for i, (n1, n2) in enumerate(detector_pairs):
    pp.subplot(2, 2, 2 + i)
    j = i + 1
    if j >= len(detector_pairs):
        j = 0
    pp.hexbin(
        delays[i] * 1000,
        delays[j] * 1000,
        cmap='plasma_r',
        mincnt=1,
        lw=0
    )
    pp.grid()
    pp.xlabel(f'{"".join(detector_pairs[i])} time delay [ms]')
    pp.ylabel(f'{"".join(detector_pairs[j])} time delay [ms]')

pp.suptitle(sys.argv[1])
pp.tight_layout()
pp.savefig(cli_args.output_file, dpi=150)
