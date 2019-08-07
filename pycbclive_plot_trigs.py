#!/usr/bin/env python

"""A tool for detailed inspection of single-detector triggers from PyCBC Live."""

import argparse
import numpy as np
import pylab as pl
import h5py
from pycbc.results import ifo_color


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--trigger-files', type=str, nargs='+', required=True)
parser.add_argument('--highlight-times', type=float, nargs='+')
parser.add_argument('--output-plot', type=str, required=True)
args = parser.parse_args()

detectors = 'H1 L1 V1'.split()

pl.figure(figsize=(20, 10))
ax = {}
for i, detector in enumerate(detectors):
    ax[detector] = pl.subplot(len(detectors), 1, i + 1)
    ax[detector].set_title(detector + ', color = SNR (orange <= 4.5, blue >= 10)')
    ax[detector].grid()
    ax[detector].set_ylabel('Template duration [s]')

min_time = np.inf
max_time = -np.inf

for fn in sorted(args.trigger_files):
    print(fn)
    with h5py.File(fn, 'r') as trigfile:
        for detector in detectors:
            if detector not in trigfile:
                continue
            grp = trigfile[detector]
            if 'end_time' not in grp or len(grp['end_time']) == 0:
                continue
            if min(grp['end_time']) < min_time:
                min_time = min(grp['end_time'])
            if max(grp['end_time']) > max_time:
                max_time = max(grp['end_time'])
            sorter = np.argsort(grp['snr'][:])
            ax[detector].scatter(grp['end_time'][:][sorter], grp['template_duration'][:][sorter],
                                 c=grp['snr'][:][sorter], cmap='plasma_r', vmin=4.5, vmax=10)

ax[detectors[-1]].set_xlabel('GPS time')

for detector in detectors:
    ax[detector].set_xlim(min_time, max_time)
    ax[detector].set_yscale('log')
    for ht in args.highlight_times:
        ax[detector].axvline(ht, ls='--', color='green')

pl.tight_layout()
pl.savefig(args.output_plot, dpi=200)
