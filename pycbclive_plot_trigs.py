#!/usr/bin/env python

"""A tool for detailed inspection of single-detector triggers from PyCBC Live."""

import argparse
import os
import numpy as np
import pylab as pl
import h5py
from pycbc.results import ifo_color


class Autorange:
    """Utility to keep track of the total range
    spanned by values in a set of arrays.
    """
    def __init__(self):
        self.low = np.inf
        self.high = -np.inf

    def update(self, values):
        minv = min(values)
        maxv = max(values)
        if minv < self.low:
            self.low = minv
        if maxv > self.high:
            self.high = maxv

def plot_gate(ax, gate):
    g_time, g_width, g_taper = gate
    ax.axvspan(g_time - g_width, g_time + g_width,
               hatch='/', facecolor='none', edgecolor='#00ff00')


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--trigger-files', type=str, nargs='+', required=True, metavar='PATH',
                    help='List of HDF5 files containing triggers')
parser.add_argument('--highlight-times', type=float, nargs='+', metavar='GPS',
                    help='List of GPS times to mark with green dashed lines')
parser.add_argument('--gates', type=str, nargs='+', metavar='IFO,CENTER,WIDTH,TAPER',
                    help='List of gating parameters to display')
parser.add_argument('--output-plot', type=str, required=True,
                    help='Path to output plot')
args = parser.parse_args()

detectors = 'H1 L1 V1'.split()

fig = pl.figure(figsize=(20, 10))
ax = {}
n = 20
for i, detector in enumerate(detectors):
    ax[detector] = pl.subplot(len(detectors), n, (n*i + 1, n*i + n-1))
    ax[detector].grid()
    ax[detector].set_ylabel('{}\nTemplate duration [s]'.format(detector))
    if i < len(detectors) - 1:
        ax[detector].set_xticklabels([])
ax['cb'] = pl.subplot(1, n, n)

ar_time = Autorange()
ar_dur = Autorange()

for fn in sorted(args.trigger_files):
    print(fn)
    with h5py.File(fn, 'r') as trigfile:
        start_time = float(os.path.basename(fn).split('-')[-2])

        for detector in detectors:
            # show boundaries of different files
            ax[detector].axvline(start_time, ls=':', color='magenta')

            if detector not in trigfile:
                continue

            grp = trigfile[detector]

            # show gates
            if 'gates' in grp:
                for gate in grp['gates'][:]:
                    plot_gate(ax[detector], gate)

            # show triggers
            if 'end_time' not in grp or len(grp['end_time']) == 0:
                continue
            ar_time.update(grp['end_time'])
            ar_dur.update(grp['template_duration'])
            sorter = np.argsort(grp['snr'][:])
            sc = ax[detector].scatter(grp['end_time'][:][sorter], grp['template_duration'][:][sorter],
                                      c=grp['snr'][:][sorter], cmap='plasma_r', vmin=4.5, vmax=10)

ax[detectors[-1]].set_xlabel('GPS time')

for detector in detectors:
    ax[detector].set_xlim(ar_time.low, ar_time.high)
    ax[detector].set_ylim(ar_dur.low, ar_dur.high)
    ax[detector].set_yscale('log')
    for ht in (args.highlight_times or []):
        ax[detector].axvline(ht, ls='--', color='green')
    for g in (args.gates or []):
        gate = g.split(',')
        if gate[0] != detector:
            continue
        plot_gate(ax[detector], map(float, gate[1:]))

cb = fig.colorbar(sc, cax=ax['cb'], extend='both')
cb.set_label('SNR')

fig.tight_layout()
fig.savefig(args.output_plot, dpi=200)
