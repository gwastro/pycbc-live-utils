#!/usr/bin/env python

"""A tool for detailed inspection of single-detector triggers from PyCBC Live."""

# TODO missing features that used to be in `pycbclive_plot_trigs`:
# * Plotting gates from HDF files.
# * Showing the boundaries of analysis segments
#   (probably only useful when plotting a few segments only).

import argparse
import logging
import os
import tqdm
import numpy as np
import matplotlib
matplotlib.use('agg')
import pylab as pl
import h5py
import lal
import glob
from ligo.segments import segment, segmentlist
from pycbc.events.ranking import newsnr_sgveto


class Autorange:
    """Utility to keep track of the total range
    spanned by values in a set of arrays.
    """
    def __init__(self):
        self.low = np.inf
        self.high = -np.inf
        self.set = False

    def update(self, values):
        minv = min(values)
        maxv = max(values)
        if minv < self.low:
            self.low = minv
        if maxv > self.high:
            self.high = maxv
        self.set = True


def plot_gate(ax, gate):
    g_time, g_width, g_taper = gate
    ax.axvspan(
        g_time - g_width,
        g_time + g_width,
        hatch='/',
        facecolor='none',
        edgecolor='#00ff00'
    )


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    '--trigger-files-glob',
    type=str,
    required=True,
    metavar='PATH',
    help='Glob command to find HDF5 files containing triggers.'
         'Must have the wildcards escaped, using a backslash, '
         'or surrounding the argument in quotes'
)
parser.add_argument(
    '--detectors',
    type=str,
    nargs='+',
    default=['H1', 'L1', 'V1'],
    help='Which detectors to plot'
)
parser.add_argument(
    '--highlight-times',
    type=float,
    nargs='+',
    metavar='GPS',
    help='List of GPS times to mark with green dashed lines'
)
parser.add_argument(
    '--gates',
    type=str,
    nargs='+',
    metavar='IFO,CENTER,WIDTH,TAPER',
    help='List of gating parameters to display'
)
parser.add_argument(
    '--output-plot',
    type=str,
    required=True,
    help='Path to output plot'
)
args = parser.parse_args()

detectors = args.detectors

fig = pl.figure(figsize=(20, 10))
ax = {}
n = 22
for i, detector in enumerate(detectors):
    ax[detector] = pl.subplot(len(detectors), n, (n*i + 1, n*i + n-1))
    ax[detector].set_ylabel(f'{detector}\nTemplate duration [s]')
    if i < len(detectors) - 1:
        ax[detector].set_xticklabels([])
ax['cb'] = pl.subplot(1, n, n)

logging.info('Reading triggers')

file_segs = segmentlist([])
trig_segs = {d: segmentlist([]) for d in detectors}
triggers = {d: None for d in detectors}
gates = {d: [] for d in detectors}

trigger_files = glob.glob(args.trigger_files_glob)

for fn in tqdm.tqdm(sorted(trigger_files)):
    try:
        with h5py.File(fn, 'r') as trigfile:
            fn_fields = os.path.basename(fn).replace('.hdf', '').split('-')
            file_start_time = int(float(fn_fields[-2]))
            file_end_time = int(float(fn_fields[-2]) + float(fn_fields[-1])) + 1
            file_segment = segment(file_start_time, file_end_time)
            file_segs.append(file_segment)

            for detector in detectors:
                if detector not in trigfile:
                    continue

                grp = trigfile[detector]

                if 'gates' in grp:
                    for gate in grp['gates'][:]:
                        gates[detector].append(gate)

                if 'end_time' not in grp or len(grp['end_time']) == 0:
                    continue

                trig_segs[detector].append(file_segment)

                trig_times = grp['end_time'][:]
                trig_durs = grp['template_duration'][:]
                trig_ranks = newsnr_sgveto(
                    grp['snr'][:],
                    grp['chisq'][:],
                    grp['sg_chisq'][:]
                )
                if triggers[detector] is None:
                    triggers[detector] = [trig_times, trig_durs, trig_ranks]
                else:
                    triggers[detector] = [
                        np.concatenate((triggers[detector][0], trig_times)),
                        np.concatenate((triggers[detector][1], trig_durs)),
                        np.concatenate((triggers[detector][2], trig_ranks))
                    ]
    except OSError:
        logging.error(f'Failed reading {fn}, ignoring')

logging.info('Plotting')

ar_dur = Autorange()
for detector in detectors:
    if triggers[detector] is not None:
        ar_dur.update(triggers[detector][1])

file_segs.coalesce()
for detector in detectors:
    trig_segs[detector].coalesce()

ax[detectors[-1]].set_xlabel('Time')

for detector in detectors:
    axd = ax[detector]
    if triggers[detector] is None:
        axd.text(
            0.5,
            0.5,
            'No triggers',
            horizontalalignment='center',
            verticalalignment='center',
            transform=axd.transAxes
        )
        axd.set_yticks([])
        continue
    axd.grid()
    axd.set_yscale('log')
    axd.set_ylim(ar_dur.low * 0.8, ar_dur.high * 1.2)
    # plot segments
    axd.hlines(
        [ar_dur.low * 0.8] * len(file_segs),
        [s[0] for s in file_segs],
        [s[1] for s in file_segs],
        color='#ff0000',
        lw=3
    )
    axd.hlines(
        [ar_dur.low * 0.8] * len(trig_segs[detector]),
        [s[0] for s in trig_segs[detector]],
        [s[1] for s in trig_segs[detector]],
        color='#00ff00',
        lw=3
    )
    # plot triggers
    sorter = np.argsort(triggers[detector][2])
    print('Max {} rw SNR {:.2f} at {:.3f}'.format(
        detector,
        triggers[detector][2][sorter[-1]],
        triggers[detector][0][sorter[-1]]
    ))
    axd.scatter(
        triggers[detector][0][sorter],
        triggers[detector][1][sorter],
        c=triggers[detector][2][sorter],
        vmin=6,
        vmax=12,
        cmap='magma_r',
        s=4,
        lw=0
    )
    for ht in (args.highlight_times or []):
        axd.axvline(ht, ls='--', color='green')
    # plot gates
    for g in (args.gates or []):
        gate = g.split(',')
        if gate[0] != detector:
            continue
        plot_gate(axd, map(float, gate[1:]))

# make nice time ticks
file_segs_extent = file_segs.extent()
min_gps = int(file_segs_extent[0])
max_gps = int(file_segs_extent[1])
min_utc = list(lal.GPSToUTC(min_gps))
if max_gps - min_gps > 3600:
    # ticks every hour
    min_utc[4] = 0
    time_tick_delta = 3600
    time_tick_fmt = '{0:04d}\n{1:02d}-{2:02d}\n{3:02d} UTC'
else:
    # ticks every minute
    time_tick_delta = 60
    time_tick_fmt = '{0:04d}-{1:02d}-{2:02d}\n{3:02d}:{4:02d} UTC'
min_utc[5] = 0
min_gps = lal.UTCToGPS(tuple(min_utc))
time_ticks = []
time_tick_labels = []
while True:
    if min_gps in file_segs_extent:
        time_ticks.append(min_gps)
        time_tick_labels.append(
            time_tick_fmt.format(*lal.GPSToUTC(min_gps))
        )
    min_gps += time_tick_delta
    if min_gps >= max_gps:
        break
for detector in detectors:
    ax[detector].set_xlim(
        file_segs_extent[0],
        file_segs_extent[1]
    )
    ax[detector].set_xticks(time_ticks)
ax[detectors[-1]].set_xticklabels(time_tick_labels)

# add colorbar
cb = fig.colorbar(
    matplotlib.cm.ScalarMappable(
        matplotlib.colors.Normalize(vmin=6, vmax=12),
        cmap='magma_r'
    ),
    cax=ax['cb'],
    extend='both'
)
cb.set_label('$\\chi^2$-weighted SNR')

fig.tight_layout()

logging.info('Saving plot')

os.makedirs(os.path.dirname(args.output_plot), exist_ok=True)
fig.savefig(args.output_plot, dpi=150)

logging.info('Done')
