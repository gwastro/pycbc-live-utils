#!/usr/bin/env python3

"""Make a plot showing PyCBC Live's lag and number of usable detectors as a
function of time and MPI rank, for a given UTC date.
"""

# Things that would be good to implement next:
# * Save parsed data to a file next to the plot, as the old lag plotter did
# * Improve timezone correction (or avoid it by changing PyCBC Live's logging)
# * Show min/avg/max of lag for rank>0 instead of each curve (Ian's suggestion)
# * Break the curves when PyCBC Live restarts
# * Hourly plots

import argparse
import logging
import os
import glob
import time
import datetime
import numpy as np
from astropy.time import Time

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as pp


def iso_to_gps(iso_time):
    """Convert a string (or list of strings) representing ISO time to the
    corresponding GPS time.

    Do not call this in a loop over many times, as it is somewhat slow.
    Instead, give it the list of times to convert.
    """
    return Time(iso_time, format='isot').gps

def set_up_x_axis(ax, day):
    """Configure the horizontal plot axis with hourly ticks,
    a range spanning the given day, and an appropriate label.
    """
    tick_locs = []
    tick_labels = []
    for hour in range(24):
        tick_locs.append(f'{day}T{hour:02d}:00:00Z')
        tick_labels.append(f'{hour:02d}')
    tick_locs = iso_to_gps(tick_locs)
    ax.set_xticks(tick_locs)
    ax.set_xticklabels(tick_labels)
    ax.set_xlim(tick_locs[0], tick_locs[-1] + 3600)  # FIXME leap seconds
    ax.set_xlabel('UTC hour of day')

def date_argument(date_str):
    if date_str == 'today':
        return datetime.datetime.utcnow().date()
    return datetime.date(*map(int, date_str.split('-')))

def parse_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--day',
        type=date_argument,
        default='today',
        metavar='{YYYY-MM-DD, today}',
        help='Which (UTC) day we want to show, default today.'
    )
    parser.add_argument(
        '--log-glob',
        required=True,
        help='Glob expression for finding all PyCBC Live log files.'
    )
    parser.add_argument(
        '--output-path',
        required=True,
        help='Path to a directory where plots will be saved.'
    )
    parser.add_argument(
        '--psd-inverse-length',
        type=float,
        default=3.5,
        help='Should match the setting used in the analysis. '
             'Determines the bottom of the lag axis.'
    )
    parser.add_argument(
        '--verbose',
        action='store_true'
    )
    return parser.parse_args()


args = parse_cli()

logging.basicConfig(
    format='%(asctime)s %(message)s',
    level=(logging.INFO if args.verbose else logging.WARN)
)

gps_now = Time.now().gps

# read data by parsing log files
# FIXME PyCBC Live's log timestamps do not specify the timezone (or UTC
# offset). This is a mess, and forces us to reconstruct them here in order to
# ultimately convert them to GPS.  For now, this is done by first assuming they
# are UTC times and converting to GPS, then assuming that the UTC offset of
# each timestamp matches the UTC offset of the local machine now, and
# subtracting that offset from the GPS timestamps.  This will be correct most
# of the time, but not always (e.g.  if there was a daylight saving switch
# meanwhile).  Maybe the pytz module could help here, but ideally we should
# just store the UTC offset in the log, or make the log always output UTC
# times.
prev_day_str = str(args.day - datetime.timedelta(days=1))
next_day_str = str(args.day + datetime.timedelta(days=1))
times = {}
data = {}
files = sorted(glob.glob(args.log_glob))
for f in files:
    logging.info('Parsing %s', f)
    with open(f, 'r') as log_f:
        for line in log_f:
            if 'Took' not in line and 'Starting' not in line:
                continue
            if line < prev_day_str or line > next_day_str:
                # skip lines from the wrong days
                continue
            fields = line.split()
            if len(fields) not in [5, 16]:
                continue
            log_date = fields[0]
            log_time = fields[1].replace(',', '.')
            log_rank = int(fields[3])
            if 'Starting' in line:
                # adding the nan's will tell matplotlib
                # to break the curves when PyCBC Live starts
                log_lag = np.nan
                log_n_det = np.nan
            else:
                log_lag = float(fields[11])
                log_n_det = int(fields[13])
            if log_rank not in times:
                times[log_rank] = []
            if log_rank not in data:
                data[log_rank] = []
            times[log_rank].append(f'{log_date}T{log_time}Z')
            data[log_rank].append((log_n_det, log_lag))

timezone_offset = time.localtime().tm_gmtoff
logging.info(
    'Converting timestamps to GPS; assuming fixed UTC offset %s s',
    timezone_offset
)
for rank in times:
    times[rank] = iso_to_gps(times[rank])
    times[rank] -= timezone_offset

num_procs = len(data)
logging.info('%d procs', num_procs)

logging.info('Plotting')
pp.figure(figsize=(15,7))
ax_lag = pp.subplot(2, 1, 1)
ax_n_det = pp.subplot(2, 1, 2)

for rank in sorted(data):
    data[rank] = np.array(data[rank])
    if rank == 0:
        # rank-0 gives the total lag, so make it stand out
        color = '#f00'
    else:
        color = pp.cm.viridis(rank / (num_procs - 1))
    sorter = np.argsort(times[rank])
    n_det = data[rank][sorter,0]
    lag = data[rank][sorter,1]
    if rank in [0, 1, num_procs // 2, num_procs - 1]:
        label = f'Rank {rank}'
    else:
        label = None
    # rank-0 gives the total lag, so put it on top
    zorder = num_procs - rank
    ax_lag.plot(
        times[rank][sorter],
        lag,
        '.-',
        lw=0.5,
        markersize=3,
        markeredgewidth=0,
        color=color,
        label=label,
        zorder=zorder
    )
    ax_n_det.plot(
        times[rank][sorter],
        n_det,
        '.-',
        lw=0.5,
        markersize=3,
        markeredgewidth=0,
        color=color,
        zorder=zorder
    )

pp.suptitle(args.day)
ax_lag.axvspan(
    gps_now,
    gps_now + 86400,
    edgecolor='none',
    facecolor='#d0d0d0'
)
set_up_x_axis(ax_lag, args.day)
ax_lag.set_ylabel('Lag [s]')
ax_lag.set_ylim(args.psd_inverse_length, 400)
ax_lag.set_yscale('log')
ax_lag.grid(which='both')
ax_lag.legend(loc='upper left')
ax_n_det.axvspan(
    gps_now,
    gps_now + 86400,
    edgecolor='none',
    facecolor='#d0d0d0'
)
set_up_x_axis(ax_n_det, args.day)
ax_n_det.set_ylabel('Number of usable detectors')
ax_n_det.set_yticks([0, 1, 2, 3])
ax_n_det.grid()

pp.tight_layout()

out_path = os.path.join(
    args.output_path,
    f'{args.day.year:04d}',
    f'{args.day.month:02d}',
    f'{args.day.day:02d}',
    f'{args.day}_lag_over_time.png'
)
logging.info('Saving plot to %s', out_path)
os.makedirs(os.path.dirname(out_path), exist_ok=True)
pp.savefig(out_path, dpi=200)

for hour in range(0, 24):
    pp.suptitle(f'{args.day}T{hour:02d}')
    tick_locs = []
    tick_labels = []
    for minute in range(60):
        tick_locs.append(f'{args.day}T{hour:02d}:{minute:02d}:00Z')
        tick_labels.append(f'{minute:02d}')
    tick_locs = iso_to_gps(tick_locs)
    if tick_locs[0] > gps_now:
        break
    for ax in [ax_lag, ax_n_det]:
        ax.set_xticks(tick_locs)
        ax.set_xticklabels(tick_labels)
        ax.set_xlim(tick_locs[0], tick_locs[-1] + 60)  # FIXME leap seconds
        ax.set_xlabel('Minute')
    pp.tight_layout()
    out_path = os.path.join(
        args.output_path,
        f'{args.day.year:04d}',
        f'{args.day.month:02d}',
        f'{args.day.day:02d}',
        f'{args.day}T{hour:02d}_lag_over_time.png'
    )
    logging.info('Saving plot to %s', out_path)
    pp.savefig(out_path, dpi=200)

logging.info('Done')
