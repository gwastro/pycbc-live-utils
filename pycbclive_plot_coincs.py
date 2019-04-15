#!/usr/bin/env python

"""Plot coincident triggers from PyCBC Live."""

import argparse
import glob
import h5py
import numpy as np
import matplotlib as mpl
mpl.use('agg')
import pylab as pl


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--trigger-glob', type=str, required=True)
parser.add_argument('--output-file', type=str, required=True)
args = parser.parse_args()

ifars = []
stats = []
ifos = []
times = []
tdurs = []

for fn in glob.glob(args.trigger_glob):
    with h5py.File(fn, 'r') as f:
        if 'foreground' not in f:
            continue
        stat = f['foreground/stat'][0]
        ifos = f['foreground/type'][()]
        ifos = sorted(ifos.split('-'))
        time = np.mean([f['foreground/' + ifo + '/end_time'][()] for ifo in ifos])
        tdur = f['foreground/' + ifos[0] + '/template_duration'][()]
        stats.append(stat)
        ifos.append(ifos)
        times.append(time)
        tdurs.append(tdur)

stats = np.array(stats)
times = np.array(times)
tdurs = np.array(tdurs)
sorter = np.argsort(stats)

pl.rcParams['font.size'] = 8

pl.scatter(times[sorter], tdurs[sorter], c=stats[sorter], s=3, cmap='viridis_r')
pl.yscale('log')
pl.xlabel('GPS time [s]')
pl.ylabel('Template duration [s]')
cb = pl.colorbar()
cb.set_label('Ranking statistic')
pl.grid()
pl.title('Coincident triggers from PyCBC Live\n' + args.trigger_glob)

pl.tight_layout()
pl.savefig(args.output_file, dpi=200)
