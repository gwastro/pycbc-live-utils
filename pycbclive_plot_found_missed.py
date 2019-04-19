#!/usr/bin/env python

import argparse
import glob
import tqdm
import numpy as np
import pylab as pl
import lal
from matplotlib.colors import LogNorm
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import ligolw, table, lsctables


class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
    pass

lsctables.use_in(LIGOLWContentHandler)

parser = argparse.ArgumentParser()
parser.add_argument('--injection-file', type=str, required=True)
parser.add_argument('--trigger-glob', type=str, required=True)
parser.add_argument('--plot-file', type=str, required=True)
parser.add_argument('--sens-plot-file', type=str)
parser.add_argument('--x-axis', type=str, choices=['time', 'mchirp'],
                    required=True)
args = parser.parse_args()

# read injections

doc = ligolw_utils.load_filename(args.injection_file, False,
                                 contenthandler=LIGOLWContentHandler)
sim_table = table.get_table(doc, lsctables.SimInspiralTable.tableName)
injections = []
for sim in sim_table:
    injections.append((float(sim.get_time_geocent()),
                       sim.mchirp,
                       sim.alpha1,
                       sim.alpha2,
                       sim.alpha3))
del doc
print('{} injections'.format(len(injections)))

# read triggers

triggers = []
for file_path in tqdm.tqdm(glob.glob(args.trigger_glob)):
    doc = ligolw_utils.load_filename(file_path, False, contenthandler=LIGOLWContentHandler)
    coinc_table = table.get_table(doc, lsctables.CoincInspiralTable.tableName)
    end_time, far, mchirp = coinc_table[0].end_time + coinc_table[0].end_time_ns * 1e-9, coinc_table[0].combined_far, coinc_table[0].mchirp
    triggers.append([end_time, far, mchirp])
    del doc
triggers = np.array(triggers)
print('{} triggers'.format(triggers.shape[0]))

# associate

found = []
for injt, mchirp, osnr_h, osnr_l, osnr_v in injections:
    delta_t = abs(triggers[:,0] - injt)
    delta_mchirp = abs(triggers[:,2] - mchirp) / mchirp
    i = np.argmin(delta_t)
    far = triggers[i,1] if (delta_t[i] < 0.5 and delta_mchirp[i] < 0.5) else np.nan
    found.append([injt, mchirp, osnr_h, osnr_l, osnr_v, far])
found = np.array(found)

decisive_snr = np.array([sorted(found[i,2:5])[1] for i in range(found.shape[0])])
found_mask = np.isfinite(found[:,5])
missed_mask = np.logical_not(found_mask)

title = '{}/{} injections found'.format(sum(found_mask), len(found_mask))

if args.x_axis == 'time':
    x_quantity = found[:,0] - found[0,0]
    x_label = 'Injection time (origin at first injection) [s]'
    x_log = False
elif args.x_axis == 'mchirp':
    x_quantity = found[:,1]
    x_label = 'Injection chirp mass [$M_\\odot$]'
    x_log = True

pl.figure(figsize=(14,7))

pl.plot(x_quantity[missed_mask], decisive_snr[missed_mask], 'xr')
pl.scatter(x_quantity[found_mask], decisive_snr[found_mask],
           c=found[found_mask,5], norm=LogNorm(vmin=1e-9, vmax=1e-4),
           vmin=1e-9, vmax=1e-4)
pl.axhline(4.5, ls='--', color='k')
if x_log:
    pl.xscale('log')
pl.yscale('log')
pl.xlabel(x_label)
pl.ylabel('Decisive optimal SNR')
cb = pl.colorbar(extend='both')
cb.set_label('FAR [Hz]')
pl.title(title)

pl.tight_layout()
pl.savefig(args.plot_file)

if args.sens_plot_file is not None:
    pl.figure()

    ifar = np.sort(1. / found[found_mask,5])
    count = np.arange(len(ifar))[::-1] + 1

    pl.step(ifar / lal.YRJUL_SI, count)
    pl.xscale('log')
    pl.xlim(1e-4, 1e3)
    pl.yscale('log')
    pl.ylim(1, 1e3)
    pl.grid()
    pl.xlabel('Inverse FAR [yr]')
    pl.ylabel('Cumulative number of injections')

    pl.tight_layout()
    pl.savefig(args.sens_plot_file)
