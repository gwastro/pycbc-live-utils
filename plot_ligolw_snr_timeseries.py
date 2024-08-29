#!/usr/bin/env python3

"""Plot triggers and SNR time series from a CBC event in LIGOLW format."""

import sys
import os
from ligo.skymap.io import LigoLWEventSource
import numpy as np
import pylab as pl
from pycbc.results import ifo_color


coinc_file = sys.argv[1]
output_base = sys.argv[2]

event, = LigoLWEventSource(coinc_file, psd_file=coinc_file, coinc_def=None).values()

end_time = np.mean([single.time for single in event.singles])

pl.figure(figsize=(10,8))

for single in event.singles:
    series = single.snr_series
    print(series.deltaT)
    series_time = np.arange(len(series.data.data)) * series.deltaT + float(series.epoch)
    series_abs = np.abs(series.data.data)
    pl.plot(series_time - end_time, series_abs, '-',
            color=ifo_color(single.detector), label=single.detector, zorder=2)
    if single.snr is not None:
        pl.plot(single.time - end_time, single.snr, 'o', color='w',
                markeredgecolor=ifo_color(single.detector), zorder=1)

pl.grid(ls='solid', alpha=0.25)
pl.xlabel('GPS time from {:.3f} [s]'.format(end_time))
pl.ylabel('SNR')
pl.legend(loc='upper left')
pl.title(os.path.basename(coinc_file))
pl.savefig(output_base + '{0:.3f}.png'.format(end_time), dpi=150)
