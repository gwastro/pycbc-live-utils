#!/usr/bin/env python

import sys
import numpy as np
import h5py
import pylab as pl
import pycbc
from pycbc.results import ifo_color


pl.figure(figsize=(10,8))

with h5py.File(sys.argv[1], 'r') as hf:
    for ifo in ['H1', 'L1', 'V1']:
        if ifo + '/psd' not in hf:
            continue
        df = hf[ifo + '/psd'].attrs['delta_f']
        asd = hf[ifo + '/psd'][:] ** 0.5 / pycbc.DYN_RANGE_FAC
        f = np.arange(len(asd)) * df
        pl.loglog(f, asd, '-', label=ifo, color=ifo_color(ifo))

pl.title(sys.argv[1])
pl.xlim(10, 1024)
pl.ylim(1e-24, 1e-19)
pl.xlabel('Frequency [Hz]')
pl.ylabel('Amplitude spectral density')
pl.grid()
pl.legend()
pl.tight_layout()
pl.savefig(sys.argv[2])
