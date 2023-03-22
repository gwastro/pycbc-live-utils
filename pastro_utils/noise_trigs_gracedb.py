#!/usr/bin/env python
from ligo.gracedb.rest import GraceDb
from matplotlib import use
use('agg')
from matplotlib import pyplot as pl
import numpy
import h5py
from tqdm import tqdm
from pycbc import conversions


try:
    prev = h5py.File('o3_plive_fargt1pwk_nosnrmax.hdf', 'r')
    snrs = prev['snr'][:]
    mchirps = prev['mchirp'][:]
    fars = prev['far'][:]
except:
    g = GraceDb()
# only get high FAR events, > 1/week as most likely noise
# expect 2885 returns based on 2022/10/14 web interface search (not excluding snr max duplicates)
    fars = []
    mchirps = []
    snrs = []
    query = ('pipeline: pycbc submitter: "PyCBC Live" far > 1.6e-6')
    for ev in tqdm(g.events(query)):
        # needed to check for SNR maximization followup events
        try:  # try the first sngl
            chisqdof = ev['extra_attributes']['SingleInspiral'][0]['chisq_dof']
        except:  # first sngl could be the followup ifo, try the second :/
            chisqdof = ev['extra_attributes']['SingleInspiral'][1]['chisq_dof']
        if chisqdof == 1:
            if ev['extra_attributes']['CoincInspiral']['snr'] < 6.5:
                print('Low maximized SNR!', ev['graceid'])
            #print('SNR max, skipping')
            continue

        fars.append(ev['far'])
        mchirp_coinc = ev['extra_attributes']['CoincInspiral']['mchirp']
        mchirp_sngl = ev['extra_attributes']['SingleInspiral'][0]['mchirp']
        # check for sngl vs coinc bug
        if numpy.log(mchirp_coinc/mchirp_sngl) > 0.001:
            raise ValueError("Sngl vs coinc mchirps don't match!", ev['graceid'])
        else:
            mchirps.append(mchirp_sngl)
        snrs.append(ev['extra_attributes']['CoincInspiral']['snr'])
        #if snrs[-1] > 12 and mchirps[-1] < 2.5:
        #    print('loud BNS!?', ev['graceid'])

    print('Got', str(len(fars)), 'events')

    with h5py.File('o3_plive_fargt1pwk_nosnrmax.hdf', 'w') as f:
        f.create_dataset('far', data=numpy.array(fars))
        f.create_dataset('mchirp', data=numpy.array(mchirps))
        f.create_dataset('snr', data=numpy.array(snrs))

pl.loglog(numpy.array(mchirps), numpy.array(snrs), "b+")
pl.grid(True)
pl.xlabel('mchirp')
pl.ylabel('network SNR')
pl.savefig('pycbc_live_fargt1pwk_mchirp_vs_snr_nomax.png')
pl.close()

mybins = numpy.array([0.87, 1.74, 3.48, 6.96, 13.92, 27.84, 55.68, 111.36, 500.0])

with h5py.File('/home/pycbc.live/production/O3/ER13_H1L1_OPT_FLOW_HYBRID_BANK.hdf', 'r') as o3bank:
    bank_mc = conversions.mchirp_from_mass1_mass2(o3bank['mass1'][:], o3bank['mass2'][:])

pl.hist(mchirps, bins=mybins, histtype='bar', rwidth=0.8, log=True, label='O3 events w/ FAR>1/wk')
pl.hist(bank_mc, bins=mybins, weights=numpy.ones_like(bank_mc) * len(mchirps)/len(bank_mc),
    histtype='bar', rwidth=0.6, log=True, label='O3 templates (scaled)')
pl.loglog()
pl.xlim(0.7, 500)
pl.ylim(ymin=0.025)
pl.legend()
pl.xlabel('mchirp')
pl.ylabel('Counts')
pl.savefig('pycbc_live_fargt1pwk_mchirp_hist_nomax.png')
pl.close()

print('Done!')
