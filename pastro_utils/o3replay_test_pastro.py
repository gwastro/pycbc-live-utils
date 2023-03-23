#!/usr/bin/env python
import argparse
from ligo.gracedb.rest import GraceDb
from matplotlib import use
use('agg')
from matplotlib import pyplot as pl
import numpy
import h5py
from tqdm import tqdm
from pycbc import conversions

parser = argparse.ArgumentParser()
parser.add_argument('--gps-start-time', help='Start from this time')
parser.add_argument('--gps-end-time', help='Go up to this time')
args = parser.parse_args()


try:
    prev = h5py.File('o3replay_plive_pastros.hdf', 'r')
    gids = prev['gid'][:]
    snrs = prev['snr'][:]
    mchirps = prev['mchirp'][:]
    fars = prev['far'][:]
    pastros = prev['pastro'][:]
except:
    g = GraceDb(service_url='https://gracedb-playground.ligo.org/api/')
    gids = []
    gps = []
    fars = []
    mchirps = []
    snrs = []
    pastros = []
    query = (f'group: Test pipeline: pycbc submitter: "PyCBC Live" '
             f'{args.gps_start_time} .. {args.gps_end_time}')
    print(query)
    for ev in tqdm(g.events(query)):
        gids.append(ev['graceid'])
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

        gps.append(ev['gpstime'])
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

        # get the file name
        pastro_file = [f for f in g.files(ev['graceid']).json() if f.endswith('p_astro.json')]
        assert len(pastro_file) < 2
        if len(pastro_file):
            pastro_dict = g.files(ev['graceid'], filename=pastro_file[0]).json()
            pastros.append(1. - float(pastro_dict['p_terr']))
        else:
            # Missing file :/
            pastros.append(-1.)

        # Alert if large ...
        if pastros[-1] > 0.1:
            print('Interesting!', ev['graceid'], ev['gpstime'], pastros[-1])

    print('Got', str(len(fars)), 'events')

    with h5py.File('o3replay_plive_pastros.hdf', 'w') as f:
        f.create_dataset('gid', data=numpy.array(gids, dtype=h5py.string_dtype(encoding='utf-8')))
        f.create_dataset('gpstime', data=numpy.array(gps))
        f.create_dataset('far', data=numpy.array(fars))
        f.create_dataset('mchirp', data=numpy.array(mchirps))
        f.create_dataset('snr', data=numpy.array(snrs))
        f.create_dataset('pastro', data=numpy.array(pastros))

interesting = pastros > 0.1
print('pastro, GID, FAR, mchirp')
print(pastros[interesting], gids[interesting], fars[interesting], mchirps[interesting])

