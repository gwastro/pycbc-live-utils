#!/usr/bin/env python

"""Match triggers with a list of injections and write the result as a
inspinjfind database that can be used with ligo-skymap-stats."""

import os
import argparse
import glob
import sqlite3 as sql
import numpy as np
import tqdm
import glue.ligolw.utils
import glue.ligolw.table
import glue.ligolw.ligolw
import glue.ligolw.lsctables


def coinc_index_to_id(i):
    return 'coinc_event:coinc_event_id:{}'.format(i)


class LIGOLWContentHandler(glue.ligolw.ligolw.LIGOLWContentHandler):
    pass

glue.ligolw.lsctables.use_in(LIGOLWContentHandler)

# parse args

parser = argparse.ArgumentParser()
parser.add_argument('--trig-glob', type=str, required=True)
parser.add_argument('--inj-file', type=str, required=True)
parser.add_argument('--time-tolerance', type=float, default=1)
parser.add_argument('--mchirp-tolerance', type=float, default=0.5)
parser.add_argument('--output-file', type=str, required=True)
args = parser.parse_args()

if os.path.exists(args.output_file):
    parser.error('Output database exists and I will not overwrite it')

# load injections

xmldoc = glue.ligolw.utils.load_filename(
        args.inj_file, verbose=False, contenthandler=LIGOLWContentHandler)
sim_inspiral = glue.ligolw.lsctables.SimInspiralTable.get_table(xmldoc)

# load triggers

if '*' in args.trig_glob:
    coincs = []
    for fn in tqdm.tqdm(glob.glob(args.trig_glob)):
        xmldoc = glue.ligolw.utils.load_filename(
                fn, verbose=False, contenthandler=LIGOLWContentHandler)
        coinc = glue.ligolw.lsctables.CoincInspiralTable.get_table(xmldoc)[0]
        coincs.append(coinc)
else:
    xmldoc = glue.ligolw.utils.load_filename(
            args.trig_glob, verbose=False, contenthandler=LIGOLWContentHandler)
    coincs = glue.ligolw.lsctables.CoincInspiralTable.get_table(xmldoc)

# do the matching

sim_times = np.array([si.geocent_end_time + si.geocent_end_time_ns * 1e-9 for si in sim_inspiral])
sim_mchirps = np.array([si.mchirp for si in sim_inspiral])
coinc_times = np.array([ci.end_time + ci.end_time_ns * 1e-9 for ci in coincs])
coinc_mchirps = np.array([ci.mchirp for ci in coincs])

matches = []
for i, si in enumerate(sim_inspiral):
    delta_t = abs(sim_times[i] - coinc_times)
    time_mask = delta_t < args.time_tolerance

    delta_mchirp = abs(sim_mchirps[i] - coinc_mchirps) / sim_mchirps[i]
    mchirp_mask = delta_mchirp < args.mchirp_tolerance

    found_idx = np.flatnonzero(np.logical_and(time_mask, mchirp_mask))
    if len(found_idx) == 1:
        matches.append((si.simulation_id, coincs[found_idx[0]].coinc_event_id))
    elif len(found_idx) > 1:
        raise RuntimeError('Simulation {} matches trigs {}'.format(si.simulation_id, ', '.join(map(str, found_idx))))

# write output
# here I just write the minimum amount of stuff
# that will make ligo-skymap-stats work

with sql.connect(args.output_file) as odb:
    # insert the simulations

    query = ('CREATE TABLE sim_inspiral '
             '(latitude REAL, longitude REAL,'
             'distance REAL, simulation_id VARCHAR)')
    odb.execute(query)

    for si in sim_inspiral:
        query = 'INSERT INTO sim_inspiral(latitude, longitude, distance, simulation_id) VALUES (?, ?, ?, ?)'
        odb.execute(query, (si.latitude, si.longitude, si.distance, si.simulation_id))

    # insert the triggers

    query = ('CREATE TABLE coinc_inspiral '
             '(combined_far REAL, snr REAL,'
             'coinc_event_id VARCHAR)')
    odb.execute(query)

    for ci in coincs:
        query = 'INSERT INTO coinc_inspiral(combined_far, snr, coinc_event_id) VALUES (?, ?, ?)'
        odb.execute(query, (ci.combined_far, ci.snr, ci.coinc_event_id))

    # insert the matches

    query = ('CREATE TABLE coinc_event_map (coinc_event_id VARCHAR,'
             'event_id VARCHAR, table_name VARCHAR)')
    odb.execute(query)

    for sim_id, coinc_id in matches:
        query = 'INSERT INTO coinc_event_map(coinc_event_id, event_id, table_name) VALUES (?, ?, ?)'
        odb.execute(query, (coinc_id, sim_id, 'sim_inspiral'))
        odb.execute(query, (coinc_id, coinc_id, 'coinc_event'))
