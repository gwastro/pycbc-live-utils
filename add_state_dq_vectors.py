#!/usr/bin/env python

"""Read a frame file containing strain data and add simulated state and DQ
vector channels to it, then write everything to a new frame file."""

import argparse
import numpy as np
from numpy.random import default_rng
from pycbc.types import TimeSeries
from pycbc.frame import read_frame, write_frame


parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument('--input-file', type=str, required=True)
parser.add_argument('--strain-channel', type=str, required=True)
parser.add_argument('--output-file', type=str, required=True)

parser.add_argument('--state-vector', type=str,
                    help='Name of state vector channel')
parser.add_argument('--state-vector-good', type=int,
                    help='Value of state vector indicating good data')
parser.add_argument('--state-off-segments', type=str, nargs='+',
                    metavar='START,STOP',
                    help='Segment(s) to be given an off state')

parser.add_argument('--dq-vector', type=str,
                    help='Name of DQ vector channel')
parser.add_argument('--dq-vector-good', type=int,
                    help='Value of DQ vector indicating good data')
parser.add_argument('--dq-bad-times', type=float, nargs='+',
                    metavar='TIME', help='Center time(s) of bad DQ epoch(s)')
parser.add_argument('--dq-bad-pad', type=float,
                    help='Duration of bad DQ epoch(s)')

parser.add_argument('--idq-channel', type=str, 
                    help='Name of idq channel')
parser.add_argument('--random-seed', type=int,
                    help='Random seed used to generate fake idq data')
parser.add_argument('--idq-bad-times', type=float, nargs='+',
                    help='Center time(s) of peaks in iDQ timeseries')
parser.add_argument('--idq-bad-pad', type=float,
                    help='Duration of peaks in iDQ data')

args = parser.parse_args()

# load frame file

strain = read_frame(args.input_file, args.strain_channel)

out_channel_names = [args.strain_channel]
out_timeseries = [strain]

# add state vector

if args.state_vector is not None:
    state_dt = 1. / 16.
    state_size = int(strain.duration / state_dt)
    state_data = np.zeros(state_size, dtype=np.uint32)
    state_data[:] = args.state_vector_good

    state_ts = TimeSeries(state_data, delta_t=state_dt,
                          epoch=strain.start_time)
    state_ts_times = state_ts.sample_times.numpy()

    for ss in args.state_off_segments:
        start, end = map(float, ss.split(','))
        fnz = np.flatnonzero(np.logical_and(state_ts_times >= start,
                                            state_ts_times < end))
        if len(fnz):
            state_ts[fnz[0]:fnz[-1]+1] = 0

    out_channel_names.append(args.state_vector)
    out_timeseries.append(state_ts)

# add DQ vector

if args.dq_vector is not None:
    # generate a fake DQ vector with random occasional vetoes
    dq_dt = 1. / 64.
    dq_size = int(strain.duration / dq_dt)
    dq_data = np.zeros(dq_size, dtype=np.uint32)
    dq_data[:] = args.dq_vector_good

    if args.dq_vector_good == 0:
        # Virgo DQ stream style
        dq_bad = 1
    else:
        # LIGO DQ vector style
        dq_bad = 0

    dq_ts = TimeSeries(dq_data, delta_t=dq_dt,
                       epoch=strain.start_time)
    dt_ts_times = dq_ts.sample_times.numpy()

    for dqt in args.dq_bad_times:
        fnz = np.flatnonzero(abs(dt_ts_times - dqt) < args.dq_bad_pad)
        if len(fnz):
            dq_ts[fnz[0]:fnz[-1]+1] = dq_bad

    out_channel_names.append(args.dq_vector)
    out_timeseries.append(dq_ts)

# add iDQ data

if args.idq_channel is not None:
    #generate a fake idq timeseries
    idq_dt = 1. / 128.
    idq_size = int(strain.duration / idq_dt)
    rng = default_rng(args.random_seed)
    idq_data = rng.standard_normal(idq_size)-1
    
    idq_ts = TimeSeries(idq_data, delta_t=idq_dt,
                        epoch = strain.start_time)
    idq_ts_times = idq_ts.sample_times.numpy()
    
    for idqt in args.idq_bad_times:
        fnz = np.flatnonzero(abs(idq_ts_times - idqt) < args.idq_bad_pad)
        if len(fnz):
            idq_ts[fnz[0]:fnz[-1]+1] += 6
    
    out_channel_names.append(args.idq_channel)
    out_timeseries.append(idq_ts)

# write frame file

write_frame(args.output_file, out_channel_names, out_timeseries)
