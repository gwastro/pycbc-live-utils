#!/usr/bin/env python
'''
Code to fetch data and split it into short duration frames for pycbc_live testing.
python split_frames.py \
        --channel-name H1:DCS-CALIB_STRAIN_C02 \
        --other-channel-names H1:DCS-CALIB_STATE_VECTOR_C02  \
        --gps-start-time 1186987282 --gps-end-time 1187012482 \
        --frame-type H1_HOFT_C00 --frame-duration 4 --outdir H1/
'''
import os
import logging
import argparse
import tqdm
from pycbc import frame
import pycbc.strain
from pycbc.types import float32, float64

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--other-channel-names",
                        type=str, nargs="+",
                        help="list of channels other than strain channel")
parser.add_argument("--gps-start-time",
                        help="The gps start time of the data "
                                "(integer seconds)", type=int)
parser.add_argument("--gps-end-time",
                        help="The gps end time of the data "
                                " (integer seconds)", type=int)
parser.add_argument('--frame-duration', type=int, default=1,
                    help='Split all data into smaller frame files of the given duration if specified.')
parser.add_argument('--outdir', type=str, help='prefix to save strain files', default='.')
parser.add_argument('--output-precision', type=str,
                    choices=['single', 'double'], default='double',
                    help='Precision of output strain, %(default)s by default')
parser.add_argument('--dyn-range-factor', action='store_true',
                    help='Scale the output strain by a large factor (%f) '
                         'to avoid underflows in subsequent '
                         'calculations' % pycbc.DYN_RANGE_FAC)
parser.add_argument('--low-frequency-cutoff', type=float,
                    help='Provide a low-frequency-cutoff for fake strain. '
                         'This is only needed if fake-strain or '
                         'fake-strain-from-file is used')
pycbc.strain.insert_strain_option_group(parser, gps_times=False)
## Strain channel is included as "--channel-name" in strain option
## data_reading_group.add_argument("--channel-name", type=str,
##                help="The channel containing the gravitational strain data")
## parser.add_argument("--frame-type",
##                         type=str,
##                         help="Use datafind to get the needed frame file(s) of this type.")

args = parser.parse_args()

logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)
start = args.gps_start_time
stop = args.gps_end_time
step = args.frame_duration

filename = '{0}-{1}-{2}-{3}.gwf'

if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

# read and condition strain as pycbc_inspiral would do
out_strain = pycbc.strain.from_cli(args, dyn_range_fac=pycbc.DYN_RANGE_FAC)

# force strain precision to be as requested
out_strain = out_strain.astype(
        float32 if args.output_precision == 'single' else float64)

# unless asked otherwise, revert the dynamic range factor
if not args.dyn_range_factor:
    out_strain /= pycbc.DYN_RANGE_FAC

# Adding strain channel info
data = {}
det_channels = {}
logging.info("Adding strain channel named {} ...".format(args.channel_name))
det_channels[args.channel_name.split(':')[0]] = [args.channel_name]
data[args.channel_name.split(':')[0]] = [out_strain]

if args.other_channel_names:
    for channel in args.other_channel_names:
        logging.info("Reading channel {} ...".format(channel))
        det = channel.split(':')[0]
        det_channels[det].append(channel)
        data[det].append(frame.query_and_read_frame(args.frame_type, channel, start, stop))

logging.info("Writing frames, each of duration {} sec".format(step))
for d in data.keys():
    for s in tqdm.trange(start, stop, step):
        fname = os.path.join(args.outdir,
                filename.format(d, args.frame_type, s, step if s+step < stop else stop-s))
        frame.write_frame(fname, det_channels[d], [ts.time_slice(s, s+step if s+step < stop else stop) for ts in data[d]])
