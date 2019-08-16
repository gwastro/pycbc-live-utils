
'''
Code to fetch data and split it into short duration frames for pycbc_live testing.
python split_frames.py \
        --channel-names H1:DCS-CALIB_STATE_VECTOR_C02 H1:DCS-CALIB_STRAIN_C02 \
        --gps-start-time 1186987282 --gps-end-time 1187012482 \
        --frame-type H1_HOFT_C02 --frame-duration 4 --outdir H1/

Then to run live offline:
mpirun \
    --verbose \
    -np 32 \
    pycbc_live \
    --size-override 1000 \
    --analysis-chunk 8 \
    --fftw-planning-limit 60 \
    --bank-file ${BANK_FILE} \
    --enable-bank-start-frequency \
    --sample-rate ${SR} \
    --max-length ${BS} \
    --pvalue-combination-livetime ${PVAL_LIVETIME} \
    --start-time ${RUN_START_TIME} \
    --end-time ${RUN_END_TIME} \
    --sync \
    --low-frequency-cutoff ${INJ_FLOW} \
    --approximant "SPAtmplt:mtotal<4" "SEOBNRv4_ROM:else" \
    --chisq-bins "0.9 * get_freq('fSEOBNRv4Peak',params.mass1,params.mass2,params.spin1z,params.spin2z) ** (2.0 / 3.0)" \
    --snr-abort-threshold 5000 \
    --snr-threshold 5.0 \
    --newsnr-threshold 5.0 \
    --max-triggers-in-batch 50 \
    --store-loudest-index 50 \
    --autogating-threshold 50 \
    --autogating-pad 0.5 \
    --autogating-cluster 0.5 \
    --autogating-width 0.25 \
    --autogating-taper 0.25 \
    --highpass-frequency 13 \
    --highpass-bandwidth 5 \
    --highpass-reduction 200 \
    --psd-samples 30 \
    --max-psd-abort-distance 6000 \
    --min-psd-abort-distance 20 \
    --psd-abort-difference .15 \
    --psd-recalculate-difference .01 \
    --psd-inverse-length 3.5 \
    --psd-segment-length 4 \
    --trim-padding .5 \
    --channel-name H1:DCS-CALIB_STRAIN_C02 L1:DCS-CALIB_STRAIN_C02 V1:Hrec_hoft_16384Hz \
    --state-channel H1:DSC-CALIB_STATE_VECTOR_C02 L1:DSC-CALIB_STATE_VECTOR_C02 V1:DQ_ANALYSIS_STATE_VECTOR \
    --increment-update-cache H1:./H1/ L1:./L1/  V1:./V1/ \
    --frame-src H1:./H1/* L1:./L1/* V1:./V1/* \
    --processing-scheme cpu:2 \
    --fftw-measure-level 0 \
    --fftw-threads-backend openmp \
    --increment 8 \
    --single-newsnr-threshold 8 \
    --single-duration-threshold 1.0 \
    --single-fixed-ifar .06 \
    --single-reduced-chisq-threshold 4.0 \
    --background-statistic phasetd_newsnr_sgveto \
    --background-statistic-files \
        /work/bhooshan.gadre/share/H1L1-PHASE_TIME_AMP_v1.hdf \
        /work/bhooshan.gadre/share/H1V1-PHASE_TIME_AMP_TITO.hdf \
        /work/bhooshan.gadre/share/L1V1-PHASE_TIME_AMP_TITO.hdf \
    --analyze-flags SCIENCE_INTENT \
    --verbose \
    --max-batch-size 16777216 \
    --frame-read-timeout 10 \
    --output-path ./data \
    --day-hour-output-prefix \
    --enable-background-estimation \
    --background-ifar-limit 100 \
    --timeslide-interval 0.1 \
    --ifar-upload-threshold .0003 \
    --ifar-double-followup-threshold 0.0003 \
    --round-start-time 4 \
    --enable-single-detector-background --verbose \
    --sgchisq-snr-threshold 4 \
    --sgchisq-locations "mtotal>40:20-30,20-45,20-60,20-75,20-90,20-105,20-120"

'''
import os
import logging
import argparse
import tqdm
from pycbc import frame

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--channel-names",
                        type=str, nargs="+",
                        help="list of frame files")
parser.add_argument("--gps-start-time",
                        help="The gps start time of the data "
                                "(integer seconds)", type=int)
parser.add_argument("--gps-end-time",
                        help="The gps end time of the data "
                                " (integer seconds)", type=int)
parser.add_argument("--frame-type",
                        type=str,
                        help="(optional), replaces frame-files. Use datafind "
                                "to get the needed frame file(s) of this type.")
parser.add_argument('--frame-duration', type=int, default=4,
                    help='Split all data into smaller frame files of the given duration if specified.')
parser.add_argument('--outdir', type=str, help='prefix to save strain files', default='.')
args = parser.parse_args()

logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)
start = args.gps_start_time
stop = args.gps_end_time
step = args.frame_duration

filename = '{0}-{1}-{2}-{3}.gwf'

if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

data = {}
det_channels = {}
for channel in args.channel_names:
    det_channels[channel.split(':')[0]] = []
    data[channel.split(':')[0]] = []

for channel in args.channel_names:
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
