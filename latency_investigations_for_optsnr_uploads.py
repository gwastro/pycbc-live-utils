#!/usr/bin/env python

"""
Code used for investigating usefulness of optimised SNR uploads
"""

import argparse
import matplotlib
import pickle
import json
import numpy as np
from matplotlib import pyplot as plt
from ligo.gracedb.rest import GraceDb as gdb
from astropy.time import Time

matplotlib.use('agg')

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--collect-data', action='store_true', default=False,
                    help='It is convenient to query and save all the trig data only once. \
Run using this flag to query and store the data first, and then without it to generate the plots from saved data.')
parser.add_argument('--latest-only', action='store_true', default=True,
                    help='Consider events after 10/12/22.')
parser.add_argument('--datafolder', type=str, default='.',
                    help='Path to Data folder.')
parser.add_argument('--plotsfolder', type=str, default='.',
                    help='Path to Plots folder.')
args = parser.parse_args()

collect_keys = ['graceid', 'pipeline', 'gpstime', 
                'reporting_latency', 'instruments', 'nevents', 
                'far', 'likelihood', 'superevent', 'search']

# Query and save the data
if collect_data:

    playground_server = 'https://gracedb-playground.ligo.org/api/'
    client = gdb(service_url=playground_server)
    client.ping()

    events = client.events("group: CBC pipeline: pycbc")

    collected_data = []
    i=0
    for e in events:

        if i%100 == 0: print(i)

        collected_data.append(e)

        i += 1
        if i == 100: break

    with open(f'{args.datafolder}/pycbclive_events_fulldata.pkl', 'wb') as f:
        pickle.dump(collected_data, f, protocol=5)

    exit()


# Assuming the data has already been saved, read the saved data
data = {}
for k in collect_keys:
    data[k] = []

data['coinc_min_duration'] = []
data['sngl_template_duration'] = []
data['chisq_dof'] = []
data['chisq'] = []
data['index'] = []

with open(f'{args.datafolder}/pycbclive_events_fulldata.pkl', 'rb') as f:
    tmp_data = pickle.load(f)

for i, d in enumerate(tmp_data):
    if i%10000 == 0: print(i)

    # skip the unwanted events
    if 'SingleInspiral' not in d['extra_attributes']:
        continue

    if (d['extra_attributes']['SingleInspiral'][0]['channel'] not in ['GDS-CALIB_STRAIN_INJ1_O3Replay', 'Hrec_hoft_16384Hz_INJ1_O3Replay']) \
        or (d['offline'] == True) \
        or (tmp_data[i]['submitter'] != 'pycbclive'):
        continue

    for j, singl_attr in enumerate(d['extra_attributes']['SingleInspiral']):
        if 'chisq' in singl_attr:

            data['index'].append(i)

            for k in collect_keys:
                data[k].append(d[k])
            data['coinc_min_duration'].append(d['extra_attributes']['CoincInspiral']['minimum_duration'])
            data['sngl_template_duration'].append(d['extra_attributes']['SingleInspiral'][0]['template_duration'])

            data['chisq_dof'].append(d['extra_attributes']['SingleInspiral'][j]['chisq_dof'])
            data['chisq'].append(d['extra_attributes']['SingleInspiral'][j]['chisq'])
            break

for k in data:
    data[k] = np.array(data[k])


# checking if we are looking at the correct events 
# first five events with max latency
isort = np.argsort(data['reporting_latency'])
i5max = data['index'][isort][-5:]

for i in i5max:
    print(json.dumps(tmp_data[i], sort_keys=True, indent=4))

for i in data['index']:
    # just to be sure if there are any non-pycbclive events
    if tmp_data[i]['submitter'] != 'pycbclive':
        print(json.dumps(tmp_data[i], sort_keys=True, indent=4))
        print(i)

# Plots
if latest_only:
    # consider events after 10 Dec '22
    suff = '_after10dec22'
    s = Time("2022-12-10T00:00:00", format='isot', scale='utc')
else:
    suff = ''
    s = Time("2000-01-01T00:00:00", format='isot', scale='utc')

# sanity check that there are no ridiculously late events
tmax = 24*3600 # 1 day
assert sum(data['reporting_latency'] > tmax) == 0

# catch optimised SNR events
i_opt_SNR = (data['chisq']==1) & (data['chisq_dof']==1) & (data['gpstime'] > s.gps)
# early warning ones
i_earlywarn = (data['search'] == 'EarlyWarning') & (data['gpstime'] > s.gps)
# Allsky ones
i_allsky = (data['search'] == 'AllSky') & (data['gpstime'] > s.gps)

# i_opt_SNR is a whole subset of i_allsky
i_allsky_excl = np.logical_and(i_allsky, np.logical_not(i_opt_SNR))

assert sum(i_opt_SNR) == sum(np.logical_and(i_allsky, i_opt_SNR))

Nabove = sum(data['reporting_latency'][i_opt_SNR]>270)
Nbelow = sum(data['reporting_latency'][i_opt_SNR]<=270)

# Cumulative histograms
plt.figure(figsize=(10,8))
# plt.loglog(np.sort(data['reporting_latency'][i_opt_SNR]), y, 'x', markersize=4)
plt.step(np.sort(data['reporting_latency'][i_allsky_excl]),
         (np.arange(sum(i_allsky_excl))[::-1]+1)/sum(i_allsky_excl), 
         label='AllSky_excl')
plt.step(np.sort(data['reporting_latency'][i_allsky]),
         (np.arange(sum(i_allsky))[::-1]+1)/sum(i_allsky), 
         label='AllSky')
plt.step(np.sort(data['reporting_latency'][i_opt_SNR]),
         (np.arange(sum(i_opt_SNR))[::-1]+1)/sum(i_opt_SNR), 
         label='Optimised SNR events')
plt.step(np.sort(data['reporting_latency'][i_earlywarn]),
         (np.arange(sum(i_earlywarn))[::-1]+1)/sum(i_earlywarn), 
         label='EarlyWarning')
plt.axvline(x=270, ymin=0, ymax=1, color='r')
plt.xscale('log')
plt.yscale('log')
plt.xlim([1,3000])
plt.xlabel('Reporting Latency (sec)')
plt.ylabel('Cumulative fraction of events')
plt.legend()
plt.title(f'Optimised events: {Nabove} with latency > 270 sec and {Nbelow} with latency <= 270 sec')
plt.savefig(f'{args.plotsfolder}/pycbc_live_event_latency_fullrange{suff}.png')

# Scatter: reporting latency vs template length
plt.figure(figsize=(10,10))
# plt.loglog(np.sort(data['reporting_latency'][i_opt_SNR]), y, 'x', markersize=4)
plt.step(data['sngl_template_duration'][i_allsky_excl],
         data['reporting_latency'][i_allsky_excl], 'co', alpha=0.3, markersize=2, 
         label='AllSky_excl')
plt.step(data['sngl_template_duration'][i_earlywarn],
         data['reporting_latency'][i_earlywarn], 'bo', alpha=0.3, markersize=2, 
         label='EarlyWarning')
plt.step(data['sngl_template_duration'][i_opt_SNR],
         data['reporting_latency'][i_opt_SNR], 'rx', markersize=2, 
         label='Optimised SNR events')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Template length (sec)')
plt.ylabel('Reporting Latency (sec)')
plt.legend()
# plt.title(f'Opt events: {Nabove} events with latency > 270 sec and {Nbelow} events with latency <= 270 sec')
plt.savefig(f'{args.plotsfolder}/pycbc_live_event_templen_latency{suff}.png')


# only early warning ones
i_earlywarn_zm = (data['search'] == 'EarlyWarning') & (data['gpstime'] > s.gps) #& (data['reporting_latency'] <= 0)

# Cumulative hist
plt.figure(figsize=(10,8))
plt.step(np.sort(data['reporting_latency'][i_earlywarn_zm]),
         (np.arange(sum(i_earlywarn_zm))+1)/sum(i_earlywarn_zm), 
         label='EarlyWarning')
# plt.axvline(x=270, ymin=0, ymax=1, color='r')
# plt.xscale('log')
plt.yscale('log')
plt.xlim([-50,30])
plt.xlabel('Reporting Latency (sec)')
plt.ylabel('Cumulative fraction of events')
plt.legend(loc='upper left')
plt.grid(which='both')
plt.savefig(f'{args.plotsfolder}/pycbc_live_event_latency_fullrange_earlywarning_zoom{suff}.png')

# Scatter: Reporting Latency vs Template length
plt.figure(figsize=(8,8))
plt.step(data['sngl_template_duration'][i_earlywarn_zm],
         data['reporting_latency'][i_earlywarn_zm], 'bo', alpha=0.3, markersize=2, 
         label='EarlyWarning')
# plt.xscale('log')
# plt.yscale('log')
plt.xlabel('Template length (sec)')
plt.ylabel('Reporting Latency (sec)')
plt.legend()
# plt.title(f'Opt events: {Nabove} events with latency > 270 sec and {Nbelow} events with latency <= 270 sec')
plt.savefig(f'{args.plotsfolder}/pycbc_live_event_templen_latency_earlywarn{suff}.png')

