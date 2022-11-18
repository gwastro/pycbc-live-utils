#!/usr/bin/env python

"""
A plotting script for superevents in GraceDB, showing early warning and allsky events from all pipelines and their latency compared to the event time
"""

import argparse
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
from lal.gpstime import gps_to_utc, utc_to_gps
from ligo.gracedb.rest import GraceDb as gdb
import os
import copy
import numpy as np
import datetime as dt
from datetime import datetime as dtdt
import logging

# This is the only pycbc import at the moment - it might be nice to get this into a
# LAL (or similar) function, so the folks at gwcelery are more likely to adopt this
from pycbc.waveform.spa_tmplt import spa_length_in_time
from pycbc import init_logging
# Set up dictionaries to use for the different pipelines and searches

# Colours chosen to be colourblind friendly.
# Some editing to do to make it grayscale-friendly
pipelinecolours = {
    'cwb':     '#AA4499',
    'gstlal':  '#F0E442',
    'mbta':    '#D55E00',
    'pycbc':   '#56B4E9',
    'spiir':   '#009E73',
}


_marker = {
    'standard': 'o',
    'highlight': '*',
    'comment': '^',
    'file': 's',
}

default_color = '#D0D0D0'
pmt_color = '#707070'
pmt_linestyle = '--'


pipelinenames = {
    'cwb':     'cWB',
    'gstlal':  'GstLAL',
    'mbta':    'MBTA',
    'pycbc':   'PyCBC',
    'spiir':   'SPIIR',
}

searchnames = {
    'AllSky': 'Full Bandwidth',
    'EarlyWarning': 'Early Warning',
}

searchedges = {
    'AllSky': 'face',
    'EarlyWarning': 'k',
}

parser = argparse.ArgumentParser()
playground_server = 'https://gracedb-playground.ligo.org/api/'
parser.add_argument('--gracedb-server', default=playground_server,
                    help="Server of the gracedb instance to use for getting "
                         "superevent details. Default = playground server")
parser.add_argument('--superevent-id',
                    help="ID of the superevent to be plotted")
parser.add_argument('--event-id',
                    help="ID of the event to plot nearby events")
parser.add_argument('--event-search-window', type=float, default=0.2,
                    help="If the --event-id option is used, how big "
                         "is the window (seconds) to search for nearby "
                         "events? Default 0.2")
parser.add_argument('--output-dir', default=os.getcwd(),
                    help="Directory to output the plots into. "
                         "Default = current directory")
parser.add_argument('--latency-limit', type=float, default=100,
                    help="Do not plot anything after this point in latency "
                         "(seconds). Default=100")
parser.add_argument("--include-test", action='store_true',
                    help="Include events marked as test events in the query")
parser.add_argument("--include-mdc", action='store_true',
                    help="Include events marked as MDC events in the query")
parser.add_argument("--pipeline-only", choices=list(pipelinenames.keys()),
                    help="Only plot events form this pipeline")
parser.add_argument("--search-only", choices=list(searchnames.keys()),
                    help="Only plot events form this search")
parser.add_argument("--verbose", action='store_true',
                    help="Print logging statements")
args = parser.parse_args()

init_logging(args.verbose)

if (not args.superevent_id) == (not args.event_id):
    # The ID arguments are either both given or both not given
    parser.error("One of --superevent-id or --event-id must be given, "
                 "but not both.")

client = gdb(service_url=args.gracedb_server)
logging.info("Pinging GDB server")
client.ping()

# Get the time difference between the final merger frequency and when this
# template ended
def premerger_time(e_event, f_final=None):
    snglinsps = e_event['extra_attributes']['SingleInspiral']
    for sngl in snglinsps:
        if 'f_final' in sngl:
            snglinsp = sngl
            break
    else:
        # not an early warning candidate, return zero:
        return 0
    f_final = f_final if f_final else snglinsp['f_final']
    pm_t = spa_length_in_time(mass1=snglinsp['mass1'],
                              mass2=snglinsp['mass2'],
                              f_lower=f_final,
                              phase_order=-1)
    return pm_t

def get_event_info(event, central_time):
    g = event['graceid']
    log_times = {k: [] for k in ['file', 'comment']}
    for log in client.logs(event['graceid']).json()['log']:
        dt_log = dtdt.strptime(log['created'], "%Y-%m-%d %H:%M:%S %Z")
        tlog = float(utc_to_gps(dt_log) - central_time)
        # Original upload / creation do not get plotted
        if any([substr in log['comment']
                for substr in ["Original", "Created"]]):
            continue
        if log['filename'] == '':
            log_times['comment'].append(tlog)
        else:
            log_times['file'].append(tlog)
    latency = event['reporting_latency']
    snr = event['extra_attributes']['CoincInspiral']['snr']
    prem_time = -premerger_time(event)
    return (g, latency, snr, prem_time, log_times)



if args.superevent_id:
    response = client.superevent(args.superevent_id).json()
    logging.info("Getting highlight event info")
    g_highlight = response['preferred_event']
    highlight_e = client.event(g_highlight).json()
    central_time = pref_coinc_insp['end_time'] \
                       + pref_coinc_insp['end_time_ns'] * 1e-9
    pref_event_info = get_event_info(highlight_e, central_time)
    highlight_pipeline = highlight_e['pipeline'].lower()
    highlight_search = highlight_e['search']

    # Get the list of g-events associated with the superevent
    gevent_list = response['gw_events']
    gevent_list.delete(g_highlight)

else:
    # Get the event time of the given event
    g_highlight = args.event_id
    highlight_e = client.event(args.event_id).json()
    # Use the original event time as central
    central_time = highlight_e['gpstime']
    pref_event_info = get_event_info(highlight_e, central_time)
    highlight_pipeline = highlight_e['pipeline'].lower()
    highlight_search = highlight_e['search']
    query = f"{central_time - args.event_search_window} .. " + \
                f"{central_time + args.event_search_window}"
    gevent_list = [r['graceid'] for r in client.events(query=query)]
    if args.include_test:
        test_query = "group: test " + query
        gevent_list += [r['graceid'] for r in client.events(query=test_query)]
    if args.include_mdc:
        mdc_query = "group: mdc " + query
        gevent_list += [r['graceid'] for r in client.events(query=mdc_query)]

    # filter out the event if interest
    gevent_list.remove(g_highlight)

logging.info(f"Found {len(gevent_list)} events")

logging.info("Getting nearby event info")

# Save all events into a dictionary, keyed on pipeline and search
all_events = {pip: {s: [] for s in searchnames.keys()}
              for pip in pipelinenames.keys()}

counter = 0
for g in gevent_list:
    e = client.event(g).json()
    pipeline, search = e['pipeline'].lower(), e['search']
    if args.pipeline_only and not pipeline == args.pipeline_only: continue
    if args.search_only and not search == args.search_only: continue
    counter += 1
    gevent_info = get_event_info(e, central_time)
    all_events[pipeline][search].append(gevent_info)

logging.info(f"{counter} events to plot")

logging.info("Plotting events")
# plot all events, with colours and edges to define the source
fig, ax = plt.subplots(1, figsize=(9.5,4.5))

def add_to_plot(ax, event_info, pipeline, search, xlim, highlight=False):
    marker = _marker['highlight'] if highlight else _marker['standard']
    s = 100 if highlight else 50

    g, latency, snr, prem_time, log_times = event_info
    latency_above = latency > args.latency_limit

    latency_to_plot = min(latency, args.latency_limit)
    ax.plot([latency_to_plot, prem_time], [snr, snr],
            linestyle=pmt_linestyle, c=pmt_color)
    ax.scatter(latency_to_plot, snr,
               color=pipelinecolours[pipeline],
               edgecolors=searchedges[search],
               marker=marker, s=s, zorder=30)
    xlim_new = [min(xlim[0], prem_time),
                max(xlim[1], latency_to_plot)]
    if latency_above:
        arr_start = latency_to_plot
        snr_above = snr
        arr_end = arr_start + 5
        ax.plot([arr_start, arr_end], [snr_above, snr_above],
                lw=2, color=pipelinecolours[pipeline])
        ax.scatter([arr_end],[snr_above],
                   color=pipelinecolours[pipeline],
                   marker='>', zorder=100)
    for k in ['comment', 'file']:
        xlim_new = [xlim_new[0],
                    max(xlim_new[1], max(log_times[k]))]
        ax.scatter(log_times[k], snr * np.ones_like(log_times[k]),
                   color=pipelinecolours[pipeline],
                   marker=_marker[k], s=10,
                   edgecolors=searchedges[search],)
    return xlim_new


xlim = [0, 0]

for pipeline in all_events.keys():
    for search in all_events[pipeline]:
        if all_events[pipeline][search] == []: continue
        events_info = tuple(zip(*all_events[pipeline][search]))
        for event_info in zip(*events_info):
            xlim = add_to_plot(ax, event_info, pipeline, search, xlim)

# Plot highlight scatter point
xlim = add_to_plot(ax, pref_event_info, highlight_pipeline,
                   highlight_search, xlim, highlight=True)

xlim = [xlim[0] - 6, min(xlim[1], args.latency_limit) + 6]
ylim_orig = ax.get_ylim()

# Want to make the legend general - plot some default colours off the
# edge of the plot for making the legends
leg1_lines = []
leg1_labels = []
leg2_lines = []
leg2_labels = []

# Always plot all pipelines in the legend, even if not contributed to event
# Make these points bigger for clarity
for pipeline in pipelinenames.keys():
    if args.pipeline_only and not pipeline == args.pipeline_only: continue
    leg1_lines.append(ax.scatter([0], [-20], s=80,
                                 color=pipelinecolours[pipeline]))
    leg1_labels.append(f"{pipelinenames[pipeline]}")


leg2_lines.append(ax.scatter([0], [-20],
                             color=default_color,
                             marker=_marker['highlight'], s=150))
# Add to the modify legend
if args.superevent_id:
    leg2_labels.append("Preferred Event")
else:
    leg2_labels.append(args.event_id)

# Show how the markers are modified for early-warning (and highlight) events
for search in searchnames.keys():
    if args.search_only and not search == args.search_only: continue

    leg2_lines.append(ax.scatter([0],[-20], color=default_color, s=80,
                                 edgecolors=searchedges[search]))
    leg2_labels.append(f"{searchnames[search]}")

for k in ['comment', 'file']:
    leg2_lines.append(ax.scatter([0],[-20],
                      color=default_color, marker=_marker[k]))
    leg2_labels.append(f"{k} added")

pmt_plot, = ax.plot([0,0], [-20,-20],
                    linestyle=pmt_linestyle, c=pmt_color)
leg2_lines.append(pmt_plot)
leg2_labels.append("Time from template end")

# Cut off the points used for adding to the legend
ax.set_ylim(bottom=max(4, ylim_orig[0]), top=ylim_orig[1])
ax.set_xlim(xlim)

# Set up the legends - want to be outside the plot so it never overlaps
# with the data
leg2 = plt.legend(leg2_lines, leg2_labels,
                  loc='lower left',
                  bbox_to_anchor=(1.01, 0))
leg1 = plt.legend(leg1_lines, leg1_labels,
                  loc='upper left',
                  bbox_to_anchor=(1.01, 1))
ax.add_artist(leg2)

# Indicate the time of the event (always zero on this scale)
ax.axvline(0, color='r',linestyle='--')

# Informative title and axes
if args.superevent_id:
    ax.set_title(f"{args.superevent_id} upload timeline")
else:
    ax.set_title(f"Upload timeline for events around {args.event_id}")
central_utc = gps_to_utc(central_time)
central_time_str = central_utc.strftime("%Y-%m-%d %H:%M:%S.") + \
    f"{central_utc.microsecond // 1000:03d}"
ax.set_xlabel(f"Time from {central_time_str}")
ax.set_ylabel('Network SNR')

# Add dotted grid
ax.grid(visible=True, which='major', linestyle='--', zorder=-50)

# Need to adjust the size of the plot on the figure to get the legends in
fig.subplots_adjust(left=0.1, right=0.75)


id_str = args.superevent_id if args.superevent_id else args.event_id
filetrunk = os.path.join(args.output_dir, f"{id_str}_timeline")
fig.savefig(filetrunk + '.pdf')
fig.savefig(filetrunk + '.png')
