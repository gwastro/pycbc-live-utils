import argparse
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
from ligo.gracedb.rest import GraceDb as gdb
import os
import numpy as np

# This is the only pycbc import at the moment - it might be nice to get this into a
# LAL (or similar) function, so the folks at gwcelery are more likely to adopt this
from pycbc.waveform.spa_tmplt import spa_length_in_time

parser = argparse.ArgumentParser()
playground_server = 'https://gracedb-playground.ligo.org/api/'
parser.add_argument('--gracedb-server', default=playground_server,
                    help="Server of the gracedb instance to use for getting "
                         "superevent details. Default = playground server")
parser.add_argument('--superevent-id', required=True,
                    help="ID of the superevent to be plotted")
parser.add_argument('--output-dir', default=os.getcwd(),
                    help="Directory to output the plots into, defult = current directory")
parser.add_argument('--ifar-thresh', type=float, default=0,
                    help="Threshold to not plot events with FAR less than this IFAR")
parser.add_argument('--latency-limit', type=float, default=100,
                    help="Do not plot anything after this point in latency "
                          "if given (seconds). Default=100")
args = parser.parse_args()

client = gdb(service_url=args.gracedb_server)
print("Pinging GDB server")
client.ping()
response = client.superevent(args.superevent_id).json()

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

print("Getting preferred event info")
# Get the information on all g_events and the preferred event
g_preferred = response['preferred_event']
gevent_list = response['gw_events']

# Get preferred event information
preferred_e = client.event(g_preferred).json()
pref_coinc_insp = preferred_e['extra_attributes']['CoincInspiral']
preferred_time = pref_coinc_insp['end_time'] \
                   + pref_coinc_insp['end_time_ns'] * 1e-9
preferred_snr = pref_coinc_insp['snr']
preferred_latency = preferred_e['reporting_latency']
preferred_pipeline = preferred_e['pipeline'].lower()
preferred_search = preferred_e['search']

preferred_ifar = (1. / preferred_e['far']) / (86400 * 365.25)
if preferred_ifar < 1e-3:
    pifar_str = f"{preferred_ifar:.3e}"
else:
    pifar_str = f"{preferred_ifar:.3f}"

if preferred_ifar < args.ifar_thresh:
    print(f"IFAR of {pifar_str} is less than " +
          f"threshold {args.ifar_thresh}. Exiting")
    exit()

print("Getting template end time information")
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

preferred_pmt = -premerger_time(preferred_e)

print("Getting all other events")
# Save all (non-preferred) events into a dictionary, keyed on pipeline and search
all_events = {}
for g in gevent_list:
    # Do not deal with preferred event in this loop
    if g == g_preferred: continue
    e = client.event(g).json()
    pipeline, search = e['pipeline'].lower(), e['search']
    latency = e['reporting_latency']
    snr = e['extra_attributes']['CoincInspiral']['snr']
    prem_time = -premerger_time(e) 
    if pipeline not in all_events:
        all_events[pipeline] = {}
    if search not in all_events[pipeline]:
        all_events[pipeline][search] = [(g, latency, snr, prem_time)]
    else:
        all_events[pipeline][search].append((g, latency, snr, prem_time))

pmt_color = '#707070'
pmt_linestyle = '--'

print("Plotting events")
# plot all events, with colours and edges to define the source
fig, ax = plt.subplots(1, figsize=(9.5,4.5))
for pipeline in all_events.keys():
    for search in all_events[pipeline]:
        g, latency, snr, prem_time = zip(*all_events[pipeline][search])
        g = np.array(g)
        latency = np.array(latency)
        snr = np.array(snr)
        prem_time = np.array(prem_time)

        latency_above = latency > args.latency_limit
        latency_to_plot = latency[:]
        latency_to_plot[latency_above] = args.latency_limit
        ax.plot([latency_to_plot, prem_time], [snr, snr],
                linestyle=pmt_linestyle, c=pmt_color)
        ax.scatter(latency_to_plot, snr,
                   color=pipelinecolours[pipeline],
                   edgecolors=searchedges[search],
                   marker='o', s=50, zorder=30)
        if any(latency_above):
            arr_start = latency_to_plot[latency_above]
            snr_above = snr[latency_above]
            arr_end = arr_start + 5
            ax.plot([arr_start, arr_end], [snr_above, snr_above],
                    lw=2, color=pipelinecolours[pipeline])
            ax.scatter([arr_end],[snr_above],
                       color=pipelinecolours[pipeline],
                       marker='>', zorder=100)

# Want to make the legend general - plot some default colours off the
# edge of the plot for making the legends
leg1_lines = []
leg1_labels = []
# Always plot all pipelines in the legend, even if not contributed to event
# Make these points bigger for clarity
for pipeline in pipelinenames.keys():
    leg1_lines.append(ax.scatter([0], [-20], s=80,
                                 color=pipelinecolours[pipeline]))
    leg1_labels.append(f"{pipelinenames[pipeline]}")

# Show how the markers are modified for early-warning and preferred events
default_color = '#D0D0D0'
leg2_lines = []
leg2_labels = []
for search in searchnames.keys():
    leg2_lines.append(ax.scatter([0],[-20], color=default_color, s=80,
                                 edgecolors=searchedges[search]))
    leg2_labels.append(f"{searchnames[search]}")

pmt_plot, = ax.plot([0,0], [-20,-20],
                    linestyle=pmt_linestyle, c=pmt_color)
leg2_lines.append(pmt_plot)
leg2_labels.append("Time from template end")

# Plot preferred scatter point
ax.plot([preferred_latency, preferred_pmt],
        [preferred_snr, preferred_snr],
        linestyle=pmt_linestyle, c=pmt_color, zorder=-30)
ax.scatter(preferred_latency, preferred_snr,
           color=pipelinecolours[preferred_pipeline],
           edgecolors=searchedges[preferred_search],
           marker='*', s=150)

# Add to the modify legend
leg2_lines.append(ax.scatter([0], [-20],
                             color=default_color,
                             marker='*', s=150))
leg2_labels.append("Preferred Event")

# Cut off the points used for adding to the legend
ax.set_ylim(bottom=4)

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
ax.set_title(f"{args.superevent_id} upload timeline")
ax.set_xlabel(f"Time from {preferred_time:.3f} [s]")
ax.set_ylabel('Network SNR')

# Add dotted grid
ax.grid(visible=True, which='major', linestyle='--', zorder=-50)

# Need to adjust the size of the plot on the figure to get the legends in
fig.subplots_adjust(left=0.1, right=0.75)

filetrunk = os.path.join(args.output_dir, f"{args.superevent_id}_timeline")
fig.savefig(filetrunk + '.pdf')
fig.savefig(filetrunk + '.png')
