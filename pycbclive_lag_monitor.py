#!/usr/bin/env python

import argparse
import os
import glob
import time
import datetime as dt
import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as pl
import matplotlib.dates as md


def get_local_tz_name():
    """Return the system's timezone name."""
    utc_now = dt.datetime.now(dt.timezone.utc)
    local_now = utc_now.astimezone()
    return local_now.tzname()


parser = argparse.ArgumentParser()
parser.add_argument("--log-glob", required=True)
parser.add_argument("--output-path", required=True)
parser.add_argument("--process-rank", default="0")
args = parser.parse_args()

pl.rcParams["figure.figsize"] = 15, 7
pl.rcParams["font.size"] = 12

today = time.strftime("%Y-%m-%d")

# the time zone is not given in PyCBC Live's log,
# so assume it matches the current one
tz = get_local_tz_name()

data = []

for input_file in glob.glob(args.log_glob):
    for line in open(input_file, "r"):
        if "duty" not in line:
            # only interested in lines reporting the stats
            continue
        fields = line.split()
        if len(fields) != 16 or len(fields[0]) != 10:
            # a few safety checks
            continue
        if fields[3] != args.process_rank:
            # only look at lines for the given process
            continue
        if fields[0] < today:
            # only consider today's log
            continue

        timestamp = "{} {}".format(
            fields[0], fields[1].split(",")[0]
        )  ## timestamp -- Date and Time

        duty = float(fields[9].replace(",", ""))
        lag = float(fields[11])
        ndet = int(fields[13])

        data.append([timestamp, duty, lag, ndet])  ## Array for plotting

if not data:
    exit()

data = sorted(data)  ##sorting time-stamp wise

lags = np.array([lag for _, _, lag, _ in data])
ndets = np.array([ndet for _, _, _, ndet in data])

tstamp = np.array([timestamp for timestamp, _, _, _ in data])

day_time_array = np.array(
    [
        dt.datetime.strptime(date_time_str, "%Y-%m-%d %H:%M:%S")
        for date_time_str in tstamp
    ]
)
daysarray = np.array([i.date() for i in day_time_array])

no_days = max(daysarray) - min(daysarray)

date_list = [min(daysarray) + dt.timedelta(days=x) for x in range(0, no_days.days + 2)]
datetime_list_for_ticks = [dt.datetime.combine(i, dt.time(00)) for i in date_list]

color = ["red", "orange", "blue", "green"]

for i in range(0, len(datetime_list_for_ticks) - 1):
    index_list_datewise = np.argwhere(
        (day_time_array >= datetime_list_for_ticks[i])
        & (day_time_array <= datetime_list_for_ticks[i + 1])
    ).flatten()

    datetime_array_daywise = day_time_array[index_list_datewise]
    lags_daywise = lags[index_list_datewise]
    ndets_daywise = ndets[index_list_datewise]

    fig, ax = pl.subplots(2, 1, sharex=True)

    ax[0].plot(datetime_array_daywise, lags_daywise, zorder=1, lw=0.5, color="#808080")

    for ndet in np.unique(ndets):
        l = lags_daywise[ndets_daywise == ndet]
        day_time_array_det = datetime_array_daywise[ndets_daywise == ndet]
        ax[0].scatter(
            day_time_array_det,
            l,
            label="{} live detector(s)".format(ndet),
            marker=".",
            s=3,
            color=color[ndet],
            zorder=2,
        )

    xfmt = md.DateFormatter("%Y-%m-%d\n%H:%M:%S " + tz)
    ax[0].xaxis.set_major_formatter(xfmt)
    ax[0].tick_params(axis="x")

    title = date_list[i].isoformat()
    ax[0].set_xlim(datetime_array_daywise[0], datetime_list_for_ticks[i + 1])
    ax[0].set_ylim(1, 400)
    ax[0].set_yscale("log")
    ax[0].grid(True)
    ax[0].set_ylabel("PyCBC Live triggers-to-disk lag [s]")
    lgnd = ax[0].legend(ncol=2, loc="best")
    ax[0].set_title(title)

    ax[1].plot(datetime_array_daywise, ndets_daywise)
    xfmt = md.DateFormatter("%Y-%m-%d\n%H:%M:%S " + tz)
    ax[1].xaxis.set_major_formatter(xfmt)
    ax[1].tick_params(axis="x")
    ax[1].set_yticks([0, 1, 2, 3])
    ax[1].set_ylabel("Live detectors' number")
    ax[1].grid(True)

    pl.tight_layout()

    # figure out and prepare the output directory
    out_base = os.path.join(args.output_path, date_list[i].strftime("%Y/%m/%d"))
    if not os.path.exists(out_base):
        os.makedirs(out_base)

    # save the plot
    out_file_name = "{}_lag_over_time.png".format(date_list[i])
    out_path = os.path.join(out_base, out_file_name)
    fig.savefig(out_path, dpi=200)

    # save the data to an ascii file in case we need it later
    out_file_name = "{}_lag_over_time.txt".format(date_list[i])
    data_txt_path = os.path.join(out_base, out_file_name)
    with open(data_txt_path, "w") as df:
        df.write("# date time_{} duty_factor lag num_active_detectors\n".format(tz.lower()))
        for entry in data:
            df.write("{} {:.2f} {:.2f} {}\n".format(*entry))
