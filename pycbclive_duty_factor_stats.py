#!/usr/bin/env python

"""Read the duty factors from a PyCBC Live log and show some stats about them."""

import argparse
import numpy as np
import matplotlib.pyplot as pp


def plot_cdf_nicely(x, **kwa):
    """A nicer version of pp.hist([...] cumulative=-1): it does not bin and
    does not plot a spurious vertical line at the lowest sample.
    """
    sorter = np.argsort(x)
    fraction = np.arange(len(x))[::-1] / len(x)
    pp.step(np.array(x)[sorter], fraction, **kwa)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--input-log', required=True)
parser.add_argument('--out-plot')
args = parser.parse_args()

duty_factors = {}

with open(args.input_log, 'r') as lf:
    for line in lf:
        if 'Took' not in line:
            continue
        pieces = line.split()
        rank = int(pieces[3])
        duty_factor = float(pieces[9].replace(',', ''))
        if rank not in duty_factors:
            duty_factors[rank] = []
        duty_factors[rank].append(duty_factor)

print('Duty factors:')
for rank in sorted(duty_factors):
    print('{:04d}: {:d} samples, median {:.2f}, 90th percentile {:.2f}'.format(
        rank,
        len(duty_factors[rank]),
        np.median(duty_factors[rank]),
        np.percentile(duty_factors[rank], 90)
    ))

if args.out_plot is not None:
    df_max = 0
    df_min = np.inf
    for rank in sorted(duty_factors):
        color = pp.cm.viridis(rank / max(duty_factors))
        lw = 2 if rank == 0 else 0.5
        plot_cdf_nicely(duty_factors[rank], color=color, lw=lw)
        df_max = max(df_max, max(duty_factors[rank]))
        df_min = min(df_min, min(duty_factors[rank]))

    pp.axvspan(
        1, df_max, hatch='/', facecolor='none', edgecolor='gray'
    )
    pp.xscale('log')
    pp.xlim(max(df_min, 0.1), df_max)
    pp.ylim(0, 1)
    pp.xlabel('Duty factor')
    pp.ylabel('Cumulative fraction')
    pp.title('Median at rank 0: {:.2f} {}'.format(
        np.median(duty_factors[0]),
        '\N{CHECK MARK}' if np.median(duty_factors[0]) < 1 else '\N{WARNING SIGN}'
    ))

    pp.tight_layout()
    pp.savefig(args.out_plot, dpi=200)
