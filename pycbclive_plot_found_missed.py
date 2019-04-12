import sys
import glob
import tqdm
import numpy as np
import pylab as pl
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import ligolw, table, lsctables


class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
    pass

lsctables.use_in(LIGOLWContentHandler)

# read injections

sim_path = '/home/tito/moving_from_atlas/projects/pycbc/bayestar_snr_timeseries/fake_strain_injections_live_paper/H1L1V1-INJECTIONS_ALL.xml.gz'
doc = ligolw_utils.load_filename(sim_path, False, contenthandler=LIGOLWContentHandler)
sim_table = table.get_table(doc, lsctables.SimInspiralTable.tableName)
injections = []
for sim in sim_table:
    injections.append((float(sim.get_time_geocent()),
                       sim.mchirp,
                       sim.alpha1,
                       sim.alpha2,
                       sim.alpha3))
del doc
print('{} injections'.format(len(injections)))

# read triggers

triggers = []
for file_path in tqdm.tqdm(glob.glob(sys.argv[1])):
    doc = ligolw_utils.load_filename(file_path, False, contenthandler=LIGOLWContentHandler)
    coinc_table = table.get_table(doc, lsctables.CoincInspiralTable.tableName)
    end_time, far, mchirp = coinc_table[0].end_time + coinc_table[0].end_time_ns * 1e-9, coinc_table[0].combined_far, coinc_table[0].mchirp
    triggers.append([end_time, far, mchirp])
    del doc
triggers = np.array(triggers)
print('{} triggers'.format(triggers.shape[0]))

# associate

found = []
for injt, mchirp, osnr_h, osnr_l, osnr_v in injections:
    delta_t = abs(triggers[:,0] - injt)
    delta_mchirp = abs(triggers[:,2] - mchirp) / mchirp
    i = np.argmin(delta_t)
    far = triggers[i,1] if (delta_t[i] < 0.5 and delta_mchirp[i] < 0.5) else np.nan
    found.append([injt, osnr_h, osnr_l, osnr_v, far])
found = np.array(found)

network_snr = np.sum(found[:,1:4] ** 2, axis=1) ** 0.5
decisive_snr = np.array([sorted(found[i,1:4])[1] for i in range(found.shape[0])])
found_mask = np.isfinite(found[:,4])
missed_mask = np.logical_not(found_mask)

title = '{}/{} found'.format(sum(found_mask), len(found_mask))

pl.figure(figsize=(14,7))

pl.plot(found[missed_mask,0] - found[0,0],
        decisive_snr[missed_mask], 'xr')
pl.scatter(found[found_mask,0] - found[0,0],
           decisive_snr[found_mask],
           c=np.log10(found[found_mask,4]), vmin=-9, vmax=-4)
pl.axhline(4.5, ls='--', color='k')
pl.yscale('log')
pl.xlabel('Time')
pl.ylabel('Decisive optimal SNR')
cb = pl.colorbar(extend='both')
cb.set_label('log10(FAR [Hz])')
pl.title(title)

pl.tight_layout()
pl.savefig('/home/tito/public_html/pycbc/live/o3_initial_review/test.png')
