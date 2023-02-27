import numpy as np
from math import isnan

# id commonName version catalog.shortName GPS reference jsonurl mass_1_source mass_1_source_lower mass_1_source_upper mass_2_source mass_2_source_lower mass_2_source_upper network_matched_filter_snr network_matched_filter_snr_lower network_matched_filter_snr_upper luminosity_distance luminosity_distance_lower luminosity_distance_upper chi_eff chi_eff_lower chi_eff_upper total_mass_source total_mass_source_lower total_mass_source_upper chirp_mass_source chirp_mass_source_lower chirp_mass_source_upper chirp_mass chirp_mass_lower chirp_mass_upper redshift redshift_lower redshift_upper far far_lower far_upper p_astro p_astro_lower p_astro_upper final_mass_source final_mass_source_lower final_mass_source_upper

# .txt file was obtained directly from GWOSC website
data = np.genfromtxt('gwosc_allevents_150722.txt', usecols=[0, 1, 2, 3, 13, 25, 28, 31, 34],
    dtype="U15,U15,f8,U20,f8,f8,f8,f8,f8", names=True, deletechars="'")

print(data.dtype)

evdata = {}
for row in data:
    gwname = row['commonName']
    if not 'GW' in gwname:
        print('Excluding %s' % gwname)
        continue
    # exclude if needed values do not exist
    if isnan(row['network_matched_filter_snr']) or \
        (isnan(row['chirp_mass']) and (isnan(row['chirp_mass_source']) or isnan(row['redshift']))):
        print('Not enough values for %s %s' % (gwname, row['catalog.shortName']))
        continue
    if gwname not in evdata:
        evdata[gwname] = [
            row['network_matched_filter_snr'],
            row['chirp_mass'],
            row['chirp_mass_source'],
            row['redshift'],
            row['catalog.shortName']
        ]
        print(gwname, evdata[gwname])
    # if there was a NaN in the SNR column, or
    # if the new event is a duplicate with higher SNR, or updating GWTC2 to 2.1, replace
    elif isnan(evdata[gwname][0]) or \
            row['network_matched_filter_snr'] > evdata[gwname][0] or \
            'GWTC-2.1' in row['catalog.shortName'] and evdata[gwname][-1] == 'GWTC-2':
        evdata[gwname] = [
            row['network_matched_filter_snr'],
            row['chirp_mass'],
            row['chirp_mass_source'],
            row['redshift'],
            row['catalog.shortName']
        ]
        print('Replacing')
        print(gwname, evdata[gwname])

print(len(evdata), 'events found')

for gwname in evdata.copy():  # don't alter the thing being iterated over
    if isnan(evdata[gwname][1]):
        print('No detector frame mchirp for', gwname)
        if isnan(evdata[gwname][2]) or isnan(evdata[gwname][3]):
            print("Can't calculate it either!")
            print("Removing", gwname)
            del evdata[gwname]
            continue
        else:
            # mchirp_det = mchirp_src * (1 + z)
            evdata[gwname][1] = evdata[gwname][2] * (1 + evdata[gwname][3])
            print('Source frame', evdata[gwname][2], 'gives detector frame', evdata[gwname][1])
    else:
        print('Detector frame mchirp for', gwname)

print(len(evdata), 'events with mchirp det')

netsnr = []
mchirp_det = []
names = []
for name, properties in evdata.items():
    netsnr.append(properties[0])
    mchirp_det.append(properties[1])
    names.append(name)
netsnr = np.array(netsnr)
mchirp_det = np.array(mchirp_det)
snrsort = np.argsort(netsnr)
mcdsort = np.argsort(mchirp_det)
print('SNR ascending:')
for idx in snrsort:
    print(netsnr[idx], mchirp_det[idx], names[idx])
print('Detector frame mchirp ascending:')
for idx in mcdsort:
    print(mchirp_det[idx], names[idx], netsnr[idx])

thr = 10.
above = netsnr >= thr
bin_edges = [0.87, 1.74, 3.48, 6.96, 13.92, 27.84, 55.68, 111.36, 500.0]
print(np.histogram(mchirp_det[above], bins=bin_edges))

# See eg https://ldas-jobs.ligo.caltech.edu/~detchar/summary/1238166018-1259193618/range/
o3_bns_ranges = [135, 110, 50]
# Representative network range ~ quadrature sum over ifos
# 2.26 factor to convert to horizon
o3_net_horizon = 2.26 * sum([r ** 2 for r in o3_bns_ranges]) ** 0.5
print(o3_net_horizon)

