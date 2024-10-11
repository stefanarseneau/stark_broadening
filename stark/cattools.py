from astropy.io import ascii, fits
from astropy.table import Table
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.constants import c
import astropy.units as u

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

import hashlib
import urllib.request
import os
import re

def air2vac(wv):
    _tl=1.e4/np.array(wv)
    return (np.array(wv)*(1.+6.4328e-5+2.94981e-2/\
                          (146.-_tl**2)+2.5540e-4/(41.-_tl**2)))

def fetch_objfile():
    objects = Table.read('http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/fits.gz?J/A+A/638/A131/objects.dat')
    objects['FileName'] = [re.sub(r'dat', 'dat.gz', s).replace(' ', '') for s in objects['FileName']]
    objects['Name'] = [s.replace(' ', '') for s in objects['Name']]
    objects['ra'] = np.ones(len(objects)) * np.nan
    objects['de'] = np.ones(len(objects)) * np.nan
    objects['decimal_date'] = np.ones(len(objects)) * np.nan
    objects['helio_corr'] = np.ones(len(objects)) * np.nan
    return objects

def read_spectrum(objects, dirpath, download_files = False):
    dirpath = os.path.abspath(dirpath)
    info_table, spec_table = {}, {}

    for i, file in enumerate(tqdm(objects['FileName'])):
        # find, download, or skip the file
        path = os.path.join(dirpath, file)
        if not os.path.isfile(path):
            print(f'Err!! Cannot find file {file}')
            continue
        # read data
        table = ascii.read(path)
        table_meta = table.meta['comments']
        find_index = lambda string : list(filter(lambda x: string in x, table_meta))
        objects['ra'][i] = float(re.findall(r'\d+\.\d+', find_index('rekta')[0])[0])
        objects['de'][i] = float(re.findall(r'\d+\.\d+', find_index('dekli')[0])[0])
        objects['decimal_date'][i] = float(re.findall(r'\d+\.\d+', find_index('norm_date')[0])[0])
        # calculate heliocentric correction
        t = Time(objects['decimal_date'][i], format='decimalyear')
        sc = SkyCoord(objects['ra'][i]*u.deg, objects['de'][i]*u.deg)
        loc = EarthLocation.of_site('lasilla')
        objects['helio_corr'][i] = sc.radial_velocity_correction(kind='heliocentric', obstime=t, location=loc).to(u.km/u.s).value
        # pull the spectrum elements
        snr = 20
        wl = air2vac(table['Table'].data)
        #wl = table['Table'].data
        fl = table[':'].data
        ivar =  snr**2 / (table[':'].data + 1e-6)**2
        spec_table[file] = (wl, fl, ivar)
    return spec_table, objects

def plot_ensemble(result_dict):
    fig, (ax1, ax2) = plt.subplots(nrows = 2, figsize=(10,10))

    for key in result_dict.keys():
        ax1.errorbar(result_dict[key]['windows'], result_dict[key]['differentials'], 
                        yerr = result_dict[key]['e_rvs'], markersize=8, capsize=10)
        ax2.scatter(result_dict[key]['windows'], result_dict[key]['significance'])
    
    ax1.axhline(y = 0, c = 'k', ls = '--')
    ax1.set_ylabel('$RV - RV_0$ [$km/s]$')

    ax2.axhline(y = 0, c = 'k', ls = '--')
    ax2.set_ylabel('$\\chi$ Agreement')
    ax2.set_xlabel(r'Fit Window Size [$\AA$]')

    fig, ax = plt.subplots(nrows = 1, figsize=(10,6))

    rvs, wls = [], []
    for key in result_dict.keys():
        wls.append(result_dict[key]['windows'])
        rvs.append(result_dict[key]['differentials'])
    rvs = np.array(rvs, dtype=float)

    wl = np.array(wls)[0,:]
    mean_rv = np.nanmean(rvs, axis=0)
    median_rv = np.nanmedian(rvs, axis=0)
    std_rv = np.nanstd(rvs, axis=0)

    scatter = {'s' : 40}
    #ax.errorbar(wl, mean_rv, yerr = std_rv, label = f'Mean (n = {len(result_dict.keys())})', c = 'k', fmt = 'o', markersize=8, capsize=10)
    ax.scatter(wl, mean_rv, label = f'Mean (n = {len(result_dict.keys())})', c = 'k', s = 60, zorder=10)
    ax.scatter(wl, median_rv, label = f'Median (n = {len(result_dict.keys())})', c = 'blue', s = 60, zorder=10)

    xmin, xmax = ax.get_xlim()
    ax.fill_between([xmin,xmax], -1, 1, color='green', alpha=0.5, zorder=0)
    ax.set_xlim(xmin, xmax)
    ax.axhline(y = 0, c = 'k', ls = '--')
    
    ax.set_xlabel(r'Fit Window Size [$\AA$]')
    ax.set_ylabel('$v_r - v_{r0}$ $[km /s]$')
    ax.legend(framealpha=0)