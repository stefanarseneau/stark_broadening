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

# Functions for processing the SPY data
# -------

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
        wl = air2vac(table['Table'].data)
        fl = table[':'].data
        mask = (5260 < wl) * (wl < 5280) # continuum region
        snr = np.nanmean(fl[mask]) / np.nanstd(fl[mask])
        ivar =  snr**2 / (table[':'].data + 1e-6)**2
        spec_table[file] = (wl, fl, ivar)
    return spec_table, obj

# Handle the analyzed data
# --------

def tabularize(data):
    # Create a DataFrame directly from the JSON
    df = pd.DataFrame([
        {'file': key, 
        'reference_window': values['windows'][0], 
        'reference_rv': values['rvs'][0],
        'reference_erv': values['e_rvs'][0],
        'reference_teff': np.nanmean(np.array(values['teff'])),
        'reference_logg': np.nanmean(np.array(values['logg'])),
        'windows': np.array(values['windows']), 
        'rvs': np.array(values['rvs']), 
        'e_rvs': np.array(values['e_rvs']), 
        'differentials': np.array(values['differentials']), 
        'significance': np.array(values['significance']),
        'redchisqr': np.array(values['redchisqr']), 
        'teff': np.array(values['teff']), 
        'e_teff': np.array(values['e_teff']), 
        'logg': np.array(values['logg']), 
        'e_logg': np.array(values['e_logg']), } 
        for key, values in data.items()
    ])
    return df

def apply_rules(data, sigma = 3 , rv_over_err = 3, chisqr = 10):
    mean = np.nanmean(np.concatenate(data.rvs.values))
    std = np.nanstd(np.concatenate(data.rvs.values))

    data['mask_arr'] = data.apply(lambda row: np.logical_and(np.logical_and(
                                np.abs(np.array(row['rvs']) / np.array(row['e_rvs'])) > rv_over_err,
                                (np.abs(np.array(row['rvs']) - mean) <= sigma * std)),
                                np.array(row['redchisqr'] < chisqr)), 
                        axis=1)     
    return data

def apply_mask(data, rows = ''):
    rows = [rows] if type(rows) != list else rows

    def fill_masked_with_nan(data, value_col):
        # Function to apply to each row
        def apply_mask(row):
            return np.where(row['mask_arr'], row[value_col], np.nan)
        # Apply the function to fill masked entries
        data[value_col] = data.apply(apply_mask, axis=1)
        return data[value_col]

    for row in rows:
        data[row] = fill_masked_with_nan(data, row)
    return data