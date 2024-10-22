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
import json

from . import utils
from . import measure

# Functions for processing the SPY data
# -------

class SPYHandler:
    def __init__(self, table, specpath = '../data/raw/sp'):
        self.table = table
        self.specpath = specpath

    def analyze_table(self, outfile, n = 10, lines = ['a', 'b'], resolution = 0.0637, from_cache = False):
        if from_cache:
            with open(outfile) as json_file:
                results = json.load(json_file)
        else:
            results = {}

        for i, row in tqdm(self.table.iterrows(), total=self.table.shape[0]):
            if str(row.FileName) not in results.keys():
                try:
                    wavl, flux, ivar = self.read_spectrum(row.FileName)
                    rv_data, windows = measure.test_windows(wavl, flux, ivar, n, lines, resolution, mask=False)
                    results[row.FileName] = measure.process_results(rv_data, windows, plot=False)
                    measure.write_dict_to_json(results, outfile)
                except:
                    print(f'Fit Failed: {row.FileName}')
        return results

    def read_spectrum(self, file, download_files = False):
        # find, download, or skip the file
        path = os.path.join(self.specpath, file)
        #if not os.path.isfile(path):
        #    print(f'Err!! Cannot find file {file}')   
        # read data
        table = ascii.read(path)
        table_meta = table.meta['comments']
        find_index = lambda string : list(filter(lambda x: string in x, table_meta))
        #table['ra'][i] = float(re.findall(r'\d+\.\d+', find_index('rekta')[0])[0])
        #table['de'][i] = float(re.findall(r'\d+\.\d+', find_index('dekli')[0])[0])
        #table['decimal_date'][i] = float(re.findall(r'\d+\.\d+', find_index('norm_date')[0])[0])
        # calculate heliocentric correction
        #t = Time(self.table['decimal_date'][i], format='decimalyear')
        #sc = SkyCoord(self.table['ra'][i]*u.deg, self.table['de'][i]*u.deg)
        #loc = EarthLocation.of_site('lasilla')
        #self.table['helio_corr'][i] = sc.radial_velocity_correction(kind='heliocentric', obstime=t, location=loc).to(u.km/u.s).value
        # pull the spectrum elements
        wl = utils.air2vac(table['Table'].data)
        fl = table[':'].data
        mask = (5260 < wl) * (wl < 5280) # continuum region
        snr = np.nanmean(fl[mask]) / np.nanstd(fl[mask])
        snr = snr if ~np.isnan(snr) else 1
        ivar =  snr**2 / (table[':'].data + 1e-6)**2
        return wl, fl, ivar

def fetch_objfile(path = 'http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/fits.gz?J/A+A/638/A131/objects.dat'):
    objects = Table.read(path)
    objects['FileName'] = [re.sub(r'dat', 'dat.gz', s).replace(' ', '') for s in objects['FileName']]
    objects['Name'] = [s.replace(' ', '') for s in objects['Name']]
    objects['ra'] = np.ones(len(objects)) * np.nan
    objects['de'] = np.ones(len(objects)) * np.nan
    objects['decimal_date'] = np.ones(len(objects)) * np.nan
    objects['helio_corr'] = np.ones(len(objects)) * np.nan
    return objects.to_pandas()