import os
import json

from astropy.io import fits
from astroquery.sdss import SDSS
import matplotlib.pyplot as plt
from tqdm import tqdm

from . import measure

class SDSSHandler:
    def __init__(self, table, source_key, sdss5_path, apo_path):
        self.table = table
        self.source_key = source_key

        self.sdss5_path = sdss5_path
        self.apo_path = apo_path
        self.specfinder = {'sdss4' : self.fetch_sdss4, 'sdss5' : self.fetch_sdss5, 'apo' : self.fetch_apo}

    def analyze_table(self, outfile, n = 10, lines = ['a', 'b'], from_cache = False):
        if from_cache:
            with open(outfile) as json_file:
                results = json.load(json_file)
        else:
            results = {}

        for i, row in tqdm(self.table.iterrows(), total=self.table.shape[0]):
            if str(row.wd_source_id) not in results.keys():
                wavl, flux, ivar = self.specfinder[row.wd_rv_from](row)
                rv_data, windows = measure.test_windows(wavl, flux, ivar, n, lines)
                results[row.wd_source_id] = measure.process_results(rv_data, windows, plot=False)
                measure.write_dict_to_json(results, outfile)
        return results

    def fetch_sdss4(self, row):
        plate = row.wd_plate
        mjd = row.wd_mjd
        fiberid = row.wd_fiberid
        spec = SDSS.get_spectra(plate=plate, mjd=mjd, fiberID=fiberid)[0]

        wavl = 10**spec[1].data['LOGLAM']
        flux = spec[1].data['FLUX']
        ivar = spec[1].data['IVAR']
        return wavl, flux, ivar
    
    def fetch_sdss5(self, row):
        filepath = os.path.join(self.sdss5_path, row.wd_filepath)
        spec = fits.open(filepath)

        wavl = 10**spec[1].data['LOGLAM']
        flux = spec[1].data['FLUX']
        ivar = spec[1].data['IVAR']
        return wavl, flux, ivar
    
    def fetch_apo(self, row):
        filepath = os.path.join(self.apo_path, row.wd_filepath)
        spec = fits.open(filepath)

        wavl = spec[0].data[0,:]
        flux = spec[0].data[1,:]
        ivar = spec[0].data[2,:]
        return wavl, flux, ivar

