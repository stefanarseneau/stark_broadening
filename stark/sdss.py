import os
import json

from astropy.io import fits
from astroquery.sdss import SDSS
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np

from . import measure
from . import interpolator as int

class SDSSHandler:
    def __init__(self, table, source_key, sdss5_path, apo_path):
        self.table = table
        self.source_key = source_key

        self.sdss5_path = sdss5_path
        self.apo_path = apo_path
        self.specfinder = {'sdss4' : self.fetch_sdss4, 'sdss5' : self.fetch_sdss5, 'apo' : self.fetch_apo}

    def analyze_table(self, outfile, n = 10, lines = ['a', 'b'], resolution = 1, from_cache = False):
        if from_cache:
            with open(outfile) as json_file:
                results = json.load(json_file)
        else:
            results = {}

        for i, row in tqdm(self.table.iterrows(), total=self.table.shape[0]):
            if str(row.wd_source_id) not in results.keys():
                wavl, flux, ivar = self.specfinder[row.wd_rv_from](row)
                rv_data, windows = measure.test_windows(wavl, flux, ivar, n, lines, resolution, mask=False)
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
    
def correct_gband(bp, rp, astrometric_params_solved, phot_g_mean_mag):
    bp_rp = bp - rp

    if np.isscalar(bp_rp) or np.isscalar(astrometric_params_solved) or np.isscalar(phot_g_mean_mag):
        bp_rp = np.float64(bp_rp)
        astrometric_params_solved = np.int64(astrometric_params_solved)
        phot_g_mean_mag = np.float64(phot_g_mean_mag)
    
    if not (bp_rp.shape == astrometric_params_solved.shape == phot_g_mean_mag.shape):
        raise ValueError('Function parameters must be of the same shape!')
    
    do_not_correct = np.isnan(bp_rp) | (phot_g_mean_mag<13) | (astrometric_params_solved == 31)
    bright_correct = np.logical_not(do_not_correct) & (phot_g_mean_mag>=13) & (phot_g_mean_mag<=16)
    faint_correct = np.logical_not(do_not_correct) & (phot_g_mean_mag>16)
    bp_rp_c = np.clip(bp_rp, 0.25, 3.0)
    
    correction_factor = np.ones_like(phot_g_mean_mag)
    correction_factor[faint_correct] = 1.00525 - 0.02323*bp_rp_c[faint_correct] + \
        0.01740*np.power(bp_rp_c[faint_correct],2) - 0.00253*np.power(bp_rp_c[faint_correct],3)
    correction_factor[bright_correct] = 1.00876 - 0.02540*bp_rp_c[bright_correct] + \
        0.01747*np.power(bp_rp_c[bright_correct],2) - 0.00277*np.power(bp_rp_c[bright_correct],3)
    
    gmag_corrected = phot_g_mean_mag - 2.5*np.log10(correction_factor)
    return gmag_corrected

def fit_radius(row, interpolator, **kwargs):
    g_mag, bp_mag, rp_mag = row.wd_phot_g_mean_mag, row.wd_phot_bp_mean_mag, row.wd_phot_rp_mean_mag
    fluxes = np.array([row.wd_phot_g_mean_flux, row.wd_phot_bp_mean_flux, row.wd_phot_rp_mean_flux])
    e_fluxes = np.array([row.wd_phot_g_mean_flux_error, row.wd_phot_bp_mean_flux_error, row.wd_phot_rp_mean_flux_error])
    e_gmag, e_bpmag, e_rpmag = fluxes / (1.09 * e_fluxes)

    astrometric_params = row.wd_astrometric_params_solved
    g_mag = correct_gband(bp_mag, rp_mag, astrometric_params, g_mag)

    obs_mag = np.array([g_mag, bp_mag, rp_mag])
    e_obs_mag = np.array([e_gmag, e_bpmag, e_rpmag])
    distance = 1000 / row.ms_parallax

    bands = ['Gaia_G', 'Gaia_BP', 'Gaia_RP']
    interp = interpolator(bands, **kwargs)
    engine = int.CoarseEngine(interp)
    return engine(obs_mag, e_obs_mag, distance, p0=[10000, 8, 0.01])


