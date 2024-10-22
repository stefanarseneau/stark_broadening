import corv
from tqdm import tqdm
import numpy as np

import matplotlib.pyplot as plt

def test_windows(wl, fl, ivar, n = 10, lines = ['a', 'b'], resolution = 0.0637, mask = True, plot_rvs = False):
    centres =  dict(a = 6564.61, b = 4862.68, g = 4341.68, d = 4102.89)
    window = dict(a = 15, b = 15, g = 15, d = 15)
    edges = dict(a = 0, b = 0, g = 0, d = 0)
    rv_data = []
    windows = []

    corvmodel = corv.models.WarwickDAModel(model_name='1d_da_nlte', names = lines, resolution = resolution, windows=window, edges=edges)
    ref_rv, ref_e_rv, ref_redchi, ref_param_res = corv.fit.fit_corv(wl, fl, ivar, corvmodel.model)
    rv_data.append([ref_rv, ref_e_rv, ref_param_res.redchi, ref_param_res.params['teff'].value, ref_param_res.params['teff'].stderr,
                        ref_param_res.params['logg'].value, ref_param_res.params['logg'].stderr])
    windows.append([window['a'], window['b']])

    if plot_rvs:
        corv.utils.lineplot(wl, fl, ivar, corvmodel.model, ref_param_res.params)

    if mask:
        mask_ha = ((centres['a'] - 3) < wl) * (wl < (centres['a'] + 3))
        mask_hb = ((centres['b'] - 3) < wl) * (wl < (centres['b'] + 3))
        mask = np.logical_or(mask_ha, mask_hb)
    else:
        mask = np.zeros(len(wl), dtype=bool)

    window['a'] += 15
    window['b'] += 15
    window['g'] += 15
    window['d'] += 15

    steps = np.linspace(0, 70, n)
    for step in steps:
        temp_window = window.copy()
        temp_window['a'] += step
        temp_window['b'] += step
        temp_window['g'] += step * 0.786 # to end at 85A after n steps
        temp_window['d'] += step * 0.643 # to end at 75A after n steps
        
        corvmodel = corv.models.WarwickDAModel(model_name='1d_da_nlte', names = lines, resolution = 0.0637, windows=temp_window, edges=edges)
        rv, e_rv, redchi, param_res = corv.fit.fit_corv(wl[~mask], fl[~mask], ivar[~mask], corvmodel.model)
        
        if plot_rvs:
            corv.utils.lineplot(wl[~mask], fl[~mask], ivar[~mask], corvmodel.model, param_res.params)
        rv_data.append([rv, e_rv, param_res.redchi, param_res.params['teff'].value, param_res.params['teff'].stderr,
                        param_res.params['logg'].value, param_res.params['logg'].stderr])
        windows.append([temp_window['a'], temp_window['b']])    
    return np.array(rv_data), np.array(windows)

def process_results(parameters, windows, plot=True, **kwargs):
    rvs, e_rvs, chisqr, teff, e_teff, logg, e_logg = parameters[:,0], parameters[:,1], parameters[:,2], parameters[:,3], parameters[:,4], parameters[:,5], parameters[:,6]
    windows = windows[:,0]
    diff = rvs - rvs[0]
    significance = lambda rv, e_rv : (rv - rvs[0]) / np.sqrt(e_rv**2 + e_rvs[0]**2)
    sigs = [significance(r, e) for r, e in zip(rvs, e_rvs)]

    if plot:
        fig, (ax1, ax2) = plt.subplots(nrows = 2, figsize=(10,10), **kwargs)
        ax1.errorbar(windows, diff, yerr = e_rvs, fmt='o', c = 'k', markersize=8, capsize=10)
        ax1.axhline(y=0, c='k', ls='--')
        ax1.set_ylabel('$RV - RV_0$ [$km/s]$')

        ax2.scatter(windows, sigs, s = 40, c = 'k')
        ax2.axhline(y=0, c='k', ls='--')
        
        ymin, ymax = ax2.get_ylim()
        if (ymax < 2):
            ax2.set_ylim(ymin,2)
        ymin, ymax = ax2.get_ylim()
        if (ymin > -2):
            ax2.set_ylim(-2, ymax)

        xmin, xmax = ax2.get_xlim()
        ax2.fill_between([xmin,xmax], -1, 1, color='green', alpha=0.5, zorder=0)
        ax2.set_xlim(xmin, xmax)
        ax2.set_ylabel('$\\chi$ Agreement')
        ax2.set_xlabel(r'Fit Window Size [$\AA$]')

    res_dict = {
        'windows' : windows,
        'rvs' : rvs,
        'e_rvs' : e_rvs,
        'differentials' : diff,
        'significance' : sigs,
        'redchisqr' : chisqr,
        'teff' : teff,
        'e_teff' : e_teff,
        'logg' : logg,
        'e_logg' : e_logg,
    }
    return res_dict

import json

def write_dict_to_json(dictionary, filename):
    def convert_ndarray(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        
    with open(filename, 'w') as json_file:
        json.dump(dictionary, json_file, indent=4, default=convert_ndarray)

def analyze_spectra(spec_dict, outfile, lines = ['a','b'], from_cache=False):
    if from_cache:
        with open(outfile) as json_file:
            results = json.load(json_file)
    else:
        results = {}
    for key in spec_dict.keys():
        if from_cache and (key in list(results.keys())):
            continue
        print(key)
        wl, fl, ivar = spec_dict[key]
        try:
            parameters, windows = test_windows(wl, fl, ivar, lines = lines)
        except:
            print('Fit failed!')
            pass
        results[key] = process_results(parameters, windows, plot=False)
        write_dict_to_json(results, outfile)
    return results

