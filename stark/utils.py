import numpy as np
import pandas as pd

def bin_temperatures(data, edges):
    temp_bn_to_range = lambda n: (edges[n-1], edges[n]) if (n != len(edges) and (n != 0)) else ((edges[n-1], np.infty) if n == len(edges) else (-np.infty, edges[n]))
    temp_bn_to_center = lambda n: np.mean(temp_bn_to_range(n))

    data['teff_bin'] = np.digitize(data['reference_teff'], edges)
    return data, temp_bn_to_range, temp_bn_to_center

def air2vac(wv):
    _tl=1.e4/np.array(wv)
    return (np.array(wv)*(1.+6.4328e-5+2.94981e-2/\
                          (146.-_tl**2)+2.5540e-4/(41.-_tl**2)))

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