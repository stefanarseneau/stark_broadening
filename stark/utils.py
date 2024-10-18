import numpy as np

def bin_temperatures(data, edges):
    temp_bn_to_range = lambda n: (edges[n-1], edges[n]) if (n != len(edges) and (n != 0)) else ((edges[n-1], np.infty) if n == len(edges) else (-np.infty, edges[n]))
    temp_bn_to_center = lambda n: np.mean(temp_bn_to_range(n))

    data['teff_bin'] = np.digitize(data['reference_teff'], edges)
    return data, temp_bn_to_range, temp_bn_to_center

def air2vac(wv):
    _tl=1.e4/np.array(wv)
    return (np.array(wv)*(1.+6.4328e-5+2.94981e-2/\
                          (146.-_tl**2)+2.5540e-4/(41.-_tl**2)))