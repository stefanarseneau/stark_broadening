import sys
sys.path.append('../')

from stark import spy
from stark import sdss

import pandas as pd
import argparse

def fit_spy(path, lines, cache):
    obj = spy.fetch_objfile()
    sp = spy.SPYHandler(obj)
    result = sp.analyze_table(path, lines=lines, from_cache=cache)

def fit_sdss(path, lines, cache):
    catalog = pd.read_csv('../data/raw/Arseneau_2024.csv')
    catalog = catalog.drop_duplicates(subset=['wd_source_id'])
    subset = catalog.query("snr > 15 & R_chance_align < 0.1 & wd_rv_from != 'falcon'")

    sp = sdss.SDSSHandler(subset, 'wd_rv_from', '../data/raw/sdss5', '../data/raw/apo')
    result = sp.analyze_table(path, lines = lines, resolution = 1, from_cache=cache)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('dataset', nargs='?')
    parser.add_argument('lines', nargs='?')
    parser.add_argument('path', nargs='?')
    parser.add_argument('--from_cache', action='store_true')
    args = parser.parse_args()

    dataset = args.dataset
    lines = list(args.lines)
    path = args.path
    cache = args.from_cache

    datasets = {'sdss' : fit_sdss, 'spy' : fit_spy}
    assert dataset in datasets.keys(), f"dataset must be in {datasets.keys()}"

    datasets[dataset](path, lines, cache)
