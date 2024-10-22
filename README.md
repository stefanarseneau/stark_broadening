## Asymmetric Stark Broadening in White Dwarf Photospheres
---

The Stark effect causes asymmetric broadening of the Balmer lines in white dwarf photospheres, which can induce a bias in radial velocity measurements. This repo is measuring the magnitude of that effect.

Script usage:
`python scripts/fitlines.py [dataset] [lines] [filepath] --from_cache`
* `dataset` must either be `sdss`, meaning to fit the SDSS+APO sample from Arseneau+2024, or `spy` meaning to fit the high-resolution SPY spectra from Napiwotski+2020.
* `lines` should be the list of Balmer lines to fit on (e.g. `abgd` for the first four lines)
* `filepath` is the path to the JSON file to write data to
* `--from_cache` should be used if you've already computed part of the dataset you're interested in and don't want to overwrite everything.

Example usage:
`python scripts/fitlines.py sdss abgd ./data/processed/arseneau2024/Habgd_v2.json --from_cache`
