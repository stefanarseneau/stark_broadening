## Asymmetric Stark Broadening in White Dwarf Photospheres
---

The Stark effect causes asymmetric broadening of the Balmer lines in white dwarf photospheres, which can induce a bias in radial velocity measurements. This repo is measuring the magnitude of that effect.


**Data Files**

The data for this project are divided into two parts: `data/processed/spy/` and `data/processed/arseneau2024/`. The former comes from the SPY survey by Napiwotski+2020, the latter comes from the mass-radius measurement from SDSS+APO of Arseneau+2024. Within the Arseneau+2024 dataset, files with the prefix `Hab` are fitted using only the first two Balmer lines, and files that have the prefixes `Habgd` are fitted with all four. Files with the suffix `_voigt` are fitted using a model spectrum consisting of two Voigt profiles per line, and files without are fit using the `corv` `1d_da_nlte` models.

**Script usage**

`python scripts/fitlines.py [dataset] [lines] [filepath] --from_cache`
* `dataset` must either be `sdss`, meaning to fit the SDSS+APO sample from Arseneau+2024, or `spy` meaning to fit the high-resolution SPY spectra from Napiwotski+2020.
* `lines` should be the list of Balmer lines to fit on (e.g. `abgd` for the first four lines)
* `filepath` is the path to the JSON file to write data to
* `--from_cache` should be used if you've already computed part of the dataset you're interested in and don't want to overwrite everything.

Example usage:

To fit the first four Balmer lines to the Arseneau+2024 sample while taking into account what has already been computed in the file `./data/processed/arseneau2024/Habgd_v2.json`, you would run the following command:

`python scripts/fitlines.py sdss abgd ./data/processed/arseneau2024/Habgd_v2.json --from_cache`

however since that file already contains data for all of those spectra, this will just terminate without actually computing anything.