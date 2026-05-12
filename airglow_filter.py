# Airglow filter script.
# For each star in the search list, checks whether the Doppler-shifted
# central wavelength of the observed feature falls within a known airglow
# emission line and flags it as a likely Earth-atmosphere false positive.
#
# Bugs fixed vs. the original version:
#   1. Added missing doppler() and find_closest_index_numpy() to
#      optical_seti_functions.py and import them here.
#   2. Moved airglow_hits list initialisation BEFORE the loop so hits
#      accumulate across all stars instead of being reset every iteration.
#   3. Added .strip() to all values read from text files so that trailing
#      newlines / whitespace do not cause float() conversion failures.
#   4. Replaced bare fits.open() with a context manager (with …) so the
#      file handle is always closed even if an exception is raised.
#   5. Removed the unused earth_shifted_wavelengths list.

from astropy.io import fits
import numpy as np
from optical_seti_functions import doppler, find_closest_index_numpy

# ── Load search list ──────────────────────────────────────────────────────────
search_file = "/Users/blfields/Desktop/LaserInjector/thebigsearchstage2.txt"
rows = [line.strip() for line in open(search_file).readlines()]

stars          = [x.split(',')[0] for x in rows]
spectral_types = [x.split(',')[1] for x in rows]
files          = [x.split(',')[2] for x in rows]
lambstarts     = [x.split(',')[5] for x in rows]
lambends       = [x.split(',')[6] for x in rows]

# ── Load airglow reference table ──────────────────────────────────────────────
airglow_file = "/Users/blfields/Desktop/LaserInjector/airglow_table_comma_delimited.txt"
airglow_rows = [line.strip() for line in open(airglow_file).readlines()]

airglow_lines = np.array([x.split(',')[2] for x in airglow_rows]).astype(float)
airglow_fwhm  = np.array([x.split(',')[3] for x in airglow_rows]).astype(float)

# ── Per-star airglow check ────────────────────────────────────────────────────
filter_multiplier = 1000

# Initialise BEFORE the loop so hits accumulate across all stars.
airglow_hits = []

for file, lamb1, lamb2 in zip(files, lambstarts, lambends):
    central_wavelength = (float(lamb1) + float(lamb2)) / 2

    with fits.open(file) as specfits:
        berv = specfits[0].header["HIERARCH ESO DRS BERV"]

    earth_shifted_wavelength = doppler(central_wavelength, berv)

    idx              = find_closest_index_numpy(earth_shifted_wavelength, airglow_lines)
    fwhm             = airglow_fwhm[idx]
    airglow_difference = abs(airglow_lines[idx] - earth_shifted_wavelength)

    print(earth_shifted_wavelength)
    print(airglow_lines[idx])
    print(airglow_difference)
    print(fwhm)

    if airglow_difference < (filter_multiplier * fwhm):
        airglow_hits.append(earth_shifted_wavelength)

# ── Summary ───────────────────────────────────────────────────────────────────
print(airglow_hits)
print(len(airglow_hits))
print(len(stars) - len(airglow_hits))
