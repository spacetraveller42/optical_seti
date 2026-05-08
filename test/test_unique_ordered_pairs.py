# Test for unique_ordered_pairs() in seti_catalog_functions.
#
# This test verifies that unique_ordered_pairs() correctly removes duplicate
# values while preserving the original order and the correspondence between
# the two lists.

import sys
from pathlib import Path
cwd = Path(__file__).parent.resolve()
optical_seti_dir = cwd.parent.resolve()
sys.path.append(str(optical_seti_dir))

import seti_catalog_functions

# --- Test 1: Basic deduplication ---
# Two parallel lists that mirror each other: every repeated star name also
# has a repeated spectral type, so after deduplication both lists are clean.
stars         = ["HD1234",  "GJ551",  "HD1234",  "Tau Ceti", "GJ551",  "Proxima"]
spectral_types = ["G2V",    "M5.5Ve", "G2V",     "G8V",      "M5.5Ve", "M5.5Ve"]

unique_stars, unique_types = seti_catalog_functions.unique_ordered_pairs(stars, spectral_types)

print("=== Test 1: Basic deduplication ===")
print("Input stars:          ", stars)
print("Input spectral types: ", spectral_types)
print("Unique stars:         ", unique_stars)
print("Unique spectral types:", unique_types)

# Order must be preserved (first occurrences only)
assert unique_stars  == ["HD1234", "GJ551", "Tau Ceti", "Proxima"], \
    f"Wrong unique_stars: {unique_stars}"
# Correspondence must be maintained (index 0 of unique_stars pairs with index 0 of unique_types)
assert unique_types == ["G2V", "M5.5Ve", "G8V", "M5.5Ve"], \
    f"Wrong unique_types: {unique_types}"

print("PASSED: order and correspondence are correct.\n")

# --- Test 2: No duplicates in input ---
# If there are no duplicates the output must equal the input exactly.
stars2          = ["Alpha Cen A", "Alpha Cen B", "Barnard's Star"]
distances2      = [1.34,          1.34,           1.83]

unique_stars2, unique_dist2 = seti_catalog_functions.unique_ordered_pairs(stars2, distances2)

print("=== Test 2: No duplicates ===")
print("Input stars:    ", stars2)
print("Unique stars:   ", unique_stars2)
assert unique_stars2 == stars2,     f"unique_stars2 changed: {unique_stars2}"
assert unique_dist2  == distances2, f"unique_dist2  changed: {unique_dist2}"
print("PASSED: unchanged when no duplicates.\n")

# --- Test 3: All duplicates ---
# A list where every entry is the same: output should have exactly one element.
stars3  = ["Sirius", "Sirius", "Sirius"]
temps3  = [9940,     9940,     9940]

unique_stars3, unique_temps3 = seti_catalog_functions.unique_ordered_pairs(stars3, temps3)

print("=== Test 3: All duplicates ===")
print("Input stars:  ", stars3)
print("Unique stars: ", unique_stars3)
assert unique_stars3 == ["Sirius"], f"unique_stars3: {unique_stars3}"
assert unique_temps3 == [9940],     f"unique_temps3: {unique_temps3}"
print("PASSED: single element when all entries are the same.\n")

# --- Test 4: Realistic CSV-style data (tab-separated, matching OSETI_targets.txt) ---
# Simulate what parse_star_list() returns for a small catalog where one star
# is repeated (multiple observations in the predownloader list).
target_list   = ["HD127423", "BD+053829", "BD+053829", "BD+072351"]
spectral_list = ["G0V",      "G5",        "G5",        "G5"]

unique_targets, unique_spectra = seti_catalog_functions.unique_ordered_pairs(target_list, spectral_list)

print("=== Test 4: Realistic catalog data ===")
for star, sp in zip(unique_targets, unique_spectra):
    print(f"  {star} -> {sp}")

assert unique_targets == ["HD127423", "BD+053829", "BD+072351"]
assert unique_spectra == ["G0V",      "G5",        "G5"]
print("PASSED: duplicate star removed, remaining stars and types intact.\n")

print("All tests passed.")
