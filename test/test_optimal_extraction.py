# Test to verify that optimal_extraction_harps correctly extracts spectra
# and filters out cosmic rays while preserving real spectral signal.

import sys
from pathlib import Path
cwd = Path(__file__).parent.resolve()
optical_seti_dir = cwd.parent.resolve()
sys.path.append(str(optical_seti_dir))

import numpy as np
import optical_seti_functions

# Create a synthetic CCD image with one spectral order
nrow, ncol = 100, 500
img = np.zeros((nrow, ncol), dtype=float)

# Add a spectral order: a horizontal band centered at row 50
# with a Gaussian spatial profile (simulating the telescope PSF / slit function)
order_center = 50
order_width = 5
for col in range(ncol):
    for row in range(order_center - order_width, order_center + order_width + 1):
        distance = abs(row - order_center)
        img[row, col] = 1000.0 * np.exp(-0.5 * (distance / 2.0)**2)

# Add a cosmic ray hit at a single pixel
cosmic_ray_col = 250
cosmic_ray_row = order_center
img[cosmic_ray_row, cosmic_ray_col] += 50000.0

# Define order traces (polynomial coefficients, highest degree first)
order_traces = np.array([
    [0, 0, 0, 0, order_center - 20],  # extension order below
    [0, 0, 0, 0, order_center],        # real order
    [0, 0, 0, 0, order_center + 20],  # extension order above
])

# Run extraction
spectrum, uncertainties, slitfunction = optical_seti_functions.optimal_extraction_harps(
    img, order_traces=order_traces, extraction_width=order_width
)

# Verify output shapes
assert spectrum.shape == (3, ncol), "Spectrum shape mismatch: expected (3, {}), got {}".format(ncol, spectrum.shape)
assert uncertainties.shape == (3, ncol), "Uncertainties shape mismatch"

# Verify cosmic ray rejection: the extracted value at the cosmic ray location
# should be close to the value at neighboring columns (within 50%)
real_spec = spectrum[1]
ratio = real_spec[cosmic_ray_col] / real_spec[cosmic_ray_col + 10]
assert abs(ratio - 1) < 0.5, "Cosmic ray not rejected: ratio = {}".format(ratio)

print("All tests passed.")
print("  Spectrum shape:", spectrum.shape)
print("  Cosmic ray rejection ratio:", ratio)
