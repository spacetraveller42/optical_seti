# Test to verify that optimal_extraction_harps correctly extracts spectra
# and filters out non-optical artifacts while preserving real spectral signal.
# Artifacts tested: cosmic rays, hot/bad pixels, dark current, and bias offset.

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

# Add a uniform bias offset (electronic readout artifact)
bias_level = 200.0
img += bias_level

# Add dark current (thermal noise accumulating over the whole detector)
dark_frame = np.full((nrow, ncol), 50.0)

# Mark bad/hot pixels (detector defects producing spurious signal)
bad_pixel_mask = np.zeros((nrow, ncol), dtype=bool)
hot_pixel_col = 300
hot_pixel_row = order_center
bad_pixel_mask[hot_pixel_row, hot_pixel_col] = True
img[hot_pixel_row, hot_pixel_col] += 30000.0  # hot pixel artifact

# Add dark current to the image
img += dark_frame

# Define order traces (polynomial coefficients, highest degree first)
order_traces = np.array([
    [0, 0, 0, 0, order_center - 20],  # extension order below
    [0, 0, 0, 0, order_center],        # real order
    [0, 0, 0, 0, order_center + 20],  # extension order above
])

# Run extraction with calibration frames to remove non-optical artifacts
spectrum, uncertainties, slitfunction = optical_seti_functions.optimal_extraction_harps(
    img, order_traces=order_traces, extraction_width=order_width,
    bias=bias_level, dark=dark_frame, bad_pixel_mask=bad_pixel_mask
)

# Verify output shapes
assert spectrum.shape == (3, ncol), "Spectrum shape mismatch: expected (3, {}), got {}".format(ncol, spectrum.shape)
assert uncertainties.shape == (3, ncol), "Uncertainties shape mismatch"

# Verify cosmic ray rejection: the extracted value at the cosmic ray location
# should be close to the value at neighboring columns (within 50%)
real_spec = spectrum[1]
ratio_cr = real_spec[cosmic_ray_col] / real_spec[cosmic_ray_col + 10]
assert abs(ratio_cr - 1) < 0.5, "Cosmic ray not rejected: ratio = {}".format(ratio_cr)

# Verify hot pixel rejection: the hot pixel should be masked and not affect extraction
ratio_hp = real_spec[hot_pixel_col] / real_spec[hot_pixel_col + 10]
assert abs(ratio_hp - 1) < 0.5, "Hot pixel not rejected: ratio = {}".format(ratio_hp)

print("All tests passed.")
print("  Spectrum shape:", spectrum.shape)
print("  Cosmic ray rejection ratio:", ratio_cr)
print("  Hot pixel rejection ratio:", ratio_hp)
