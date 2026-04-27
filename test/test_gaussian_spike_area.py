# Test gaussian_spike_area on synthetic data.
#
# Constructs a spectrum with a known Gaussian spike, fits it, and checks that
# the recovered area is close to the true area (A * sigma * sqrt(2*pi)).

import sys
from pathlib import Path
cwd = Path(__file__).parent.resolve()
optical_seti_dir = cwd.parent.resolve()
sys.path.append(str(optical_seti_dir))

import numpy as np
import optical_seti_functions

# Build a synthetic spectrum: flat continuum + Gaussian spike
np.random.seed(42)
n = 1000
wave = np.linspace(5000.0, 5100.0, n)  # wavelengths in Angstroms

true_amplitude = 50.0
true_mean = 5050.0
true_stddev = 0.5   # Angstroms
true_area = true_amplitude * true_stddev * np.sqrt(2 * np.pi)

continuum = 100.0
flux = continuum + true_amplitude * np.exp(-0.5 * ((wave - true_mean) / true_stddev) ** 2)

# Identify the spike region (pixels within 3-sigma of the peak)
spike_mask = np.abs(wave - true_mean) < 3 * true_stddev
spike_indices = np.where(spike_mask)[0]
hits_start = spike_indices[0]
hits_end = spike_indices[-1]

area = optical_seti_functions.gaussian_spike_area(wave, flux, hits_start, hits_end)

print(f"True area:      {true_area:.4f}")
print(f"Recovered area: {area:.4f}")
print(f"Relative error: {abs(area - true_area) / true_area * 100:.2f}%")

assert abs(area - true_area) / true_area < 0.01, (
    f"Area {area:.4f} differs from true area {true_area:.4f} by more than 1%"
)
print("PASSED")
