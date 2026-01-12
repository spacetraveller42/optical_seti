# Test script for Gaussian generation and addition functions
# This tests the generate_gaussian and add_gaussian_to_array functions

import sys
from pathlib import Path
cwd = Path(__file__).parent.resolve() 
optical_seti_dir = cwd.parent.resolve() 
sys.path.append(str(optical_seti_dir))

import numpy as np
import matplotlib.pyplot as plt
import optical_seti_functions

# Test 1: Generate a simple Gaussian and verify its properties
print("Test 1: Generate Gaussian with FWHM=10, amplitude=100, center=50, length=100")
fwhm = 10.0
amplitude = 100.0
center = 50.0
array_length = 100

gaussian = optical_seti_functions.generate_gaussian(fwhm, amplitude, center, array_length)

# Verify peak is at the correct position and has the correct amplitude
peak_index = np.argmax(gaussian)
peak_value = np.max(gaussian)
print(f"Peak position: {peak_index} (expected: {center})")
print(f"Peak amplitude: {peak_value:.2f} (expected: {amplitude})")

# Verify FWHM by finding points where value = amplitude/2
half_max = amplitude / 2.0
indices_above_half = np.where(gaussian >= half_max)[0]
measured_fwhm = indices_above_half[-1] - indices_above_half[0]
print(f"Measured FWHM: {measured_fwhm:.2f} (expected: ~{fwhm})")

# Test 2: Add Gaussian to an existing array
print("\nTest 2: Add Gaussian to a baseline array")
baseline = np.ones(100) * 50.0  # Baseline of 50
data_with_gaussian = optical_seti_functions.add_gaussian_to_array(baseline, fwhm, amplitude, center)

# Verify the peak includes both baseline and gaussian
expected_peak = baseline[0] + amplitude
measured_peak = np.max(data_with_gaussian)
print(f"Peak with baseline: {measured_peak:.2f} (expected: {expected_peak:.2f})")

# Test 3: Create a plot showing the Gaussian is visible
print("\nTest 3: Create visualization plots")

# Create figure with three subplots
fig, axes = plt.subplots(3, 1, figsize=(10, 12))

# Subplot 1: Pure Gaussian
axes[0].plot(gaussian, 'b-', linewidth=2)
axes[0].axhline(y=amplitude/2, color='r', linestyle='--', label='Half Maximum')
axes[0].set_title(f'Generated Gaussian (FWHM={fwhm}, Amplitude={amplitude})')
axes[0].set_xlabel('Array Index')
axes[0].set_ylabel('Value')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# Subplot 2: Baseline data
axes[1].plot(baseline, 'g-', linewidth=2)
axes[1].set_title('Baseline Array (constant value)')
axes[1].set_xlabel('Array Index')
axes[1].set_ylabel('Value')
axes[1].grid(True, alpha=0.3)

# Subplot 3: Combined data showing Gaussian added to baseline
axes[2].plot(data_with_gaussian, 'r-', linewidth=2, label='Data + Gaussian')
axes[2].plot(baseline, 'g--', linewidth=1, alpha=0.5, label='Original Baseline')
axes[2].set_title('Gaussian Added to Baseline Array (clearly visible)')
axes[2].set_xlabel('Array Index')
axes[2].set_ylabel('Value')
axes[2].legend()
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('test_gaussian_generation.png', dpi=150)
print("Plot saved as 'test_gaussian_generation.png'")
plt.show()

# Test 4: Multiple Gaussians
print("\nTest 4: Add multiple Gaussians to array")
noisy_data = np.random.randn(200) * 5 + 100  # Noisy baseline around 100
result = optical_seti_functions.add_gaussian_to_array(noisy_data, 8.0, 80.0, 50)
result = optical_seti_functions.add_gaussian_to_array(result, 12.0, 60.0, 120)
result = optical_seti_functions.add_gaussian_to_array(result, 6.0, 100.0, 170)

plt.figure(figsize=(12, 6))
plt.plot(noisy_data, 'gray', alpha=0.5, linewidth=1, label='Original Noisy Data')
plt.plot(result, 'b-', linewidth=2, label='With Added Gaussians')
plt.title('Multiple Gaussians Added to Noisy Data')
plt.xlabel('Array Index')
plt.ylabel('Value')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig('test_multiple_gaussians.png', dpi=150)
print("Plot saved as 'test_multiple_gaussians.png'")
plt.show()

print("\nAll tests completed successfully!")
