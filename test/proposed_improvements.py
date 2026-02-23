"""
Proposed improvements to prevent edge injection issues.

Key changes:
1. Better error messages for out-of-range wavelengths
2. Warning when wavelength seems to be in wrong units
3. Clearer documentation about units
"""

import numpy as np

def generate_gaussian_improved(fwhm, amplitude=None, center=None, array_length=None, 
                               area=None, wavelength=None, wavelengths=None):
    """
    Improved version with better unit detection and error messages.
    """
    
    # ... (previous validation code) ...
    
    # Handle wavelength-based center specification
    if wavelength is not None:
        if wavelengths is None:
            raise ValueError("'wavelengths' array must be provided when using 'wavelength' parameter")
        
        wavelengths = np.asarray(wavelengths)
        
        # Validate wavelengths array
        if len(wavelengths) != array_length:
            raise ValueError(f"wavelengths array length ({len(wavelengths)}) must match array_length ({array_length})")
        
        # Check if wavelengths array is monotonic
        if not (np.all(np.diff(wavelengths) > 0) or np.all(np.diff(wavelengths) < 0)):
            raise ValueError("wavelengths array must be monotonically increasing or decreasing")
        
        # Check if target wavelength is within range
        wave_min, wave_max = np.min(wavelengths), np.max(wavelengths)
        
        # IMPROVED: Detect likely unit mismatch
        if wavelength > wave_max:
            # Check if dividing by 1000 would put it in range (nm -> microns)
            if wave_min <= wavelength/1000 <= wave_max:
                raise ValueError(
                    f"wavelength {wavelength} is outside wavelengths range [{wave_min}, {wave_max}].\n"
                    f"HINT: Your wavelength looks like it's in NANOMETERS ({wavelength} nm),\n"
                    f"      but your wavelengths array appears to be in MICRONS ({wave_min}-{wave_max}).\n"
                    f"      Try using wavelength={wavelength/1000:.6f} instead (converted to microns)."
                )
            else:
                raise ValueError(
                    f"wavelength {wavelength} is outside wavelengths range [{wave_min}, {wave_max}].\n"
                    f"Check that both wavelength and wavelengths array use the same units (nm or microns)."
                )
        
        if wavelength < wave_min:
            # Check if multiplying by 1000 would put it in range (microns -> nm)
            if wave_min <= wavelength*1000 <= wave_max:
                raise ValueError(
                    f"wavelength {wavelength} is outside wavelengths range [{wave_min}, {wave_max}].\n"
                    f"HINT: Your wavelength looks like it's in MICRONS ({wavelength} μm),\n"
                    f"      but your wavelengths array appears to be in NANOMETERS ({wave_min}-{wave_max}).\n"
                    f"      Try using wavelength={wavelength*1000:.1f} instead (converted to nanometers)."
                )
            else:
                raise ValueError(
                    f"wavelength {wavelength} is outside wavelengths range [{wave_min}, {wave_max}].\n"
                    f"Check that both wavelength and wavelengths array use the same units (nm or microns)."
                )
        
        # Convert wavelength to fractional array index using interpolation
        indices = np.arange(array_length)
        center = float(np.interp(wavelength, wavelengths, indices))
    
    # Rest of function...
    pass


# Example improved error messages:

"""
OLD ERROR:
ValueError: wavelength 532 is outside wavelengths range [0.2, 1.995336705]

NEW ERROR:
ValueError: wavelength 532 is outside wavelengths range [0.2, 1.995336705].
HINT: Your wavelength looks like it's in NANOMETERS (532 nm),
      but your wavelengths array appears to be in MICRONS (0.2-1.995).
      Try using wavelength=0.532 instead (converted to microns).
"""
