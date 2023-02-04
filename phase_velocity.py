# phase_velocity.py
# ----------------------------------------------------------------------------------------------------------
# Finds the first root of Bancroft's equation using the bisection method, for a defined Poisson's ratio,
# and over a defined range of normalised wavelength (d/L). The result is the normalised phase velocity, c_p/c_0,
# which corresponds to the first mode of propagation for longitudinal waves in an elastic cylindrical bar.
# Normalised wavelengths are also converted to normalised frequency, fa/c_0.

# Normalised phase velocities are then used to calculate Tyas and Wilson's factors m1 and m2,
# which account for wavelength-dependent radial variations in strain and Young's modulus, respectively.

# The relationships between normalised frequency and phase velocity, m1 and m2 can be used to account
# for first-mode dispersion effects in pressure bar measurements, and are saved in the file format used by the
# open-source python algorithm 'Process_SHPB', referenced below (see Van Lerberghe and Barr (2023)).
# Both 'Process_SHPB' and 'Phase_velocity' algorithms, are inspired by MATLAB scripts created by Barr (2016 & 2023).

# INPUTS:
# - nu: Poisson's ratio of bar material used for split-Hopkinson pressure bar tests
# - l_ratios: Normalised wavelength range to calculate the first root of Bancroft's (1941) equation.

# OUTPUTS:
# - A folder titled 'dispersion_factors', with 4 .pickle files containing 'm1', 'm2', 'norm_freqs' and 'v_ratios'.

# REFERENCES:
# - Bancroft, D. (1941) The Velocity of Longitudinal Waves in Cylindrical Bars. Physical Review, 59, 588-593.
# - Tyas, A., Wilson, A. J. (2001) An investigation of frequency domain dispersion correction of pressure bar signals.
# International Journal of Impact Engineering, 25, 87-101.

# MATLAB SOFTWARE:
# - Barr, A. D. (2016) dispersion.m - A MATLAB script for phase angle and amplitude correction of pressure bar signals.
# University of Sheffield.
# Software ORDA link: (https://doi.org/10.15131/shef.data.3996876.v1)
	[Google](https://www.google.com) - _Google | Youtube | Gmail | Maps | PlayStore | GoogleDrive_

# - Barr, A. D. (2023) phasevelocity.m - A MATLAB script to calculate the frequency-dependent phase velocity and
# radial variation of elastic waves in cylindrical bars. University of Sheffield.
# Software ORDA [link](https://doi.org/10.15131/shef.data.21982604.v1)

# PYTHON SOFTWARE:
# - Van Lerberghe, A., Barr, A. D. (2023) Process_SHPB, an open-source python algorithm for stress wave dispersion
# correction in split-Hopkinson pressure bar experiments. University of Sheffield.
# Software ORDA link: (https://doi.org/10.15131/shef.data.21973325)
# Software GitHub link: [ADD GitHub Link]

# AUTHORS:
# Arthur Van Lerberghe (<avanlerberghe1@sheffield.ac.uk>) & Andrew D. Barr (<a.barr@sheffield.ac.uk>).
# ----------------------------------------------------------------------------------------------------------
# Imported modules:
import matplotlib.pyplot as plt
from scipy.special import jn
from pathlib import Path
import numpy as np
import warnings

# Disregard Warnings:
warnings.filterwarnings("ignore")

# ----------------------------------------------------------------
# VARIABLES
# ----------------------------------------------------------------
nu = 0.29  # Poisson's ratio of bar material.
l_ratios = np.arange(0, 2.01, 0.01).astype(complex)  # Normalised wavelength range to calculate.

# ----------------------------------------------------------------
# BANCROFT'S PHASE VELOCITY EQUATION: (see Bancroft, 1941)
# ----------------------------------------------------------------


def bancroft(v_ratio, l_ratio, nu):
    a = 1
    L = (2 * a) / l_ratio
    x = (v_ratio ** 2) * (1 + nu)
    beta = (1 - 2 * nu) / (1 - nu)
    h = (2 * np.pi / L) * np.sqrt(beta * x - 1)
    k = (2 * np.pi / L) * np.sqrt(2 * x - 1)
    Jha = h * a * jn(0, h * a) / jn(1, h * a)
    Jka = k * a * jn(0, k * a) / jn(1, k * a)

    return (x - 1) ** 2 * Jha - (beta * x - 1) * (x - Jka)  # out = 0 for valid roots


# ----------------------------------------------------------------
# CALCULATE PHASE VELOCITIES USING BISECTION METHOD:
# ----------------------------------------------------------------
v_ratios = np.ones(len(l_ratios)).astype(complex)  # Normalised phase velocity.

# For each normalised wavelength...
for i in range(1, l_ratios.shape[0]):
    l_ratio = l_ratios[i]  # Ratio of diameter to wavelength: d/L.
    high = v_ratios[i - 1]  # Initial guess interval high limit.
    low = 0  # Initial guess interval low limit.
    v_ratio = (low + high) / 2  # Initial v_ratio guess.
    tolerance = 1e-6  # Allowable tolerance.
    result = bancroft(v_ratio, l_ratio, nu)  # Initial result from Bancroft's equation.

    # ... use bisection method until root found (within tolerance):
    while abs(result) > tolerance:
        if result > 0:
            high = v_ratio
        else:
            low = v_ratio
        v_ratio = (low + high) / 2
        result = bancroft(v_ratio, l_ratio, nu)

    v_ratios[i] = v_ratio

norm_freqs = (l_ratios / 2) * v_ratios  # Normalised frequencies, f*a/c0.

# ----------------------------------------------------------------
# CALCULATE FACTORS m1 AND m1:  (see Tyas and Watson, 2001)
# ----------------------------------------------------------------
S = (1 - 2 * nu) / (1 - nu)
Z = (1 + nu) * (v_ratios ** 2)
ha = (l_ratios / 2) * (2 * np.pi) * np.sqrt(S * Z - 1)
ka = (l_ratios / 2) * (2 * np.pi) * np.sqrt(2 * Z - 1)
Jha = ha * jn(0, ha) / jn(1, ha)
Jka = ka * jn(0, ka) / jn(1, ka)

m1 = np.real((2 * (1 + (1 - S * Z) / (Z - 1))) / (Jha + ((1 - S * Z) / (Z - 1)) * Jka))  # Factor M_1
m1[0] = 1  # Above equation contains a division by zero at l_ratio = 0. Assume M_1(1) = 1.
m2 = v_ratios ** 2  # Normalised factor M_2 (M_2/E).

# ----------------------------------------------------------------
# PLOT RESULTS:
# ----------------------------------------------------------------
plt.Figure()
plt.rcParams.update({'font.family': 'Times'})

# Plot:
plt.xlabel('Normalised frequency, fa/c_0')
plt.title(f'\u03BD = {nu}')
plt.tick_params(which='both', direction='in', right=True, top=True)
plt.tick_params(axis='x', pad=8)
plt.plot(norm_freqs, v_ratios, label='$c_{p}$/$c_{0}$', linewidth=1.2)
plt.plot(norm_freqs, m1, label='m1', linewidth=1.2)
plt.plot(norm_freqs, m2, label='m2/E', linewidth=1.2)
plt.tight_layout(pad=1.5)
plt.legend(loc='upper right', edgecolor='black')
plt.ylim(-2, 3)
plt.xlim(0, 0.6)

# Save figure:
plt.savefig('phase_velocity.pdf')

# ----------------------------------------------------------------
# SAVE DATA:
# ----------------------------------------------------------------
print('-' * 20 + ' PROCESSING RESULTS ' + '-' * 20)
print(f'For \u03BD = {nu}')

# File format:
end = '.pickle'

# Folder to save data:
Path('dispersion_factors/').mkdir(parents=True, exist_ok=True)

# Save dispersion factors:
v_ratios_results = np.transpose(v_ratios.real)
norm_freqs_results = np.transpose(norm_freqs.real)
m1_results = np.transpose(m1.real)
m2_results = np.transpose(m2.real)

# Create designated files for dispersion factors:
filename_v_ratios_results = "v_ratios" + f'_\u03BD={nu}' + end
filename_norm_freqs_results = "norm_freqs" + f'_\u03BD={nu}' + end
filename_m1_results = "m1" + f'_\u03BD={nu}' + end
filename_m2_results = "m2" + f'_\u03BD={nu}' + end

# Create new filepaths for dispersion factors:
filepath_v_ratios_results = Path('dispersion_factors/' + filename_v_ratios_results)
filepath_norm_freqs_results = Path('dispersion_factors/' + filename_norm_freqs_results)
filepath_m1_results = Path('dispersion_factors/' + filename_m1_results)
filepath_m2_results = Path('dispersion_factors/' + filename_m2_results)

# Print new filename created:
print(f"New filename: {filename_v_ratios_results}")
print(f"New filename: {filename_norm_freqs_results}")
print(f"New filename: {filename_m1_results}")
print(f"New filename: {filename_m2_results}")

# Print new filepath created:
print(f"Filepath: {filepath_v_ratios_results}")
print(f"Filepath: {filepath_norm_freqs_results}")
print(f"Filepath: {filepath_m1_results}")
print(f"Filepath: {filepath_m2_results}")

# Save data:
v_ratios_results = np.savetxt(filepath_v_ratios_results, v_ratios_results, delimiter=',', fmt='%.10f')
norm_freqs_results = np.savetxt(filepath_norm_freqs_results, norm_freqs_results, delimiter=',', fmt='%.10f')
m1_results = np.savetxt(filepath_m1_results, m1_results, delimiter=',', fmt='%.10f')
m2_results = np.savetxt(filepath_m2_results, m2_results, delimiter=',', fmt='%.10f')

print('-' * 19 + ' PROCESSING COMPLETED ' + '-' * 19 + '\n')
