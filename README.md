## Phase_velocity, an open-source python algorithm for calculating frequency-dependent phase velocity and radial variation of elastic waves in cylindrical bars.

#### DESCRIPTION: 
The correlation between normalised frequency and phase velocity, m1 and m2, can be utilised to account for first-mode dispersion effect in pressure bar measurements using Process_SHPB, open-source python algorithm. The Software is available to download on ORDA [link] & GitHub [link].

The open-source python algorithm Phase_velocity, finds the first root of Bancroft’s (1941) equation using the bisection method, for a defined Poisson’s ratio, and over a defined range of normalised wavelength (d/L). The result is the normalised wave velocity, cp/c0, which corresponds to the first mode of propagation for longitudinal waves in an elastic cylindrical bar. Normalised wavelengths are also converted to normalised frequencies, fa/c0.

Normalised phase velocities are then used to calculate Tyas and Wilson’s (2001) factors m1 and m2, which account for wavelength dependent radial fluctuations in strain and Young’s modulus respectively.

The results m1, m2, norm_freqs and v_ratios are saved in 4 separate pickle files, in a folder titled dispersion-factors, for the corresponding Poisson’s ratio selected.

#### FILES INCLUDED:
-	phase_velocity.py: Includes the main python function phase_velocity.py, with the documentation on the use of the function included in the file as comments.
-	phase_velocity.pdf: An image showing the phase velocities, the factor m1 and normalised factor m2/E.

#### REFERENCES:
-	Bancroft, D. (1941) The Velocity of Longitudinal Waves in Cylindrical Bars. Physical Review, 59, 588-593.
-	Tyas, A., Wilson, A. J. (2001) An investigation of frequency domain dispersion correction of pressure bar signals. International Journal of Impact Engineering, 25, 87-101.

#### MATLAB SOFTWARE:

#### PYTHON SOFTWARE:

#### AUTHORS:
Arthur Van Lerberghe <avanlerberghe1@sheffield.ac.uk> & Andrew D. Barr <a.barr@sheffield.ac.uk>.
