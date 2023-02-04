## phase_velocity.py, a Python algorithm for calculating frequency-dependent phase velocity and radial variation of elastic waves in cylindrical bars

#### DESCRIPTION: 
The correlation between normalised frequency and phase velocity, m1 and m2, can be utilised to account for first-mode dispersion effect in pressure bar measurements using *process_SHPB.py* (see Van Lerberghe and Barr (2023)).

The open-source python algorithm *phase_velocity.py*, finds the first root of Bancroft’s (1941) equation using the bisection method, for a defined Poisson’s ratio, and over a defined range of normalised wavelength (d/L). The result is the normalised wave velocity, cp/c0, which corresponds to the first mode of propagation for longitudinal waves in an elastic cylindrical bar. Normalised wavelengths are also converted to normalised frequencies, fa/c0.

Normalised phase velocities are then used to calculate Tyas and Wilson’s (2001) factors m1 and m2, which account for wavelength dependent radial fluctuations in strain and Young’s modulus respectively.

The results m1, m2, norm_freqs and v_ratios are saved in 4 separate pickle files, in a folder titled *dispersion-factors*, for the corresponding Poisson’s ratio selected.

Both *process_SHPB.py* and *phase_velocity.py*, open-source Python algorithms, are inspired by Matlab scripts created by Barr (2016 & 2023), see below.

#### FILES INCLUDED:
-	*phase_velocity.py*: Includes the main python function *phase_velocity.py*, with the documentation on the use of the function included in the file as comments.
-	*phase_velocity.pdf*: An image showing the phase velocities, the factor m1 and normalised factor m2/E.

#### REFERENCES:
-	Bancroft, D. (1941) The Velocity of Longitudinal Waves in Cylindrical Bars. Physical Review, 59, 588-593.
-	Tyas, A., Wilson, A. J. (2001) An investigation of frequency domain dispersion correction of pressure bar signals. International Journal of Impact Engineering, 25, 87-101.

#### MATLAB SOFTWARE:
- Barr, A. D. (2016) dispersion.m - A MATLAB script for phase angle and amplitude correction of pressure bar signals. University of Sheffield.\
Software ORDA link: [https://doi.org/10.15131/shef.data.3996876.v1]
- Barr, A. D. (2023) phasevelocity.m - A MATLAB script to calculate the frequency-dependent phase velocity and
radial variation of elastic waves in cylindrical bars. University of Sheffield.\
Software ORDA link: [https://doi.org/10.15131/shef.data.21982604.v1]

#### PYTHON SOFTWARE:
- Van Lerberghe, A., Barr, A. D. (2023) *process_SHPB.py*, an Python algorithm for stress wave dispersion correction in split-Hopkinson pressure bar experiments. University of Sheffield.\
Software ORDA link: [https://doi.org/10.15131/shef.data.21973325] \
Software GitHub link: [https://github.com/ArthurVL-maker/Process_SHPB.git]

#### AUTHORS:
Arthur Van Lerberghe <avanlerberghe1@sheffield.ac.uk> & Andrew D. Barr <a.barr@sheffield.ac.uk>.
