/*
Copyright 2024 Padraic Shafer.

Adapted from the work of Joe Dvorak, which is based on the published article
by R. Reininger et al., NIM-A 538, 760-770 (2005).
*/

#ifndef VIA_H
#define VIA_H


/// Print verbose debug output from VIA-related functions
#  ifndef VIA_DEBUG
#    define VIA_DEBUG (0)
#  else
#    ifdef DEBUG
#      define VIA_DEBUG (DEBUG)
#    else
#      define VIA_DEBUG (0)
#    endif
#  endif


/// Photon energy-wavelength product, eV * nm
#define eVnm_PRODUCT (1239.84197)

/// Convert distance units, mm per nm
#define mm_PER_nm (1e-6)


/// @brief Calculates the grating incidence angle for the specified diffraction order
/// @param energy photon energy, in eV
/// @param cCosRatio cos(beta)/cos(alpha) for grating (alpha > beta > 0), unitless
/// @param kLineDensity grating line density, per mm
/// @param mDiffractionOrder diffraction order reflected from grating, unitless
/// @return grating incidence angle, in degrees
double ruben2005eqn8m(
	double energy,
	double cCosRatio,
	double kLineDensity,
	int mDiffractionOrder);


/// @brief Calculates the optimal ratio cos(beta)/cos(alpha) for vanishing defocus in VIA monochromator
/// @param energy photon energy, in eV
/// @param kLineDensity grating line density, per mm
/// @param mDiffractionOrder diffraction order reflected from grating, unitless
/// @param raSourceDist source-to-grating distance, in mm
/// @param rbObjectDist source-to-exit-slit distance, in mm
/// @param b2Focus VLS focusing coefficient, in mm^-2
/// @return optimal grating angle ratio cos(beta)/cos(alpha) 
double ruben2005eqn9m(
	double energy,
	double kLineDensity,
	int mDiffractionOrder,
	float raSourceDist,
	float rbObjectDist,
	float b2Focus);


#endif
