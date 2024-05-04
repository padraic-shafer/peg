/*
Copyright 2024 Padraic Shafer.

Adapted from the work of Joe Dvorak, which is based on the published article
by R. Reininger et al., NIM-A 538, 760-770 (2005).
*/

#include "via.h"
#include <cmath>
#include <iostream>


double ruben2005eqn8m(double energy, double cCosRatio, double kLineDensity, int mDiffractionOrder) {

	// Returns the grating incidence angle for the specified diffraction order
	// per Eqn. 8 in R. Reininger et al., NIM-A 538, 760-770 (2005)
	// and generalized to arbitrary diffraction order

	// energy: photon energy, in eV
	// cCosRatio: cos(beta)/cos(alpha) for grating (alpha > beta > 0), unitless
	// kLineDensity: grating line density, per mm
	// mDiffractionOrder: diffraction order reflected from grating, unitless

	double scaleFactor = 1.0 / (cCosRatio*cCosRatio - 1.0);
	double mmWavelength = (eVnm_PRODUCT / energy) * mm_PER_nm;
	double incidenceAngleSin = (-mDiffractionOrder * kLineDensity * mmWavelength * scaleFactor) + sqrt(1.0 + std::pow(cCosRatio * mDiffractionOrder * kLineDensity * mmWavelength * scaleFactor, 2.0));
	double incidenceAngleDeg = asin(incidenceAngleSin) * 180.0 / M_PI;

	std::cout << "ruben2005eqn8m" << std::endl;
	std::cout << "  energy: " << energy << std::endl;
	std::cout << "  cCosRatio: " << cCosRatio << std::endl;
	std::cout << "  kLineDensity: " << kLineDensity << std::endl;
	std::cout << "  mDiffractionOrder: " << mDiffractionOrder << std::endl;
	std::cout << "  ---" << std::endl;
	std::cout << "  scaleFactor: " << scaleFactor << std::endl;
	std::cout << "  mmWavelength: " << mmWavelength << std::endl;
	std::cout << "  incidenceAngleSin: " << incidenceAngleSin << std::endl;
	std::cout << "  incidenceAngleDeg: " << incidenceAngleDeg << std::endl;
	std::cout << "-------------------" << std::endl;

	return incidenceAngleDeg;
}


double ruben2005eqn9m(double energy, double kLineDensity, int mDiffractionOrder, float raSourceDist, float rbObjectDist, float b2Focus) {

	// Returns the optimal ratio cos(beta)/cos(alpha) for vanishing defocus in VIA monochromator
	// per Eqn. 9 in R. Reininger et al., NIM-A 538, 760-770 (2005)
	// and generalized to arbitrary diffraction order

	// energy: photon energy, in eV
	// kLineDensity: grating line density, per mm
	// mDiffractionOrder: diffraction order reflected from grating, unitless
	// raSourceDist: source-to-grating distance, in mm
	// rbObjectDist: source-to-exit-slit distance, in mm
	// b2Focus: VLS focusing coefficient, in mm^-2

	double mmWavelength = (eVnm_PRODUCT / energy) * mm_PER_nm;
	double A0 = mDiffractionOrder * kLineDensity * mmWavelength;
	double A2 = -0.5 * mDiffractionOrder * kLineDensity * rbObjectDist * b2Focus;
	double rA2A0 = A2 / A0;
	double r = rbObjectDist / raSourceDist;
	double r1 = r + 1.0;

	double termA = 2*A2 + 4*rA2A0*rA2A0 + (4 + 2*A2 - A0*A0)*r;
	double termB = -4*rA2A0 * sqrt(r1*r1 + 2*A2*r1 - A0*A0*r);
	double termC = -4 + A0*A0 - 4*A2 + 4*rA2A0*rA2A0;
	double cCosRatio = sqrt((termA + termB) / termC);

	std::cout << "ruben2005eqn9m" << std::endl;
	std::cout << "  energy: " << energy << std::endl;
	std::cout << "  kLineDensity: " << kLineDensity << std::endl;
	std::cout << "  mDiffractionOrder: " << mDiffractionOrder << std::endl;
	std::cout << "  raSourceDist: " << raSourceDist << std::endl;
	std::cout << "  rbObjectDist: " << rbObjectDist << std::endl;
	std::cout << "  b2Focus: " << b2Focus << std::endl;
	std::cout << "  ---" << std::endl;
	std::cout << "  mmWavelength: " << mmWavelength << std::endl;
	std::cout << "  A0: " << A0 << std::endl;
	std::cout << "  A2: " << A2 << std::endl;
	std::cout << "  rA2A0: " << rA2A0 << std::endl;
	std::cout << "  r: " << r << std::endl;
	std::cout << "  r1: " << r1 << std::endl;
	std::cout << "  ---" << std::endl;
	std::cout << "  termA: " << termA << std::endl;
	std::cout << "  termB: " << termB << std::endl;
	std::cout << "  termC: " << termC << std::endl;
	std::cout << "  cCosRatio: " << cCosRatio << std::endl;
	std::cout << "-------------------" << std::endl;

	return cCosRatio;
}
