#ifndef POINTING_H
#define POINTING_H

#include <Eigen\Dense>
#include "site.h"

using namespace Eigen;

struct pointingTerms_t {
	double IA = 10.0;
	double IB = 10.0;
	double VD = 0.0000;
	double CA = 0.0;
	double NP = 0.0000;
	double AW = 0.001;
	double AN = 0.0005;
};


Vector3d CalculateTelescopeVector(const Vector2d &AzEl, pointingTerms_t *Pointing_terms, const  Matrix3d &M);
Vector3d CalculateBoresight(const Vector3d &telescope_vector, const Vector2d &standard_coordinates);
Matrix3d CalculateMountAttitude(siteParamaters_t *site_parameters, pointingTerms_t *pointing_terms);
Vector3d CalculateAim(const Vector2d &AzEl, const Matrix3d &mount_attitude);
Vector4d CalculateAxisPosition(const Vector3d &aimVector, const Vector3d &boresightVector, double NP);
Vector4d CalculateAxisPositionFromObservedPos(Vector2d AzEl, pointingTerms_t pointingTerms, Matrix3d M);
Matrix3d R1(double x);
Matrix3d R2(double y);
Matrix3d R3(double z);
#endif
