#include "pointing.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>

using namespace std;
using namespace Eigen;



Vector3d CalculateTelescopeVector(const Vector2d &AzEl, pointingTerms_t *Pointing_terms, const Matrix3d &M) {
	Vector2d diff_direction_1;
	Vector2d diff_direction_2;
	Vector3d telescope_vector;

	Vector2d Vertical_displacement(0.0, Pointing_terms->VD);
	Vector2d AzEl_0 = AzEl + Vertical_displacement;
	Vector2d pointing_origin(0, 0);

	Vector3d V(cos(AzEl(0)) * cos(AzEl(1)), sin(AzEl(0)) * cos(AzEl(1)), sin(AzEl(1)));
	Vector3d V_0(cos(AzEl_0(0)) * cos(AzEl_0(1)), sin(AzEl_0(0)) * cos(AzEl_0(1)), sin(AzEl_0(1)));

	Vector3d A = M * V;
	Vector3d A_0 = M * V_0;

	Vector3d telescope_0(cos(Pointing_terms->CA), sin(Pointing_terms->CA), 0);

	Vector3d undeflected_boresight = CalculateBoresight(telescope_0, pointing_origin);

	Vector4d undeflected_telescope = CalculateAxisPosition(A_0, undeflected_boresight, Pointing_terms->NP);
	Vector4d deflected_telescope = CalculateAxisPosition(A, undeflected_boresight, Pointing_terms->NP);

	diff_direction_1(0) = undeflected_telescope(0) - deflected_telescope(0);
	diff_direction_1(1) = undeflected_telescope(1) - deflected_telescope(1);

	diff_direction_2(0) = undeflected_telescope(2) - deflected_telescope(2);
	diff_direction_2(1) = undeflected_telescope(3) - deflected_telescope(3);

	telescope_vector(0) = telescope_0(0) - diff_direction_1(0) * telescope_0(1) - diff_direction_1(1)*telescope_0(2) * telescope_0(0);
	telescope_vector(1) = telescope_0(1) - diff_direction_1(0) * telescope_0(0) - diff_direction_1(1)*telescope_0(2) * telescope_0(1);
	telescope_vector(2) = telescope_0(2) + diff_direction_1(1);

	telescope_vector / sqrt(1 + diff_direction_1(0) * diff_direction_1(0) + diff_direction_1(1) * diff_direction_1(1));

	return telescope_vector;
}


Vector3d CalculateBoresight(const Vector3d &telescope_vector, const Vector2d &standard_coordinates) {
	Vector3d boresight;
	double radius = sqrt(telescope_vector(0) * telescope_vector(0) + telescope_vector(1) * telescope_vector(1));

	boresight(0) = telescope_vector(0) - (standard_coordinates(0) * telescope_vector(1) + standard_coordinates(1) * telescope_vector(2) * telescope_vector(0)) / radius;
	boresight(1) = telescope_vector(1) - (standard_coordinates(0) * telescope_vector(0) + standard_coordinates(1) * telescope_vector(2) * telescope_vector(1)) / radius;
	boresight(2) = telescope_vector(2) + standard_coordinates(1) * radius;

	boresight /= sqrt(1 + standard_coordinates(0) * standard_coordinates(0) + standard_coordinates(1) * standard_coordinates(1));

	return boresight;
}


Matrix3d CalculateMountAttitude(siteParamaters_t *site_parameters, pointingTerms_t *pointing_terms) {
	Matrix3d M_0 = R1(0.0) * R2((-90 - 90)*(M_PI / 180)) * R3(0.0);
	Matrix3d M = R1(pointing_terms->AW) * R2(-pointing_terms->AN);//* M_0;
	return M;
}


Vector3d CalculateAim(const Vector2d &AzEl, const Matrix3d &mount_attitude) {
	Vector3d aim;
	Vector3d pointing_direction_unit_vector;

	pointing_direction_unit_vector(0) = cos(AzEl(0)) * cos(AzEl(1));
	pointing_direction_unit_vector(1)= sin(AzEl(0)) * cos(AzEl(1));
	pointing_direction_unit_vector(2) = sin(AzEl(1));

	aim = mount_attitude * pointing_direction_unit_vector;
	return aim;
}

Vector4d CalculateAxisPosition(const Vector3d &aimVector, const Vector3d &boresightVector, double NP) {
	double pitch1 = atan2( (aimVector(2) + sin(NP) * boresightVector(1)) , +sqrt(aimVector(0) * aimVector(0) + aimVector(1) * aimVector(1) - boresightVector(1)*(2 * aimVector(2) * sin(NP) + boresightVector(1)) - sin(NP) * sin(NP))) - atan2(boresightVector(2), boresightVector(0));
	double pitch2 = atan2((aimVector(2) + sin(NP) * boresightVector(1)) , -sqrt(aimVector(0) * aimVector(0) + aimVector(1) * aimVector(1) - boresightVector(1)*(2 * aimVector(2) * sin(NP) + boresightVector(1)) - sin(NP) * sin(NP))) - atan2(boresightVector(2), boresightVector(0));

	Vector3d NProtation1 = R1(NP) * R2(pitch1) * boresightVector;
	Vector3d NProtation2 = R1(NP) * R2(pitch2) * boresightVector;

	double roll1 = atan2((aimVector(1) * NProtation1(0) - aimVector(0) * NProtation1(1)), (aimVector(0) * NProtation1(0) + aimVector(1) * NProtation1(1)));
	double roll2 = atan2((aimVector(1) * NProtation2(0) - aimVector(0) * NProtation2(1)), (aimVector(0) * NProtation2(0) + aimVector(1) * NProtation2(1)));

	Vector4d results(roll1, pitch1, roll2, pitch2);
	return results;
}



Matrix3d R1(double x) {
	Matrix3d m;
	m(0, 0) = 1.0;
	m(0, 1) = 0.0;
	m(0, 2) = 0.0;

	m(1, 0) = 0.0;
	m(1, 1) = cos(x);
	m(1, 2) = -sin(x);

	m(2, 0) = 0.0;
	m(2, 1) = sin(x);
	m(2, 2) = cos(x);

	return (m);
}

Matrix3d R2(double y) {
	Matrix3d m;
	m(0, 0) = cos(y);
	m(0, 1) = 0.0;
	m(0, 2) = sin(y);

	m(1, 0) = 0.0;
	m(1, 1) = 1.0;
	m(1, 2) = 0.0;

	m(2, 0) = -sin(y);
	m(2, 1) = 0.0;
	m(2, 2) = cos(y);

	return (m);
}

Matrix3d R3(double z) {
	Matrix3d m;
	m(0, 0) = cos(z);
	m(0, 1) = -sin(z);
	m(0, 2) = 0.0;

	m(1, 0) = sin(z);
	m(1, 1) = cos(z);
	m(1, 2) = 0.0;

	m(2, 0) = 0.0;
	m(2, 1) = 0.0;
	m(2, 2) = 1.0;

	return (m);
}