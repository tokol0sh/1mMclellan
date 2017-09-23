#include <sofa.h>
#include <stdio.h>
#include <time.h>
#include <ctime>
#include <sofa.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <math.h>
#include <Eigen\Dense>

using namespace std;
using namespace Eigen;
using namespace std::chrono;

void reprd(char* s, double ra, double dc, double hr = 0, double dr = 0, double LST = 0);
struct pointingTerms_t {
	double IA = 0.0;
	double IB = 0.0;
	double VD = 0.0;
	double CA = 0.0;
	double NP = 0.0;
	double AW = 0.0;
	double AN = 0.0;
}pointingTerms;

struct siteParamaters_t {
	double longitude = 0.0;
	double latitude = 0.0;
	double height = 0.0;
}siteParamaters;

Vector3d CalculateTelescopeVector(const Vector2d &AzEl, pointingTerms_t Pointing_terms, const  Matrix3d &M);
Vector3d CalculateBoresight(const Vector3d &telescope_vector, const Vector2d &standard_coordinates);
Vector4d CalculateAxisPosition(const Vector3d &aimVector, const Vector3d &boresightVector, double NP);
Matrix3d R1(double x);
Matrix3d R2(double y);
Matrix3d R3(double z);



void domeCalc() {
	double phi = 0;  // elevation of the north end of the polar axis
	double r_D = 3200; // radius of dome (mm)
	double p = 0; // separation of polar axis and declination axis
	double q = 0; // separation between telescope and declination axis
	double r = 0.0; // separation between telescope and declination axis
	double HA, DEC = 0.0; // Coordinates of telescope (mount coordinates)

	double x, y, z;
	double x_mount, y_mount, z_mount;
	double x_mount_offset, y_mount_offset, z_mount_offset;
	double x_dome, y_dome, z_dome;
	double x_s, y_s, z_s;
	double s_dt, t2_m, w, f;
	double azimuth_dome, elevation_dome;

	y = p + r * sin(DEC);
	x_mount = q * cos(HA) + y * sin(HA);
	y_mount = -q * sin(HA) + y * sin(HA);
	z_mount = r * cos (DEC);

	x_dome = x_mount_offset + x_mount;
	y_dome = y_mount_offset + y_mount * sin(phi) + z_mount_offset * cos(phi);
	z_dome = z_mount_offset - y_mount * cos(phi) + z_mount_offset * sin(phi);

	x = -sin(HA) * cos(DEC);
	y = -cos(HA) * cos(DEC);
	z = sin(DEC);
	x_s = x;
	y_s = y * sin(phi) + z * cos(phi);
	z_s = -y * cos(phi) + z * sin(phi);

	s_dt = (x_s * x_dome) + (y_s - y_dome) + (z_s * z_dome);
	t2_m = (x_dome * x_dome) + (y_dome * y_dome) + (z_dome * z_dome);
	w = (s_dt * s_dt) - t2_m + (r_D*r_D);
	f = -s_dt + sqrt(w);


	x = x_dome + f * x_s;
	y = y_dome + f * y_s;
	z = z_dome + f * z_s;

	azimuth_dome = atan2(x, y);
	elevation_dome = atan2(z, sqrt(x*x + y*y));
	
}





int main()
{


	typedef duration<int, ratio_multiply<hours::period, ratio<24> >::type> days;

	std::chrono::time_point<std::chrono::high_resolution_clock> previous_time_point_high_res;
	std::chrono::time_point<std::chrono::high_resolution_clock> current_time_point_high_res;
	system_clock::time_point current_time_point;
	system_clock::time_point current_time_ms;
	tm utc_tm;
	time_t current_time_t;
	milliseconds msec = duration_cast<milliseconds>(current_time_point - current_time_ms);
	double delta_t;

	double previous_RA, previous_DEC, previous_HA, LST;
	double delta_RA, delta_DEC, delta_HA;
	double rate_RA, rate_DEC, rate_HA;

	iauASTROM astrom;
	double latitude, longitude, height, barometric_pressure, temperature, relative_humidity, wavelength, utc1, utc2, tai1, tai2,
		tt1, tt2, ut11,ut12, xp, yp, dut1, dx, dy, rc, dc, pr, pd, px, rv,
		eo, ri, di, rca, dca, ra, da, aot, zot, hot, dot, rot,
		aob, zob, hob, dob, rob, aob_next, zob_next, hob_next, dob_next, rob_next,
		pvh[2][3], x, y, s;

	/* Minimum cos(alt) and sin(alt) for refraction purposes */
	const double CELMIN = 1e-6;
	const double SELMIN = 0.05;



	current_time_point = system_clock::now();
	current_time_t = system_clock::to_time_t(current_time_point);
	utc_tm = *gmtime(&current_time_t);
	current_time_ms = system_clock::from_time_t(current_time_t);



	/* Site longitude, latitude (radians) and height above the geoid (m). */
	iauAf2a('-', 43,33, 11.5, &latitude);
	iauAf2a('+', 172, 32, 41.575, &longitude );
	height = 25.0;

	/* Ambient pressure (HPa), temperature (C) and rel. humidity (frac). */
	barometric_pressure = 952.0;
	temperature = 18.5;
	relative_humidity = 0.83;

	/* Effective color (microns). */
	wavelength = 0.55;


	// These can be obtained form IERS Bulletin B: ftp://hpiers.obspm.fr/iers/bul/bulb_new/bulletinb.dat
	/* EOPs:  polar motion in radians, UT1-UTC in seconds. */
	xp = 209.842e-3 * DAS2R;
	yp = 409.914e-3 * DAS2R;
	dut1 = 346.5856e-3;

	/* Corrections to IAU 2000A CIP (radians). */
	dx = 0 * DAS2R;
	dy = 0 * DAS2R;


	/* Time  */
	iauDtf2d("UTC", (2017), utc_tm.tm_mon + 1, utc_tm.tm_mday, utc_tm.tm_hour, utc_tm.tm_min, (utc_tm.tm_sec + msec.count()/1000.0), &utc1, &utc2);
	iauUtctai(utc1, utc2, &tai1, &tai2);
	iauTaitt(tai1, tai2, &tt1, &tt2);
	iauUtcut1(utc1, utc2, dut1, &ut11, &ut12);




	/* Star ICRS RA,Dec (radians). */
	iauTf2a(' ', 4, 17, 53.72, &rc);
	iauAf2a('-', 33, 47,53.9, &dc);
	reprd("ICRS, epoch J2000.0:", rc, dc);




	/* Proper motion: RA/Dec derivatives, epoch J2000.0. */
	pr = atan2(-0 * DAS2R, cos(dc));
	pd = 0 * DAS2R;

	/* Parallax (arcsec) and recession speed (km/s). */
	px = 0.0;
	rv = 0.0;


	// Calculate star independent astrometry paramaters
	iauApco13(utc1, utc2, dut1, longitude, latitude, relative_humidity, xp, yp, barometric_pressure, temperature, relative_humidity, wavelength, &astrom, &eo);

	// Transform from ICRS to CIRS
	iauAtciq(rc, dc, pr, pd, px, rv, &astrom, &ri, &di);
	reprd("ICRS -> CIRS:", ri, di);

	/* CIRS to Apparent place. */
	ra = iauAnp(ri - eo);
	da = di;
	reprd("Classical Apparent place:", ra, da);

	/* CIRS to topocentric. */
	if (iauAtio13(ri, di, utc1, utc2, dut1, longitude, latitude, height, xp, yp, 0.0, 0.0, 0.0, 0.0,
		&aot, &zot, &hot, &dot, &rot)) return -1;
	reprd("CIRS -> topocentric:", astrom.eral - rot, dot);

	// CIRS to Observed
	iauAtioq(ri, di, &astrom, &aob, &zob, &hob, &dob, &rob);
	reprd("CIRS -> observed:", astrom.eral - rob, dob);




	///UTC date. 
	while (1) {
		current_time_point = system_clock::now();

		delta_t = 50.0;

		current_time_t = system_clock::to_time_t(current_time_point);
		utc_tm = *gmtime(&current_time_t);
		current_time_ms = system_clock::from_time_t(current_time_t);
		milliseconds msec = duration_cast<milliseconds>(current_time_point - current_time_ms);


		iauDtf2d("UTC", 1900 + utc_tm.tm_year, utc_tm.tm_mon + 1, utc_tm.tm_mday, utc_tm.tm_hour, utc_tm.tm_min, (utc_tm.tm_sec + msec.count() / 1000.0), &utc1, &utc2);
		iauUtctai(utc1, utc2, &tai1, &tai2);
		iauTaitt(tai1, tai2, &tt1, &tt2);

		// Calculate star independent astrometry paramaters
		iauApco13(utc1, utc2, dut1, longitude, latitude, relative_humidity, xp, yp, barometric_pressure, temperature, relative_humidity, wavelength, &astrom, &eo);
		// CIRS to Observed
		iauAtioq(ri, di, &astrom, &aob, &zob, &hob, &dob, &rob);

		LST = eo - astrom.eral;

		iauDtf2d("UTC", 1900 + utc_tm.tm_year, utc_tm.tm_mon + 1, utc_tm.tm_mday, utc_tm.tm_hour, utc_tm.tm_min, (delta_t + utc_tm.tm_sec + msec.count() / 1000.0), &utc1, &utc2);
		iauUtctai(utc1, utc2, &tai1, &tai2);
		iauTaitt(tai1, tai2, &tt1, &tt2);

		// Calculate star independent astrometry paramaters
		iauApco13(utc1, utc2, dut1, longitude, latitude, relative_humidity, xp, yp, barometric_pressure, temperature, relative_humidity, wavelength, &astrom, &eo);
		// CIRS to Observed
		iauAtioq(ri, di, &astrom, &aob_next, &zob_next, &hob_next, &dob_next, &rob_next);

		rate_RA = (rob_next - rob) / delta_t;
		rate_DEC = (dob_next - dob) / delta_t;
		rate_HA = (hob_next - hob) / delta_t;


		//reprd("ICRS -> observed:", rob, dob, rate_RA, rate_HA, LST);
		std::this_thread::sleep_for(std::chrono::milliseconds(1));

		Matrix3d M_0 = R1(siteParamaters.latitude - 90.0) * R2(0.0) * R3(0.0);
		Matrix3d M = R1(pointingTerms.AW) * R2(-pointingTerms.AN) * M_0;
		Vector2d AzEl(3.14159, 3.14159);
		Vector3d result;
		result = CalculateTelescopeVector(AzEl, pointingTerms, M);
		cout << result << endl;
	}


	return 0;
}



Vector3d CalculateTelescopeVector(const Vector2d &AzEl, pointingTerms_t Pointing_terms, const Matrix3d &M) {
	Vector2d Vertical_displacement(0.0, Pointing_terms.VD);
	Vector2d AzEl_0 = AzEl + Vertical_displacement;
	Vector2d pointing_origin(0,0);

	Vector3d V(cos(AzEl(0)) * cos(AzEl(1)), sin(AzEl(0)) * cos(AzEl(1)), sin(AzEl(1)));
	Vector3d V_0(cos(AzEl_0(0)) * cos(AzEl_0(1)), sin(AzEl_0(0)) * cos(AzEl_0(1)), sin(AzEl_0(1)));

	Vector3d A = M * V;
	Vector3d A_0 = M * V_0;

	Vector3d telescope_0(cos(Pointing_terms.CA), sin(Pointing_terms.CA), 0);

	Vector3d undeflected_boresight = CalculateBoresight(telescope_0, pointing_origin);

	Vector4d Undeflected_telescope = CalculateAxisPosition(A_0, undeflected_boresight, Pointing_terms.NP);
	Vector4d deflected_telescope = CalculateAxisPosition(A, undeflected_boresight, Pointing_terms.NP);

	return Vector3d(Undeflected_telescope(0), Undeflected_telescope(1), 0);
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

Vector4d CalculateAxisPosition(const Vector3d &aimVector, const Vector3d &boresightVector, double NP) {
	double pitch1 = atan2( (aimVector(2) + sin(NP) * boresightVector(1)), + sqrt(aimVector(0) * aimVector(0) + aimVector(1) * aimVector(1) - boresightVector(1)*(2 * aimVector(2) * sin(NP) + boresightVector(1)) - sin(NP) * sin(NP))) - atan2(boresightVector(2), boresightVector(0));
	double pitch2 = atan2( (aimVector(2) + sin(NP) * boresightVector(1)), - sqrt(aimVector(0) * aimVector(0) + aimVector(1) * aimVector(1) - boresightVector(1)*(2 * aimVector(2) * sin(NP) + boresightVector(1)) - sin(NP) * sin(NP))) - atan2(boresightVector(2), boresightVector(0));
	
	Vector3d NProtation1 = R1(NP) * R2(pitch1) * boresightVector;
	Vector3d NProtation2 = R1(NP) * R2(pitch2) * boresightVector;

	double roll1 = atan2((aimVector(1) * NProtation1(0) - aimVector(0) * NProtation1(1)), (aimVector(0) * NProtation1(0) + aimVector(1) * NProtation1(1)));
	double roll2 = atan2((aimVector(1) * NProtation2(0) - aimVector(0) * NProtation2(1)), (aimVector(0) * NProtation2(0) + aimVector(1) * NProtation2(1)));

	Vector4d results(pitch1, roll1, pitch2, roll2);
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


void reprd(char* s, double ra, double dc, double hr , double dr , double LST)
{
	char pm = '-';
	int i[4];
	printf("%25s", s);
	iauA2tf(7, ra, &pm, i);
	printf(" %2.2d %2.2d %2.2d.%7.7d", i[0], i[1], i[2], i[3]);
	iauA2af(7, dc, &pm, i);
	printf(" %2.2d %2.2d %2.2d.%7.7d", i[0], i[1], i[2], i[3]);
	iauA2af(7, hr, &pm, i);
	printf(" %2.2d %2.2d %2.2d.%7.7d", i[0], i[1], i[2], i[3]);
	iauA2af(7, dr, &pm, i);
	printf(" %2.2d %2.2d %2.2d.%7.7d", i[0], i[1], i[2], i[3]);
	iauA2tf(7, LST, &pm, i);
	printf(" %2.2d %2.2d %2.2d.%7.7d\n", i[0], i[1], i[2], i[3]);
}


