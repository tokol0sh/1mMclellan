#include <sofa.h>
#include <stdio.h>
#include <time.h>
#include <ctime>
#include <sofa.h>
#include <chrono>
#include <iostream>

void reprd(char*, double, double);



int main()
{
	using namespace std;
	using namespace std::chrono;

	typedef duration<int, ratio_multiply<hours::period, ratio<24> >::type> days;

	system_clock::time_point current_time_point;
	system_clock::time_point current_time_ms;
	tm utc_tm;
	time_t current_time_t;

	iauASTROM astrom;
	iauLDBODY b[3];
	double latitude, longitude, height, barometric_pressure, temperature, relative_humidity, wavelength, utc1, utc2, tai1, tai2,
		tt1, tt2, ut11,ut12, xp, yp, dut1, dx, dy, rc, dc, pr, pd, px, rv,
		eo, ri, di, rca, dca, ra, da, aot, zot, hot, dot, rot,
		aob, zob, hob, dob, rob,
		pvh[2][3], pvb[2][3], r[3][3], x, y, s;


	current_time_point = system_clock::now();
	current_time_t = system_clock::to_time_t(current_time_point);
	utc_tm = *gmtime(&current_time_t);
	current_time_ms = system_clock::from_time_t(current_time_t);
	milliseconds msec = duration_cast<milliseconds>(current_time_point - current_time_ms);


	/* Site longitude, latitude (radians) and height above the geoid (m). */
	iauAf2a('-', 0,0, 0, &longitude);
	iauAf2a('-', 0, 0, 0, &latitude);
	height = 625.0;

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
	iauTf2a(' ', 22, 8, 13.98473, &rc);
	iauAf2a('-', 46, 57, 39.5078, &dc);
	reprd("ICRS, epoch J2000.0:", rc, dc);

	/* Proper motion: RA/Dec derivatives, epoch J2000.0. */
	pr = atan2(-0 * DAS2R, cos(dc));
	pd = 0 * DAS2R;

	/* Parallax (arcsec) and recession speed (km/s). */
	px = 0.0;
	rv = 0.0;

	/* ICRS to CIRS (geocentric observer). */
	iauAtci13(rc, dc, pr, pd, px, rv, tt1, tt2, &ri, &di, &eo);
	reprd("catalog -> CIRS:", ri, di);

	/* CIRS to ICRS (astrometric). */
	iauAtic13(ri, di, tt1, tt2, &rca, &dca, &eo);
	reprd("CIRS -> astrometric:", rca, dca);

	/* ICRS (astrometric) to CIRS (geocentric observer). */
	iauAtci13(rca, dca, 0.0, 0.0, 0.0, 0.0, tt1, tt2, &ri, &di, &eo);
	reprd("astrometric -> CIRS:", ri, di);

	/* Apparent place. */
	ra = iauAnp(ri - eo);
	da = di;
	reprd("geocentric apparent:", ra, da);

	/* CIRS to topocentric. */
	if (iauAtio13(ri, di, utc1, utc2, dut1, longitude, latitude, height, xp, yp,
		0.0, 0.0, 0.0, 0.0,
		&aot, &zot, &hot, &dot, &rot)) return -1;
	reprd("CIRS -> topocentric:", rot, dot);

	/* CIRS to observed. */
	if (iauAtio13(ri, di, utc1, utc2, dut1, longitude, latitude, height, xp, yp,
		barometric_pressure, temperature, relative_humidity, wavelength,
		&aob, &zob, &hob, &dob, &rob)) return -1;
	reprd("CIRS -> observed:", rob, dob);

	/* ICRS to observed. */


	/* UTC date. */
	while (1) {
		current_time_point = system_clock::now();
		current_time_t = system_clock::to_time_t(current_time_point);
		utc_tm = *gmtime(&current_time_t);
		current_time_ms = system_clock::from_time_t(current_time_t);
		milliseconds msec = duration_cast<milliseconds>(current_time_point - current_time_ms);

		int test = 0;
		test = (iauDtf2d("UTC", 1900 + utc_tm.tm_year, utc_tm.tm_mon+1, utc_tm.tm_mday, utc_tm.tm_hour, utc_tm.tm_min, (utc_tm.tm_sec + msec.count() / 1000.0), &utc1, &utc2));
		//std::cout << utc_tm.tm_hour << ':' << utc_tm.tm_min << ':' << (utc_tm.tm_sec + msec.count() / 1000.0) << '\n' ;
		/* TT date. */
		if (iauUtctai(utc1, utc2, &tai1, &tai2)) return -1;
		if (iauTaitt(tai1, tai2, &tt1, &tt2)) return -1;

		if (iauAtco13(rc, dc, pr, pd, px, rv, utc1, utc2, dut1,
			longitude, latitude, height, xp, yp, barometric_pressure, temperature, relative_humidity, wavelength,
			&aob, &zob, &hob, &dob, &rob, &eo)) return -1;
		reprd("ICRS -> observed:", rob, dob);

	}
	
	return 0;
}

void reprd(char* s, double ra, double dc)
{
	char pm;
	int i[4];
	printf("%25s", s);
	iauA2tf(7, ra, &pm, i);
	printf(" %2.2d %2.2d %2.2d.%7.7d", i[0], i[1], i[2], i[3]);
	iauA2af(6, dc, &pm, i);
	printf(" %c%2.2d %2.2d %2.2d.%6.6d\n", pm, i[0], i[1], i[2], i[3]);
}
