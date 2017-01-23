#ifndef CALIBRATION_SETUP_H
#define CALIBRATION_SETUP_H


#define WINDOWS 1

#ifdef WINDOWS
#include "stdafx.h"
#endif

#include <math.h>
#include <stdlib.h>

#ifndef PI
# define PI	3.14159265358979323846264338327950288
#endif

#define SIM 1
#define BOARD (!SIM)

#define SPEED_OF_SOUND (343)
#define TRANSDUCER_OBS_SEPARATION_DISTANCE_METERS (1.3)

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

struct calibration {
	double frequency; //frequency of tone
	double amplitude; //amplitude of tone
	double eqgain; //gain needed to equalize (scalar factor, not dB)
	double rtt; //delay in seconds it took for round trip time
	double internaltime; //rtt minus airtime (signal propagation + computation time)
};
typedef struct calibration CALIB_t;

struct amplitude_delay {
	double amplitude; //amplitude of signal
	double delay; //delay in seconds it took for round trip time
};
typedef struct amplitude_delay AMPDELAY_t;


#define TEST_DUR (0.7)	//seconds for a given test signal
#define MIC_FS (8000)		//sinusoidal resolution, as well as sampling rate
#define DAC_FS (8000)

/* Global variables */

const unsigned int signal_array_length = (unsigned int) (TEST_DUR * MIC_FS); //total number of samples

const unsigned int freq_count_pos = 4;
double frequencies[freq_count_pos] = { 10, 50, 100, 200 };


const unsigned int amplitude_count = 3;
double amplitudes[amplitude_count] = { 1, 0.8, 0.6 };
/* Global variables end */

#endif //CALIBRATION_SETUP_H