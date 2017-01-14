#define WINDOWS 1

#ifdef WINDOWS
#include "stdafx.h"
#endif

#include <math.h>
#include <stdlib.h>
#include "fir_filtering.h"

#ifndef PI
# define PI	3.14159265358979323846264338327950288
#endif

#define SPEED_OF_SOUND (343)
#define TRANSDUCER_OBS_SEPARATION_DISTANCE_METERS (1.3)
#define WHITE_NOISE_STD_DEVIATION (0.05)

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



/* Global variables */
const unsigned int fs = 8000; //sinusoidal resolution, as well as sampling rate
const double duration = 1; //seconds for a given test signal

const unsigned int freq_count_pos = 4;
double frequencies[freq_count_pos] = { 10, 50, 100, 200 };
double system_gain[freq_count_pos] = { 0.9, 0.8, 0.7, 0.6 };
double system_delay[freq_count_pos] = { 0.1,0.2,0.3,0.4 };

double channel_delay_seconds = (double)TRANSDUCER_OBS_SEPARATION_DISTANCE_METERS / (double)SPEED_OF_SOUND;

const unsigned int amplitude_count = 3;
double amplitudes[amplitude_count] = { 1, 0.8, 0.6 };
/* Global variables end */
