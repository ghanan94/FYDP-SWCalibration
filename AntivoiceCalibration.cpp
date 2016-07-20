//
//Definitions: 
	//Strategy document: https://docs.google.com/document/d/10Bu1yvGFMuSSRcf0y2kXYmD8a38s4vKvYzKdGpoMt9Q/edit
	//MOCK: Simplifies a hardware/embedded problem to the best that can be done to model it using command line sw.
	//STUB: can't do anything in the command line sw to actually represent this action.
/*Pending TODOs:
	-After applying gaussian noise, amplitude detection is not as accurate
		-use BPF to improve accuracy:
			https://sestevenson.wordpress.com/implementation-of-fir-filtering-in-c-part-1/
	
	-do not need to store test_signal as you can calculate each index
	-find a way to perform obs signal analysis without storing the entire signal in memory (lower memory consumption)
*/
	
#define WINDOWS 1

#ifdef WINDOWS
	#include "stdafx.h"
#endif

#include <math.h>
#include <stdlib.h>

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
	double internaltime; //rtt minus airtime (does not include computation time)
};
typedef struct calibration CALIB_t;

struct amplitude_delay {
	double amplitude; //amplitude of signal
	double delay; //delay in seconds it took for round trip time
};
typedef struct amplitude_delay AMPDELAY_t;



/* Global variables */
const unsigned int freq_count_pos = 4;
double frequencies[freq_count_pos] = { 1, 2, 3, 4 };
double system_gain[freq_count_pos] = { 0.9, 0.8, 0.7, 0.6 };
double system_delay[freq_count_pos] = { 0.1,0.2,0.3,0.4 };

double channel_delay_seconds = (double)TRANSDUCER_OBS_SEPARATION_DISTANCE_METERS / (double)SPEED_OF_SOUND;

const unsigned int amplitude_count = 3;
double amplitudes[amplitude_count] = { 1, 0.8, 0.6 };
/* Global variables end */

/*
	Helper function for getting frequency index from Hertz
	params:
		f is the frequency in Hertz
	returns:
		i, index into frequencies array
*/
unsigned int GetFrequenciesArrayIndex(double f)
{
	unsigned int i = 0; //i is the index
	while (i < freq_count_pos)
	{
		//incase there is a small round error see if approximately the same
		double abs_difference = MAX(f,frequencies[i]) - MIN(f,frequencies[i]);
		
		if (f == frequencies[i] || abs_difference < 0.01)
		{
			break;
		}
		i++;
	}
	return i;
}

/*
	params: 
		f: frequency
		duration: in seconds
		fs: sample rate for DAC, thus signal resolution
		ampltiude: amplitude of sinusoid
		signal_array_length: length of test signal to be generated
	returns:
		an array of points where each point represents the time domain height of the signal
*/
double * TestSignalGenerator(double f, double duration, unsigned int fs, double amplitude, unsigned int signal_array_length)
{
	double * test_signal = 0;
	test_signal = (double *)malloc(signal_array_length * sizeof(double)); //allocate memory

	unsigned int i; //the sample number
	for (i = 0; i < signal_array_length; i++)
	{
		//create the test signal point for each sample i 
		test_signal[i] = amplitude * cos(2 * PI * f * (i / (double) fs));
	}

	return test_signal;
}

/*
	Helper function
	Apply @num_samples delay to @signal in-place
*/
void ApplySampleDelay(double * signal, unsigned int num_samples, unsigned int signal_array_length)
{
	if (num_samples < signal_array_length)
	{
		//Apply specified delay
		unsigned int delay_bytes = (signal_array_length - num_samples) * sizeof(double);
		memcpy(&signal[num_samples], signal, delay_bytes);

		//zero out the delayed amount
		memset(signal, 0, num_samples * sizeof(double));
	}
	else
	{
		//zero out the entire signal and show warning
		memset(signal, 0, signal_array_length * sizeof(double));
		printf("\n ****** \n [WARNING] delay is greater than signal length (duration) \n ****** \n");
	}
}

/*
	Applies magnitude and group delay transfer function to test signal at frequency f
	parameters:
		f: frequency in Hz
		test_signal: array of points, withe each point representing signal height at that sample number
	returns:
		changed_signal, with transfer function applied
*/
double * ApplyMicTransducerTransferFunction(double f, unsigned int fs, double * test_signal, unsigned int signal_array_length)
{
	unsigned int i;

	double * changed_signal = (double *)malloc(signal_array_length * sizeof(double));
	memcpy(changed_signal, test_signal, signal_array_length * sizeof(double)); //perform deep copy of test_signal to changed_signal

	//Get gain and delay amounts for this frequency f
	double gain_amount = system_gain[GetFrequenciesArrayIndex(f)];
	double seconds_delay = system_delay[GetFrequenciesArrayIndex(f)];
	unsigned int samples_delay = (unsigned int)round(seconds_delay * fs);

	for (i = 0; i < signal_array_length; i++)
	{
		//Apply gain
		changed_signal[i] *= gain_amount;
	}

	ApplySampleDelay(changed_signal, samples_delay, signal_array_length);
	
	return changed_signal;
}

/* 
	Generates additive white Gaussian Noise samples with zero mean and a standard deviation of @d. 
	This method is invoked once per sample, and the value it returns should be added to the signal.
	Author: https://www.embeddedrelated.com/showcode/311.php
*/
double AWGN_generator(double d)
{
	double temp1;
	double temp2;
	double result;
	int p;

	p = 1;

	while (p > 0)
	{
		//rand() function generates an integer between 0 and RAND_MAX,
			//which is defined in stdlib.h.
		temp2 = (rand() / ((double)RAND_MAX));

		if (temp2 == 0)
		{// temp2 is >= (RAND_MAX / 2)
			p = 1;
		}// end if
		else
		{// temp2 is < (RAND_MAX / 2)
			p = -1;
		}// end else

	}// end while()

	temp1 = cos((2.0 * (double)PI) * rand() / ((double)RAND_MAX));
	result = sqrt(-2.0 * log(temp2)) * temp1;

	return (result * d);
}

/*
	Applies channel noise to signal and air time
parameters:
	in_signal: array of points, withe each point representing signal height at that sample number
	ch_delay: time in seconds the channel delays by
returns:
	out_signal, with channel noise and delay applied
*/
double * ApplyChannel(double * in_signal, unsigned int samples_delay, unsigned int fs, unsigned int signal_array_length)
{
	unsigned int i;

	double * out_signal = in_signal; //no deep copy required since in_signal is intermediate and can be readily modified

	//Apply channel delay
	ApplySampleDelay(out_signal, samples_delay, signal_array_length);

	//Apply channel noise 
	for (i = 0; i < signal_array_length; i++)
	{
		out_signal[i] += AWGN_generator(WHITE_NOISE_STD_DEVIATION);
	}

	return out_signal;
}

/*
	Returns the amplitude of ideal cosine signal.  
		Averages all the minimums and maximums, adds them, and divide by 2.
		This algorithm underestimates the amplitude a little, improves with higher fs.
		//Optimization available possibly:
		//We know h is of the form: h[k] = (a)cos(2 * PI * f * (k / N)
		//therefore amplitude a = h[k] / cos(2 * PI * f * (k / N)
		//can also assume sinusoid is linear near zero-crossings, and therefore fit the best amplitude between 2 zeros

	params:
		h: signal level at each sample
		f: the frequency in Hz
		fs: sample frequency in Hz
	returns:
		{amplitude of sinusoid, signal delay}
*/
AMPDELAY_t DetectSinusoidalAmplitudeAndDelay(double * h, double f, unsigned int fs, unsigned int signal_array_length)
{
	AMPDELAY_t ret;
	double amplitude;

	int N = (unsigned int)ceil((double)fs / f); //Samples in period

	//Step 1: Identify max and min levels across entire signal, and their indices
	int i = 0;
	double min = 0, max = 0;
	int min_index = 0, max_index = 0;
	while(i < signal_array_length)
	{
		if (h[i] <= min)
		{
			min = h[i];
			min_index = i;
		}

		if (h[i] >= max)
		{
			max = h[i];
			max_index = i;
		}
		i++;
	}

	//Step 2: Sum max and min levels for all periods
	double minsum = 0, maxsum = 0;
	int mincount = 0, maxcount = 0;

	//min levels
	i = min_index;
	while (i < signal_array_length)
	{
		if (h[i] < 0.5 * min) //ensure not zeroed out
		{
			minsum += h[i];
			mincount += 1;
		}
		
		i += N;
	}

	if (min_index - N > 0)
	{
		i = min_index - N;
		while (i >= 0)
		{
			if (h[i] < 0.5 * min) //ensure not zeroed out
			{
				minsum += h[i];
				mincount += 1;
			}
			i -= N;
		}
	}

	//max levels
	i = max_index;
	while (i < signal_array_length)
	{
		if (h[i] > 0.5 * max) //ensure not zeroed out
		{
			maxsum += h[i];
			maxcount += 1;
		}
		i += N;
	}

	if (max_index - N > 0)
	{
		i = max_index - N;
		while (i >= 0)
		{
			if (h[i] > 0.5 * max) //ensure not zeroed out
			{
				maxsum += h[i];
				maxcount += 1;
			}
			i -= N;
		}
	}

	//Step 3: Find average of max and min
	double minavg = minsum / mincount;
	double maxavg = maxsum / maxcount;

	//Step 4: Calculate amplitude
	amplitude = (maxavg - minavg) / 2;
	ret.amplitude = amplitude;

	//Step 5: Calculate delay
	for (i = 0; i < signal_array_length; i++)
	{
		if (h[i] > 0.5 * max)
		{
			break; //at this i, cosine maxima starts therefore i is the number of samples delayed by
		}
	}
	unsigned int delay_samples = (unsigned int) i;
	double delay_seconds = (double)delay_samples / (double)fs;
	ret.delay = delay_seconds;
	
	return ret;
}

/*
	Performs calibration routine
	parameters:
		f: the frequency at which to calibrate in Hertz
		test_amplitude: the amplitude of the signal being generated
	returns:
		CALIB_t containing calibration results
*/
CALIB_t GetCalibration(double f, double test_amplitude, unsigned int fs, unsigned int channel_delay_samples)
{
	//Generate test signal
	double duration = 1; //seconds

	unsigned int signal_array_length = (unsigned int)ceil(duration * fs); //total number of samples

	double * test_signal = TestSignalGenerator(f, duration, fs, test_amplitude, signal_array_length);

	//
	//STUB: Use DAC to send to Speaker driver
	//

	//MOCK: We model the time it takes for sound to travel from transducer to observation microphone
		//v = d/t => t = d/v;  
	
	//
	//STUB: ADC of observation microphone to MCU
	//

	//MOCKS (transducer and mic transfer functions unified for this simulation):
		//MOCK: transducer magnitude response and group delay
		//MOCK: observation magnitude response and group delay
		//transducer and obs mic are bundled together as 1 transfer function
	
	double * changed_signal = ApplyMicTransducerTransferFunction(f, fs, test_signal, signal_array_length);


	//MOCK: Channel (air / thermal) noise, model as gaussian noise
	double * obs_signal = ApplyChannel(changed_signal, channel_delay_samples, fs, signal_array_length);

	//Analyze the obs signal relative to test signal to generate calibration outputs needed for antivoice algorithm
		//filter noise from obs signal by using BPF centered at f
	//TODO above^^
	
	//determine amplitude and seconds delay of signal relative to test
	AMPDELAY_t obs_ampdelay = DetectSinusoidalAmplitudeAndDelay(obs_signal, f, fs, signal_array_length);
	
	//Determine results
	CALIB_t result;
	result.frequency = f;
	result.amplitude = test_amplitude;
	result.eqgain = test_amplitude / obs_ampdelay.amplitude;
	result.rtt = obs_ampdelay.delay;
	result.internaltime = result.rtt - channel_delay_seconds;

	//free memory
	free(test_signal);
	free(obs_signal);

	return result;
}		

void PrintCalib_t(CALIB_t c)
{
	printf("%.2f | %.3f      | %.5f | %.5f | %.5f | %.9f | %.9f |\n", 
		c.frequency, c.amplitude, 
		system_gain[GetFrequenciesArrayIndex(c.frequency)], system_delay[GetFrequenciesArrayIndex(c.frequency)],
		c.eqgain, c.rtt, c.internaltime);
}

int main(int argc, char ** argv)
{
	//TODO: make frequencies be passed as an array by command line?

	//samples per second or resolution of sinusoid (highest freq * some constant)
	unsigned int fs = (unsigned int)round(frequencies[freq_count_pos - 1] * 250);

	//Calculate channel delay in terms of samples
	unsigned int channel_delay_samples = (unsigned int)round(channel_delay_seconds * fs);
	printf("\n[Warning] channel delay of %.6f seconds rounded to %.6f seconds due to nearest sample.",
		channel_delay_seconds, channel_delay_samples / (double) fs);
	channel_delay_seconds = channel_delay_samples / (double)fs;
	
	CALIB_t calib_results[freq_count_pos * amplitude_count]; //stores the calibratoin results

	printf("\nFreq | Amplitude  | Sysgain | SysDelay| eqGain  |  RTT        | Int. Delay  | \n");

	//iterate over all frequencies and collect results for each frequency
	unsigned int i,j;
	for (i = 0; i < freq_count_pos; i++)
	{
		for (j = 0; j < amplitude_count; j++)
		{
			calib_results[i + j] = GetCalibration(frequencies[i], amplitudes[j], fs, channel_delay_samples);
			PrintCalib_t(calib_results[i + j]);
		}
	}

    return 0;
}

