#ifndef PHYSICAL_MOCKING_H
#define PHYSICAL_MOCKING_H


/*
	This file contains all functions which work only in a command line program
	and not usable on the board since they mimick the physical environment.
*/

#include "calibration_setup.h"

#define WHITE_NOISE_STD_DEVIATION (0.05)

double system_gain[freq_count_pos] = { 0.9, 0.8, 0.7, 0.6 };
double system_delay[freq_count_pos] = { 0.1,0.2,0.3,0.4 };
double channel_delay_seconds = (double)TRANSDUCER_OBS_SEPARATION_DISTANCE_METERS / (double)SPEED_OF_SOUND;
double generated_signal[signal_array_length];



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
		double abs_difference = MAX(f, frequencies[i]) - MIN(f, frequencies[i]);

		if (f == frequencies[i] || abs_difference < 0.01)
		{
			break;
		}
		i++;
	}
	return i;
}



/*
Helper function
Apply @num_samples delay to @signal in-place
*/
void ApplySampleDelay(double * signal, unsigned int num_samples)
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
test_signal: array of points, with each point representing signal height at that sample number
returns:
changed_signal, with transfer function applied
*/
double * ApplyMicTransducerTransferFunction(double f, unsigned int fs, double * test_signal)
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

	ApplySampleDelay(changed_signal, samples_delay);

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
double * ApplyChannel(double * in_signal, unsigned int fs)
{
	unsigned int i;

	double * out_signal = in_signal; //no deep copy required since in_signal is intermediate and can be readily modified

									 //Apply channel delay
	unsigned int channel_delay_samples = (unsigned int)round(channel_delay_seconds * fs);
	ApplySampleDelay(out_signal, channel_delay_samples);

	//Apply channel noise 
	for (i = 0; i < signal_array_length; i++)
	{
		out_signal[i] += AWGN_generator(WHITE_NOISE_STD_DEVIATION);
	}

	return out_signal;
}

void PrintCalib_t(CALIB_t c)
{
	printf("%.2f | %.3f      | %.5f | %.5f | %.5f | %.9f | %.9f |\n",
		c.frequency, c.amplitude,
		system_gain[GetFrequenciesArrayIndex(c.frequency)], system_delay[GetFrequenciesArrayIndex(c.frequency)],
		c.eqgain, c.rtt, c.internaltime);
}

#endif //PHYSICAL_MOCKING_H